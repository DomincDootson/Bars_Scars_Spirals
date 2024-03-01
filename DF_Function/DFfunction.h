#ifndef DFFUNCTION
#define DFFUNCTION

#include <vector>
#include <complex>
#include <string>
#include <iostream>

#include <Eigen/Dense>

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/DFClass.h"
#include "../Volterra_Solver/VolterraSolver.h"

class DFfunction
{
public:
	template <class T>
	DFfunction(const ActionAngleBasisContainer & bf, const VolterraSolver & solver, const T & df) :
	m_m1Max{bf.maxFourierHarmonic()}, m_fourierHarmonic{bf.fourierHarmonic()}, m_maxRadialIndex{bf.maxRadialIndex()}, m_dfSize{bf.size(0)},
	m_numbTimeStep{solver.numbTimeSteps()}, m_timeStep{solver.timeStep()}, m_dfStep{bf.step()},
	m_omega1{df.omega1Grid(df).array(bf.size(0), bf.spacing())}, m_omega2{df.omega2Grid(df).array(bf.size(0), bf.spacing())}, 


	m_j2apoPeriJacobian(df.energyAngMomJacobain(bf.size(0), bf.spacing())),
	m_rPeri{Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize)}, m_rApo{Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize)}, 
	m_eCoords{Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize)}, m_lCoords{Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize)},
	v_dFdJ{dFdJsetUp(bf, df)}
	{
		coordSetUp(df, bf.spacing());
	}

	~DFfunction() {}

	template <class T>
	void dfEvolution(const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver, const T & df, const double r, const double phi);
	
	void df2fileApoPeri(const std::string & filename, const int skip = 50);
	void df2fileEL(const std::string & filename, const int skip = 50);

private: 
	const int m_m1Max, m_fourierHarmonic, m_maxRadialIndex, m_dfSize, m_numbTimeStep;
	const double m_timeStep, m_dfStep;

	const Eigen::ArrayXXd m_omega1, m_omega2, m_j2apoPeriJacobian;
	Eigen::ArrayXXd m_theta1, m_theta2, m_rPeri, m_rApo, m_eCoords, m_lCoords; 
	std::vector<Eigen::ArrayXXcd> v_exponential, v_thetaExponential, v_potentialTimesdFdJ;
	
	const std::vector<Eigen::ArrayXXd> v_dFdJ;
	std::vector<Eigen::ArrayXXd>  v_DF; 

	template <class T>
	std::vector<Eigen::ArrayXXd> dFdJsetUp(const ActionAngleBasisContainer & bf, const T & df);
	Eigen::ArrayXXd dFdJGrid(const int m1) {return v_dFdJ[m1 + m_m1Max];}

	int radius2index(double r) const {return static_cast<int>(round(r/m_dfSize));} // 
	template <class T>
	Eigen::ArrayXXd theta1(const double radius, const double phi, const T & df) const; // How do we deal with radii that are not between r- and r+? 

	template <class T>
	Eigen::ArrayXXd theta2(const double radius, const double phi, const T & df) const;

	template <class T>
	void getTheta(const double radius, const double phi, const T & df) {m_theta1 = theta1(radius, phi, df); m_theta2 = theta2(radius, phi, df);}

	template <class T >
	void coordSetUp(const T & df, const double spacing);


	void potentialTerm(const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver, const int timeIndex); 
	Eigen::ArrayXXcd potentialSum(const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver, const int m1, const int timeIndex); 


	void exponentialTerm(const double time); 
	Eigen::ArrayXXcd matrixExponential(const double time, const int m1);

	template <class T>
	void thetaExponential(const double radius, const double phi, const T & df);
	Eigen::ArrayXXcd thetaExponential(const int m1, const double radius) const;

	Eigen::ArrayXXcd sumOverM1(const int timeEnd, const int timeIndex);

	Eigen::ArrayXXd dfIntegrate(const int timeEnd);
	void updateVectors(const int timeEnd, const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver);

	void outApoPeri(std::ofstream & out) const;
	void outEL(std::ofstream & out) const;
};

template <class T>
std::vector<Eigen::ArrayXXd> DFfunction::dFdJsetUp(const ActionAngleBasisContainer & bf, const T & df)
{
	Eigen::ArrayXXd dFdL{df.dFdJgrid(bf.spacing(), m_omega2).array()}, dFdE{df.dFdEgrid(bf.spacing(), m_omega1).array()};
	std::vector<Eigen::ArrayXXd> dFdJ; 

	for (int m1 = -m_m1Max; m1 <= m_m1Max; ++m1){
		dFdJ.push_back(m1 * dFdE + m_fourierHarmonic * dFdL);
	}
	return dFdJ;
}




template <class T>
Eigen::ArrayXXd DFfunction::theta1(const double radius, const double phi, const T & df) const
{
	Eigen::ArrayXXd theta1 = Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize);
	for (int i = radius2index(radius); i < m_dfSize; ++i){
		for (int j = 1; j < radius2index(radius); ++j){theta1(i,j) = df.theta1(radius, i*m_dfStep, j*m_dfStep, m_omega1(i,j));}
	}
	return theta1;
}

template <class T>
Eigen::ArrayXXd DFfunction::theta2(const double radius, const double phi, const T & df) const
{
	Eigen::ArrayXXd theta2 = Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize);
	for (int i = radius2index(radius); i < m_dfSize; ++i){
		for (int j = 1; j < radius2index(radius); ++j){theta2(i,j) = (df.theta2(radius, i*m_dfStep, j*m_dfStep, m_omega2(i,j))+phi);}
	}
	return theta2;
}

template <class T>
void DFfunction::coordSetUp(const T & df, const double spacing)
{
	for (int i = 0; i < m_dfSize; ++i){
		for (int j = 1; j < i; ++j){
			m_rPeri(i,j) = j*spacing;
			m_rApo(i,j)  = i*spacing;
			m_eCoords(i,j) = df.rad2Energy(m_rApo(i,j), m_rPeri(i,j));
			m_lCoords(i,j) = df.rad2AngMom(m_rApo(i,j), m_rPeri(i,j));
		}
	}
}


void DFfunction::potentialTerm(const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver, const int timeIndex)
{
	for (int m1 = -m_m1Max; m1 <= m_m1Max; ++m1)
	{
		v_potentialTimesdFdJ.push_back(potentialSum(basisFunc, solver, m1, timeIndex));
	}
}

Eigen::ArrayXXcd DFfunction::potentialSum(const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver,  const int m1, const int timeIndex)
{
	Eigen::MatrixXcd sum = Eigen::MatrixXcd::Zero(basisFunc.size(0), basisFunc.size(0));
	for (int n = 0; n <= m_maxRadialIndex; ++n){
		Eigen::VectorXcd coeff{solver.responseCoef(timeIndex) + solver.perturbationCoef(timeIndex)};
		sum += coeff[n] * basisFunc(n, m1); 
	}
	return dFdJGrid(m1) * sum.array();
}


// Same Structure as above

void DFfunction::exponentialTerm(const double time){
	for (int m1 = -m_m1Max; m1 <= m_m1Max; ++m1){ v_exponential.push_back(matrixExponential(time, m1));}
}

Eigen::ArrayXXcd DFfunction::matrixExponential(const double time, const int m1){ // Note that time here really is t - t'
	std::complex<double> unitComplex(0,1);
	Eigen::ArrayXXcd exponetial = Eigen::ArrayXXcd::Zero(m_dfSize, m_dfSize); 

	for (int i = 0 ; i < exponetial.rows(); ++i)
	{
		 for (int j = 1; j < i; ++ j){
		 	exponetial(i,j) = exp(-unitComplex * time * (m1*m_omega1(i,j) + m_fourierHarmonic * m_omega2(i,j)));
		 }
	}
	return exponetial; 
}

template <class T>
void DFfunction::thetaExponential(const double radius, const double phi, const T & df){
	getTheta(radius, phi, df);
	for (int m1 = -m_m1Max; m1 <= m_m1Max; ++m1){
		v_thetaExponential.push_back(thetaExponential(m1, radius));
	}
}

Eigen::ArrayXXcd DFfunction::thetaExponential(const int m1, const double radius) const {
	std::complex<double> unitComplex(0,1);
	Eigen::ArrayXXcd exponetial = Eigen::ArrayXXcd::Zero(m_dfSize, m_dfSize); 

	for (int i = radius2index(radius); i < m_dfSize; ++i){
		for (int j = 1; j < radius2index(radius); ++j){exponetial(i,j) = exp(unitComplex * (m1*m_theta1(i,j) + m_fourierHarmonic*m_theta2(i,j)));}
	}
	return exponetial; 
}


// Integrand and integral terms
Eigen::ArrayXXcd DFfunction::sumOverM1(const int timeEnd, const int timeIndex){
	Eigen::ArrayXXcd integrand = Eigen::ArrayXXcd::Zero(m_dfSize, m_dfSize);

	for (int m1 = -m_m1Max; m1 <= m_m1Max; ++m1){		
		integrand += v_thetaExponential[m1+m_m1Max]*v_exponential[(2*m_m1Max+1)*(timeEnd-timeIndex) + m1+m_m1Max] * 
		v_potentialTimesdFdJ[(2*m_m1Max+1)*timeIndex + m1+m_m1Max]; 	
	}
	return integrand; 
}


// Function the integrates and returns a DF
Eigen::ArrayXXd DFfunction::dfIntegrate(const int timeEnd){
	Eigen::ArrayXXcd integral = Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize);
	for (int timeIndex = 1; timeIndex < timeEnd; ++timeIndex){
		integral += m_timeStep * sumOverM1(timeEnd, timeIndex); 
	}
	return 2*(integral + 0.5 * m_timeStep * sumOverM1(timeEnd, timeEnd)).real();
}


// Function that loops through all the time steps 
void DFfunction::updateVectors(const int timeEnd, const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver)
{
	
	exponentialTerm(timeEnd * m_timeStep);
	potentialTerm(basisFunc, solver, timeEnd); 
	if (timeEnd != 0) {v_DF.push_back(dfIntegrate(timeEnd));}
	else {v_DF.push_back(Eigen::ArrayXXd::Zero(m_dfSize, m_dfSize));}	
}

template <class T>
void DFfunction::dfEvolution(const ActionAngleBasisContainer & basisFunc, const VolterraSolver & solver, const T & df, const double radius, const double phi)
{
	thetaExponential(radius, phi, df);
	updateVectors(0, basisFunc, solver);
	for (int timeEnd = 1; timeEnd < m_numbTimeStep; ++timeEnd){
		if ((timeEnd % 50)  == 0 ){std::cout << "Fraction of DF calculated: " << timeEnd/((double) m_numbTimeStep) << '\n';} 
		updateVectors(timeEnd, basisFunc, solver);
	}
}

void DFfunction::outApoPeri(std::ofstream & out) const{
	
	for (int i =0; i < m_dfSize; ++i){
		for (int j = 1; j < i; ++j)
		{
			if (j == m_dfSize-2) {out << m_rPeri(i, j) << '\n';} 
			else {out << m_rPeri(i, j) << ',';}
		}
	}
	for (int i =0; i < m_dfSize; ++i){
		for (int j = 1; j < i; ++j)
		{
			if (j == m_dfSize-2) {out << m_rApo(i, j) << '\n';} 
			else {out << m_rApo(i, j) << ',';}
		}
	}
}
	
void DFfunction::outEL(std::ofstream & out) const {
	for (int i =0; i < m_dfSize; ++i){
		for (int j = 1; j < i; ++j)
		{
			if (j == m_dfSize-2) {out << m_eCoords(i, j) << '\n';} 
			else {out << m_eCoords(i, j) << ',';}
		}
	}
	for (int i =0; i < m_dfSize; ++i){
		for (int j = 1; j < i; ++j)
		{
			if (j == m_dfSize-2) {out << m_lCoords(i, j) << '\n';} 
			else {out << m_lCoords(i, j) << ',';}
		}
	}
}


void DFfunction::df2fileApoPeri(const std::string & filename, const int skip)
{
	std::ofstream out(filename);
	outApoPeri(out);
	for (int time = skip; time < m_numbTimeStep; time += skip){
		for (int i =0; i < m_dfSize; ++i){
			for (int j = 1; j < i; ++j)
			{
				if (j == m_dfSize-2) {out << v_DF[time](i,j)*m_j2apoPeriJacobian(i,j) << '\n';} // We include the jacobian factor here. 
				else {out << v_DF[time](i,j)*m_j2apoPeriJacobian(i,j) << ',';}
			}
		}
	}
	out.close(); 
}

void DFfunction::df2fileEL(const std::string & filename, const int skip)
{
	std::ofstream out(filename);
	outEL(out);
	for (int time = skip; time < m_numbTimeStep; time += skip){
		for (int i =0; i < m_dfSize; ++i){
			for (int j = 1; j < i; ++j){
				if (j == m_dfSize-2) {out << v_DF[time](i,j)*(1/m_omega1(i,j)) << '\n';} // We include the jacobian factor here. 
				else {out << v_DF[time](i,j)*(1/m_omega1(i,j)) << ',';}
			}
		}
	}
	out.close(); 
}



#endif