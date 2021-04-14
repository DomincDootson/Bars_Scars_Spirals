#ifndef BAR2D
#define BAR2D 
//For sake of argument we assume that it is the positive harmoic

#include <Eigen/Dense>
#include <complex>
#include <iostream>

#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"


class Bar2D
{
public: 
	template <class Tbf>
	Bar2D(const Eigen::VectorXcd expansionCoeff, const Tbf& basisFunctions, const double omega0) :
	m_theta{0}, m_omega{omega0}, m_alpha{0},
	m_expansionCoeff{expansionCoeff},
	m_momentOfInertia{1}, //m_momentOfInertia{momentOfInertia(basisFunctions)},
	m_fourierHarmonic{basisFunctions.fourierHarmonic()},
	m_scriptE{basisFunctions.getScriptE()}
	{v_theta.push_back(m_theta); v_omega.push_back(m_omega); v_alpha.push_back(m_alpha);}

	Bar2D(const Eigen::VectorXcd expansionCoeff, const double omega0) :
	m_theta{0}, m_omega{omega0}, m_alpha{0},
	m_expansionCoeff{expansionCoeff},
	m_momentOfInertia{1},
	m_fourierHarmonic{2},
	m_scriptE{Eigen::MatrixXd::Identity(m_expansionCoeff.size(), m_expansionCoeff.size())}
	{v_theta.push_back(m_theta); v_omega.push_back(m_omega); v_alpha.push_back(m_alpha);}



	~Bar2D() {} 

	Eigen::VectorXcd barCoeff() const {return m_expansionCoeff;}

	double torque(const Eigen::VectorXcd &diskCoeff) const;
	template <class Tbf>
	double momentOfInertia(const Tbf& basisFunctions) const;

	void drift(const double timeStep);
	void kick(const double timeStep, const Eigen::VectorXcd &diskCoeff, const double freelyRotating);

	void saveBarEvolution(const std::string & evolutionFilename, const int skip = 10) const;

	
private: 
	double m_theta, m_omega, m_alpha; // Note that alpha = Torque/MOI.
	Eigen::VectorXcd m_expansionCoeff;
	std::vector<double> v_theta, v_omega, v_alpha; // We use these to keep track of the evolution of the bar 

	const double m_momentOfInertia;  
	const int m_fourierHarmonic;  
	const Eigen::MatrixXd m_scriptE; 

	

	Eigen::ArrayXXd rSquaredArray(const int nGrid, const double rMax) const;
};

void Bar2D::saveBarEvolution(const std::string & evolutionFilename, const int skip) const {
	std::ofstream out(evolutionFilename);
	for (int time = 0; time < v_theta.size(); time += skip){
		out << v_theta[time] << ',' << v_omega[time] << ',' << v_alpha[time] << ',' << m_momentOfInertia * v_alpha[time] << '\n';
	}
	out.close();
}


double Bar2D::torque(const Eigen::VectorXcd &diskCoeff) const // Torque = il * (A.scriptE.conj(B) - b.scriptE.conj(A))
{
	std::complex<double>    firstTerm{ m_expansionCoeff.dot(m_scriptE *  diskCoeff) }; // So the dot product of complex vectors
	std::complex<double> secondTerm{ diskCoeff.dot(m_scriptE *  m_expansionCoeff) }; // has the conjugation built in
	std::complex<double> unitComplex(0,1);
	double l{ (double) m_fourierHarmonic};
	
	return std::real(unitComplex * ((firstTerm - secondTerm) * l)); 

}

Eigen::ArrayXXd Bar2D::rSquaredArray(const int nGrid, const double rMax) const{ 
	Eigen::ArrayXXd rSquaredArray{Eigen::ArrayXXd::Zero(nGrid, nGrid)};
	double x{}, y{}, rSquared{}, spacing{2*rMax/((double) nGrid-1)}, centre{0.5*(nGrid-1)};
	for (int i = 0; i < rSquaredArray.rows(); ++i){
		for (int j = 0; j < rSquaredArray.cols(); ++j){
			x = spacing * (i - centre); y = spacing * (j - centre); rSquared = x*x +y*y;
			rSquaredArray(i,j) = rSquared;
		}
	}
	return rSquaredArray;
}

template <class Tbf>
double Bar2D::momentOfInertia(const Tbf& basisFunctions) const { // Ditto with this, a check wouldn't hurt 
	int nGrid{801}; // Parameters that set the size of the array. 
	double rMax{20}, stepSize{2*rMax/((double) nGrid-1)}, totalMass{.01};

	Eigen::ArrayXXd density{basisFunctions.densityArrayReal(m_expansionCoeff, nGrid, rMax)};
	double normalisation{stepSize*stepSize*density.sum()};
	
	Eigen::ArrayXXd r2Array{rSquaredArray(nGrid, rMax)};
	std::cout << "Moment of Inertia " << totalMass*(stepSize*stepSize*(density*r2Array).sum())/normalisation << '\n';
	return totalMass*(stepSize*stepSize*(density*r2Array).sum())/normalisation;
}

void Bar2D::drift(const double timeStep){
	std::complex<double> unitComplex(0,1);
	double deltaTheta{m_omega*timeStep + 0.5 * m_alpha * timeStep * timeStep}; 
	
	m_expansionCoeff *= exp(-unitComplex * (deltaTheta * m_fourierHarmonic)); 
	m_theta += deltaTheta; 
	v_theta.push_back(m_theta);
}

void Bar2D::kick(const double timeStep, const Eigen::VectorXcd &diskCoeff, const double freelyRotating){
	double oldAlpha = m_alpha;
	m_alpha = torque(diskCoeff)/m_momentOfInertia; v_alpha.push_back(m_alpha); 
	
	m_omega += freelyRotating * 0.5 * (oldAlpha + m_alpha) * timeStep;  
	v_omega.push_back(m_omega);
}


#endif