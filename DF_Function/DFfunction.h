#ifndef DFFUNCTION
#define DFFUNCTION

#include <vector>
#include <complex>

#include <Eigen/Dense>

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/DFClass.h"
#include "../Volterra_Solver/VolterraSolver.h"

class DFfunction
{
public:
	DFfunction() {}
	~DFfunction() {}

private: 
	const int m_m1Max, m_fourierHarmonic, m_maxRadialIndex;
	const Eigen::ArrayXXcd m_dFdJ, m_dFdE, m_Omega1, m_Omega2; 
	std::vector<Eigen::ArrayXXcd> v_exponential, v_potentialTimesdFdJ; 


	std::vector<Eigen::ArrayXXcd> dFdJGrid(const int m1) {return m_fourierHarmonic * m_dFdJ + m1 * m_dFdE;}

	void potentialTerm(const VolterraSolver & solver, const ActionAngleBasisContainer & basisFunc); 
	Eigen::ArrayXXcd potentialSum(const int m1, const int timeIndex); 


	void exponentialTerm(const double time); 
	Eigen::ArrayXXcd matrixExponential(const double time, const int m1);

	Eigen::ArrayXXcd sumOverM1(const int timeIndex);
};


void DFfunction::potentialTerm(const VolterraSolver & solver, const ActionAngleBasisContainer & basisFunc)
{
	for (int m1 = -m_m1Max; m1 <= m_m1Max; ++m1){v_potentialTimesdFdJ.push_back(potentialSum(solver, basisFunc, m1, timeIndex));}
}

Eigen::ArrayXXcd DFfunction::potentialSum(const VolterraSolver & solver, const ActionAngleBasisContainer & basisFunc,  const int m1, const int timeIndex)
{
	Eigen::MatrixXcd sum = Eigen::MatrixXcd::Zero(basisFunc.size(0));
	for (int n = 0; n <= m_maxRadialIndex; ++n){
		sum += (solver.responseCoef(timeIndex) + solver.perturbationCoef(timeIndex)) * basisFunc(n, m1)
	}
	return dFdJGrid(m1) * sum.array();
}


// Same Structure as above

void DFfunction::exponentialTerm(const double time){
	for (int m1 = -m_m1Max; m1 < = m_m1Max; ++m1){ v_exponential.push_back(matrixExponential(time, m1));}
}

Eigen::ArrayXXcd DFfunction::matrixExponential(const double time, const int m1){
	std::complex<double> unitComplex(0,1);
	Eigen::ArrayXXcd exponetial = Eigen::ArrayXXcd::Zero(size, size); //FIGURE OUT THE SIZE

	for (int i = 0 ; i < exponetial.rows(); ++i)
	{
		 for (int j = 1; j < i; ++ j){
		 	exponetial(i,j) = exp(-unitComplex * time * (m1*m_Omega1(i,j) + m_fourierHarmonic * m_Omega2(i,j)));
		 }
	}
	return exponetial; 
}

// Integrand and integral terms
Eigen::ArrayXXcd DFfunction::sumOverM1(const int timeIndex){
	Eigen::ArrayXXcd integrand = Eigen::ArrayXXcd::Zero(v_potentialTimesdFdJ[0].cols(), v_potentialTimesdFdJ[0].rows())

	for (int m1 = -m_m1Max; m1 <= m_m1Max; ++m1){
		integrand += v_exponential[(2*m_m1Max+1)*timeIndex + m1] * v_potentialTimesdFdJ[(2*m_m1Max+1)*timeIndex + m1];
	}
	return integrand; 
}

#endif