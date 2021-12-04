#ifndef CBESOLVER
#define CBESOLVER 

#include <Eigen/Dense>

#include "PerturbedDFContainer.h"

class CBESolver
{
public:
	CBESolver() {} 
	~CBESolver() {}

	// Some function that does the evolution
	// Some function that saves the DF


private: 
	// We need a vector of perturbed DFs, each one representing a different fourier mode

	const int m_maxFourierHarmonic, m_fourierHarmonic; 

	ActionAngleBasisContainer m_potentialCoefficents; 
	ExpansionCoeff m_perturbationCoeff;

	Eigen::ArrayXXcd m_om1Grid, m_om2Grid, m_dfdEGrid, m_dfdJGrid; 
	
	std::vector<Eigen::ArrayXXcd> v_expoOmega, v_nDotDeriv; 
	std::vector<PerturbedDFContainer> v_perturbedDF; 
	const double m_timeStep; 

	std::vector<Eigen::ArrayXXcd> vectorExpoOmega();
	std::vector<Eigen::ArrayXXcd> vectornDotDeriv(); 

};


std::vector<Eigen::ArrayXXcd> CBESolver::vectorExpoOmega() {
	std::vector<Eigen::ArrayXXcd> expoOmega(2*m_maxFourierHarmonic + 1);

	for (int n = -m_maxFourierHarmonic; n <= m m_maxFourierHarmonic; ++n) {expoOmega[m+m_maxFourierHarmonic] = 
		exp(-2 * (unitComplex * m_timeStep) * (n*m_om1Grid + m_fourierHarmonic*m_om2Grid));}
	return expoOmega; 
}

std::vector<Eigen::ArrayXXcd> CBESolver::vectornDotDeriv() {
	std::vector<Eigen::ArrayXXcd> vectornDotDeriv(2*m_maxFourierHarmonic + 1);
	for (int n = -m_maxFourierHarmonic; n <= m m_maxFourierHarmonic; ++n) {vectornDotDeriv[m+m_maxFourierHarmonic] = n*m_dfdEGrid + m_fourierHarmonic * m_om2Grid;}

	return vectornDotDeriv; 
}



void CBESolver::evolution() {
	for (int time = 1; time < (m_numbTimeSteps-1); ++time) {
		for (int n = -m_maxFourierHarmonic; n <= m_maxFourierHarmonic; ++n) {takeEvolutionStep(n, time);}
	}
}

void CBESolver::takeEvolutionStep(const int n, const int time) {
	// Some function that takes an individual step 

}

#endif