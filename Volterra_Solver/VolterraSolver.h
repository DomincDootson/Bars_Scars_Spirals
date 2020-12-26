// Vector dep comes along in EvolutionKernels.h
#ifndef VOLTERRASOLVER
#define VOLTERRASOLVER 

#include <Eigen/Dense>

#include "EvolutionKernels.h"
#include "ExpansionCoeff.h"

class VolterraSolver
{
public:
	VolterraSolver(int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, double timeStep) :// We probably need to pass it a DF
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_numbTimeSteps{numbTimeSteps}, m_timeStep{timeStep},
	m_kernels(m_numbTimeSteps),
	m_responseCoef(m_numbTimeSteps), m_perturbationCoef(m_numbTimeSteps) {}

	
	VolterraSolver(const std::string & kernelFilename, int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, double timeStep) :
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_numbTimeSteps{numbTimeSteps}, m_timeStep{timeStep},
	m_kernels(kernelFilename, m_numbTimeSteps),
	m_responseCoef(m_numbTimeSteps), m_perturbationCoef(m_numbTimeSteps) {}
	
	~VolterraSolver() {}

	template <class Tdf>
	void generateKernel(const Tdf & df, const ActionAngleBasisContainer & basisFunc, const Eigen::MatrixXcd & scriptE);// Come up with some standard way of naming kernels
	void volterraSolver(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent = true);

	void evolution();

private:
	const int m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps;
	const double m_timeStep;

	EvolutionKernels m_kernels;
	ExpansionCoeff m_responseCoef, m_perturbationCoef;

	Eigen::VectorXcd timeIntegration(const int timeIndex, const double includeSelfConsistent) const;	
};

template <class Tdf>
void VolterraSolver::generateKernel(const Tdf & df, const ActionAngleBasisContainer & basisFunc, const Eigen::MatrixXcd & scriptE)
{
	m_kernels.getVolterraParams(m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_timeStep);
	m_kernels.kernelCreation(df, basisFunc, scriptE);
}

void VolterraSolver::volterraSolver(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent) 
{
	Eigen::MatrixXcd identity{Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1)};
	double includeSelfConsistent{1};
	if (!isSelfConsistent){
		includeSelfConsistent = 0;
	}

	for (int timeIndex = 1; timeIndex < m_numbTimeSteps; ++timeIndex){
		m_responseCoef(timeIndex) = m_responseCoef(0) + ((identity - includeSelfConsistent*0.5*m_kernels(timeIndex)).inverse()) 
									* timeIntegration(timeIndex, includeSelfConsistent);	
	}
	m_responseCoef.write2File(outFilename);
}

Eigen::VectorXcd VolterraSolver::timeIntegration(const int timeIndex, const double includeSelfConsistent) const // Should we make this a memeber function? 
{
	Eigen::VectorXcd integral = 0.5 * m_timeStep * m_kernels(timeIndex) * (m_perturbationCoef(0) + includeSelfConsistent*m_responseCoef(0));
	
	for (int i = 1; i < (timeIndex); ++i){
		integral += m_timeStep * m_kernels(timeIndex) * (m_perturbationCoef(i) + includeSelfConsistent*m_responseCoef(i)); 
	}
	
	integral += 0.5 * m_timeStep * m_kernels(0) * m_perturbationCoef(timeIndex);	

	return integral; 
}

#endif