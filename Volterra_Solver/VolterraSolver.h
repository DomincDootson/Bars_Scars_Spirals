// Vector dep comes along in EvolutionKernels.h
#ifndef VOLTERRASOLVER
#define VOLTERRASOLVER 

#include <Eigen/Dense>

#include "EvolutionKernels.h"
#include "ExpansionCoeff.h"
#include "../Bar2D/Bar2D.h"

class VolterraSolver
{
public:
	VolterraSolver(int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, double timeStep) :// We probably need to pass it a DF
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_numbTimeSteps{numbTimeSteps}, m_skip{10}, m_timeStep{timeStep}, m_xi{1}, 
	m_kernels(m_numbTimeSteps),
	m_responseCoef(m_numbTimeSteps, m_maxRadialIndex), m_perturbationCoef(m_numbTimeSteps, m_maxRadialIndex) {}

	
	VolterraSolver(const std::string & kernelFilename, int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, double timeStep) :
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_numbTimeSteps{numbTimeSteps}, m_skip{10}, m_timeStep{timeStep}, m_xi{1},
	m_kernels(kernelFilename, m_numbTimeSteps),
	m_responseCoef(m_numbTimeSteps, m_maxRadialIndex), m_perturbationCoef(m_numbTimeSteps, m_maxRadialIndex) {} 
	
	~VolterraSolver() {}

	template <class Tdf>
	void generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc);// Come up with some standard way of naming kernels
	void volterraSolver(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent = true);
	
	void barRotation(Bar2D & bar, const std::string & outFilename, const std::string & evolutionFilename, const bool isSelfConsistent = true);


	void activeFraction(double xi);
	void resetActiveFraction() {activeFraction(1/m_xi);}

private:
	const int m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_skip;
	const double m_timeStep;
	double m_xi; 

	EvolutionKernels m_kernels;
	ExpansionCoeff m_responseCoef, m_perturbationCoef;

	Eigen::VectorXcd timeIntegration(const int timeIndex, const double includeSelfConsistent) const;	
	
	void printTimeIndex(const int timeIndex) {if (timeIndex % (m_skip*10) ==0) {std::cout << "Time step: " << timeIndex << '\n';}}
	double selfConsistentDouble(bool isSelfConsistent) {if (isSelfConsistent) {return 1;} else {return 0;}}
};

template <class Tdf>
void VolterraSolver::generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc)
{
	m_kernels.getVolterraParams(m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_timeStep);
	m_kernels.kernelCreation(fileName, df, basisFunc);
}

void VolterraSolver::volterraSolver(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent) 
{
	m_perturbationCoef.coefficentReadIn(perturbationFilename);
	Eigen::MatrixXcd identity{Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1)};

	double includeSelfConsistent{selfConsistentDouble(isSelfConsistent)};

	for (int timeIndex = 1; timeIndex < m_numbTimeSteps; ++timeIndex){
		printTimeIndex(timeIndex);
		m_responseCoef(timeIndex) = m_responseCoef(0) + ((identity - includeSelfConsistent*0.5*m_kernels(timeIndex)).inverse()) 
									* timeIntegration(timeIndex, includeSelfConsistent);
	}
	m_responseCoef.write2File(outFilename, m_skip);
}

Eigen::VectorXcd VolterraSolver::timeIntegration(const int timeIndex, const double includeSelfConsistent) const // Should we make this a memeber function? 
{
	Eigen::VectorXcd integral = 0.5 * m_timeStep * m_kernels(timeIndex) * (m_perturbationCoef(0) + includeSelfConsistent*m_responseCoef(0));
	for (int i = 1; i < (timeIndex); ++i){
		integral += m_timeStep * m_kernels(timeIndex - i) * (m_perturbationCoef(i) + includeSelfConsistent*m_responseCoef(i)); 
	}
	integral += 0.5 * m_timeStep * m_kernels(0) * m_perturbationCoef(timeIndex);	
	return integral; 
}


void VolterraSolver::barRotation(Bar2D & bar, const std::string & outFilename, const std::string & evolutionFilename, const bool isSelfConsistent)
{
	Eigen::MatrixXcd identity{Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1)};
	double includeSelfConsistent{selfConsistentDouble(isSelfConsistent)};
	
	for (int timeIndex = 1; timeIndex < m_numbTimeSteps; ++timeIndex){
		printTimeIndex(timeIndex);

		bar.drift(m_timeStep);
		m_perturbationCoef(timeIndex) = bar.barCoeff(); // Could overload this function?? 

		m_responseCoef(timeIndex) = m_responseCoef(0) + ((identity - includeSelfConsistent*0.5*m_kernels(timeIndex)).inverse()) 
									* timeIntegration(timeIndex, includeSelfConsistent);

		bar.kick(m_timeStep, m_responseCoef(timeIndex));
	}
	m_perturbationCoef.write2File(outFilename);
	bar.saveBarEvolution(evolutionFilename);
}


void VolterraSolver::activeFraction(const double xi)
{
	std::cout << "Active fraction: " << xi*m_xi<< "\n";
	m_xi = xi; 
	for (int time = 0; time < m_numbTimeSteps; ++time){
		m_kernels(time) *= m_xi; 
	}
}


#endif