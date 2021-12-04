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
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_numbTimeSteps{numbTimeSteps}, m_skip{1}, m_timeStep{timeStep}, m_xi{1}, 
	m_kernels(m_numbTimeSteps),
	m_responseCoef(m_numbTimeSteps, m_maxRadialIndex), m_perturbationCoef(m_numbTimeSteps, m_maxRadialIndex) {}

	
	VolterraSolver(const std::string & kernelFilename, int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, double timeStep) :
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_numbTimeSteps{numbTimeSteps}, m_skip{1}, m_timeStep{timeStep}, m_xi{1},
	m_kernels(kernelFilename, m_numbTimeSteps),
	m_responseCoef(m_numbTimeSteps, m_maxRadialIndex), m_perturbationCoef(m_numbTimeSteps, m_maxRadialIndex) {} 
	 
	~VolterraSolver() {}

	Eigen::VectorXcd     responseCoef(const int timeIndex) const {return     m_responseCoef(timeIndex);}
	Eigen::VectorXcd perturbationCoef(const int timeIndex) const {return m_perturbationCoef(timeIndex);}

	int numbTimeSteps() const {return m_numbTimeSteps;}
	double timeStep() const {return m_timeStep;}

	template <class Tdf>
	void generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc);// Come up with some standard way of naming kernels
	
	void solveVolterraEquation(const std::string & perturbationFilename, const bool isSelfConsistent);
	void solveVolterraEquationDecoupled(const std::string & perturbationFilename, const bool isSelfConsistent, const int index);
	void coefficentEvolution(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent = true);
	
	template <class Tbf>
	void densityEvolution(const Tbf & bf, const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent = true);

	template <class Tbf>
	std::vector<double> energyEvolution(const Tbf & bf, const std::string & perturbationFilename, const bool isSelfConsistent = true);

	template <class Tbf>
	std::vector<double> energyEvolutionDecoupled(const Tbf & bf, const std::string & perturbationFilename, const int index, const bool isSelfConsistent = true);

	void barRotation(Bar2D & bar, const std::string & outFilename, const std::string & evolutionFilename, const bool isSelfConsistent = true, const bool isFreelyRotating = true, const bool isEvolving = false);
	template <class Tbf> 
	void writeDensity2File(const std::string & outFilename, Tbf & bf) const {m_responseCoef.write2dDensity2File(outFilename, bf, m_skip); std::cout << "Density evolution written to: " << outFilename <<'\n';}

	void activeFraction(double xi);
	void resetActiveFraction() {activeFraction(1/m_xi);}

	template <class Tbf>
	void density2dEvolution(const std::string & outFilename, const Tbf & bf, const int skip = 1) const {m_responseCoef.write2dDensity2File(outFilename, bf, skip);}

	template <class Tbf>
	void potential2dEvolution(const std::string & outFilename, const Tbf & bf, const int skip = 1) const {m_responseCoef.write2dPotential2File(outFilename, bf, skip);}

	Eigen::VectorXcd longTimeCoeff(const Eigen::VectorXcd & perturbationCoeff, bool isSelfConsistent = true) const;

	void saveResponseCoeff(const std::string & filename) const {m_responseCoef.write2File(filename, 1);}
	void saveKernelForPython(const std::string & filename) const {m_kernels.saveForPython(filename);} 

	void kernelDecouple(const int index, const bool isZero, const bool isDiag) {m_kernels.kernelDecouple(index, isZero, isDiag);}

private:
	const int m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_skip;
	const double m_timeStep;
	double m_xi; 

	EvolutionKernels m_kernels;
	ExpansionCoeff m_responseCoef, m_perturbationCoef;

	Eigen::VectorXcd timeIntegration(const int timeIndex, const double includeSelfConsistent) const;	
	Eigen::VectorXcd timeIntegrationDecoupled(const int timeIndex, const double includeSelfConsistent,const int index) const;
	
	void printTimeIndex(const int timeIndex) {if (timeIndex % (m_skip*10) == 0) {std::cout << "Time step: " << timeIndex << '\n';}}
	double selfConsistentDouble(bool isSelfConsistent) {if (isSelfConsistent) {return 1;} else {return 0;}}
	double freelyRotatingDouble(bool isFreelyRotating) {if (isFreelyRotating) {return 1;} else {return 0;}}
};

template <class Tdf>
void VolterraSolver::generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc)
{
	m_kernels.getVolterraParams(m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_timeStep);
	m_kernels.kernelCreation(fileName, df, basisFunc);
}

void VolterraSolver::solveVolterraEquation(const std::string & perturbationFilename, const bool isSelfConsistent)
{
	m_perturbationCoef.coefficentReadIn(perturbationFilename);
	Eigen::MatrixXcd identity{Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1)};
	double includeSelfConsistent{selfConsistentDouble(isSelfConsistent)};
	for (int timeIndex = 1; timeIndex < m_numbTimeSteps; ++timeIndex){
		printTimeIndex(timeIndex);
		m_responseCoef(timeIndex) = m_responseCoef(0) + ((identity - includeSelfConsistent*0.5*m_kernels(timeIndex)).inverse()) 
									* timeIntegration(timeIndex, includeSelfConsistent);
	}
}

void VolterraSolver::solveVolterraEquationDecoupled(const std::string & perturbationFilename, const bool isSelfConsistent, const int index)
{
	std::cout << "We are doing this\n";
	m_perturbationCoef.coefficentReadIn(perturbationFilename);
	Eigen::MatrixXcd identity{Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1)};
	double includeSelfConsistent{selfConsistentDouble(isSelfConsistent)};
	for (int timeIndex = 1; timeIndex < m_numbTimeSteps; ++timeIndex){
		printTimeIndex(timeIndex);
		m_responseCoef(timeIndex) = m_responseCoef(0) + ((identity - includeSelfConsistent*0.5*m_kernels(timeIndex)).inverse()) 
									* timeIntegrationDecoupled(timeIndex, includeSelfConsistent, index);
	}
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

Eigen::VectorXcd VolterraSolver::timeIntegrationDecoupled(const int timeIndex, const double includeSelfConsistent,const int index) const // Should we make this a memeber function? 
{
	Eigen::VectorXcd integral = 0.5 * m_timeStep * (m_kernels(timeIndex) * m_perturbationCoef(0) + includeSelfConsistent*m_kernels.kernelDecouple(timeIndex,index) *m_responseCoef(0));
	for (int i = 1; i < (timeIndex); ++i){
		integral += m_timeStep * (m_kernels(timeIndex - i) * m_perturbationCoef(i) + includeSelfConsistent*m_kernels.kernelDecouple(timeIndex - i,index)*m_responseCoef(i)); 
	}
	integral += 0.5 * m_timeStep * m_kernels(0) * m_perturbationCoef(timeIndex);	
	return integral; 
}


void VolterraSolver::coefficentEvolution(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent) 
{
	solveVolterraEquation(perturbationFilename, isSelfConsistent);
	m_responseCoef.write2File(outFilename, m_skip);
}

template <class Tbf>
void VolterraSolver::densityEvolution(const Tbf & bf, const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent)
{
	solveVolterraEquation(perturbationFilename, isSelfConsistent);
	m_responseCoef.writeDensity2File(outFilename, bf, m_skip);
} 


template <class Tbf>
std::vector<double> VolterraSolver::energyEvolution(const Tbf & bf, const std::string & perturbationFilename, const bool isSelfConsistent)
{
	solveVolterraEquation(perturbationFilename, isSelfConsistent);
	return m_responseCoef.energyEvolution(bf.scriptE());
}

template <class Tbf>
std::vector<double> VolterraSolver::energyEvolutionDecoupled(const Tbf & bf, const std::string & perturbationFilename, const int index, const bool isSelfConsistent) {
	solveVolterraEquationDecoupled(perturbationFilename, isSelfConsistent, index);
	return m_responseCoef.energyEvolution(bf.scriptE());		
}

//This function could definitely be tydied up 
void VolterraSolver::barRotation(Bar2D & bar, const std::string & outFilename, const std::string & evolutionFilename, const bool isSelfConsistent, const bool isFreelyRotating, const bool isEvolving)
{
	Eigen::MatrixXcd identity{Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1)};
	double includeSelfConsistent{selfConsistentDouble(isSelfConsistent)}; double freelyRotating{freelyRotatingDouble(isFreelyRotating)};
	if (isEvolving) {m_perturbationCoef(0) = bar.barCoeff(0);}
	else {m_perturbationCoef(0) = bar.barCoeff();}
	for (int timeIndex = 1; timeIndex < m_numbTimeSteps; ++timeIndex){
		printTimeIndex(timeIndex);
		bar.drift(m_timeStep, freelyRotating);

		if (isEvolving) {m_perturbationCoef(timeIndex) = bar.barCoeff(timeIndex * m_timeStep);}
		else {m_perturbationCoef(timeIndex) = bar.barCoeff();}
		m_responseCoef(timeIndex) = m_responseCoef(0) + ((identity - includeSelfConsistent*0.5*m_kernels(timeIndex)).inverse()) 
									* timeIntegration(timeIndex, includeSelfConsistent);

		bar.kick(m_timeStep, m_responseCoef(timeIndex), freelyRotating, m_timeStep * timeIndex);

	}
	m_responseCoef.write2File(outFilename);
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

Eigen::VectorXcd VolterraSolver::longTimeCoeff(const Eigen::VectorXcd & perturbationCoeff, bool isSelfConsistent) const {
	Eigen::MatrixXd id = Eigen::MatrixXd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1);
	if (isSelfConsistent){ return 2*M_PI*(id - 2*M_PI*m_kernels(0)).inverse() * (m_kernels(0) * perturbationCoeff);}
	return 2*M_PI * (m_kernels(0) * perturbationCoeff);
}

#endif