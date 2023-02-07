#ifndef VOLTERRASOLVER
#define VOLTERRASOLVER 

#include <Eigen/Dense>

#include "EvolutionKernels.h"
#include "ExpansionCoeff.h"

#include "../Bar2D/Bar2D.h"
#include "../Spiral2D/Spiral2D.h"
#include "../Wave/Wave.h"

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





	/* General Functions */
	/* ----------------- */ 

	int maxRadialIndex() const {return m_maxRadialIndex+1;} // Note this returns due to our silly convention of using <=..., better decisions could have been made in the past
	int numbTimeSteps() const {return m_numbTimeSteps;}
	double timeStep() const {return m_timeStep;}
	void activeFraction(double xi);
	void resetActiveFraction() {activeFraction(1/m_xi);}
	Eigen::MatrixXcd operator()(int timeIndex) const {return m_kernels(timeIndex);} 


	template <class Tdf>
	void generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc);

	template <class Tdf>
	void generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc, const double massRatio);


	void solveVolterraEquation(const bool isSelfConsistent, const int integrateFromIndex = 1); 
	void setInitialCondition(const Eigen::VectorXcd & ic) {m_responseCoef(0) = ic;}
	void setInitialCondition(const std::vector<Eigen::VectorXcd> & ic);

	/* Perturbation Evolution */ 
	/* ---------------------- */ 

	void solveVolterraEquationPerturbation(const std::string & perturbationFilename, const bool isSelfConsistent);
	void coefficentEvolution(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent = true);
	void setInitalPerturbation(const Eigen::VectorXcd & coeff, int time = 0) {m_perturbationCoef(time) = coeff;}


	template <class Tbf>
	void densityEvolution(const Tbf & bf, const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent = true);
	template <class Tbf>
	std::vector<double> energyEvolution(const Tbf & bf, const std::string & perturbationFilename, const bool isSelfConsistent = true);
	
	void kernelTesting(const int nudgeCoef, const bool isSelfConsistent) {m_perturbationCoef(0)[nudgeCoef] = 1; solveVolterraEquation(isSelfConsistent);}
	void kernelTesting(const std::string & filename, const int nudgeCoef, const bool isSelfConsistent) {kernelTesting(nudgeCoef,isSelfConsistent); m_responseCoef.write2File(filename);}
	
	void deltaPerturbationTest();
	void deltaPerturbationConsistent();

	/* Saving to File */
	/* -------------- */ 

	/*template <class Tbf> 
	void writeDensity2File(const std::string & outFilename, Tbf & bf) const {m_responseCoef.write2dDensity2File(outFilename, bf, m_skip); std::cout << "Density evolution written to: " << outFilename <<'\n';}*/ 
	


	template <class Tbf> 
	void density1dEvolution(const std::string & outFilename, const Tbf & bf, const int skip = 1) const {m_responseCoef.writeDensity2File(outFilename, bf, skip);}
	template <class Tbf>
	void density1dCorotating(const std::string & outFilename, const Tbf & bf, const double omegaP, const double innerR = 0, const double outerR = 15, const int skip = 1) const; 
	template <class Tbf>
	void density2dEvolution(const std::string & outFilename, const Tbf & bf, const int skip = 1, const double rMax = 10) const {m_responseCoef.write2dDensity2File(outFilename, bf, skip, rMax);}



	template <class Tbf>
	void density2dEvolution(const int timeIndex, const std::string & outFilename, const Tbf & bf, const double rMax, const int nStep) const {m_responseCoef.write2dDensity2File(timeIndex,outFilename, bf, rMax,nStep);}
	template <class Tbf>
	void potential2dEvolution(const std::string & outFilename, const Tbf & bf, const int skip = 1, const double rMax = 10) const {m_responseCoef.write2dPotential2File(outFilename, bf, skip, rMax);}

	template<class Tbf> 
	double maxDensity(const Tbf & bf, const double rMax =15, const int nGrid=201) const  {return m_responseCoef.maxDensity(bf, rMax, nGrid);}

	
	void saveResponseCoeff(const std::string & filename) const {m_responseCoef.write2File(filename, 1);}
	void saveFinalResponseCoeff(std::ofstream & out) const {m_responseCoef.write2File(m_numbTimeSteps-1, out);}



	/* Bar Evolution */
	/* ------------- */
	void barRotation(Bar2D & bar, const std::string & outFilename, const std::string & evolutionFilename, const bool isSelfConsistent = true, const bool isFreelyRotating = true, const bool isEvolving = false); 
	void barRotationUnsaving(Bar2D & bar, const bool isSelfConsistent = true, const bool isFreelyRotating = false, const bool isEvolving = true); 
	void barRotation(Bar2D & bar, const bool isSelfConsistent);

	/* Spiral */ 
	/* ------ */ 

	template <class T> 
	void spiralEvolution(T & spiral);  // This can also be used for Wave class
	

private:
	const int m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_skip;
	const double m_timeStep;
	double m_xi; 

	EvolutionKernels m_kernels;
	ExpansionCoeff m_responseCoef, m_perturbationCoef;

	Eigen::VectorXcd timeIntegration(const int timeIndex, const double includeSelfConsistent) const;	
	
	void printTimeIndex(const int timeIndex) {if (timeIndex % (m_skip*10) == 0) {std::cout << "Time step: " << timeIndex << '\n';}}
	double selfConsistentDouble(bool isSelfConsistent) {if (isSelfConsistent) {return 1;} else {return 0;}}
	double freelyRotatingDouble(bool isFreelyRotating) {if (isFreelyRotating) {return 1;} else {return 0;}}

	template <class T>
	void transferCoeff(T & spiral) {for (int time = 1; time < m_numbTimeSteps; ++time) {spiral(time) = m_responseCoef(time);}}
};


/* General Functions */
/* ----------------- */ 

void VolterraSolver::activeFraction(const double xi)
{
	std::cout << "Active fraction: " << xi*m_xi<< "\n";
	m_xi = xi; 
	for (int time = 0; time < m_numbTimeSteps; ++time){
		m_kernels(time) *= m_xi; 
	}
}

template <class Tdf>
void VolterraSolver::generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc){
	m_kernels.getVolterraParams(m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_timeStep);
	m_kernels.kernelCreation(fileName, df, basisFunc);
}

template <class Tdf>
void VolterraSolver::generateKernel(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc, const double massRatio) {
	generateKernel(fileName, df, basisFunc);
	m_kernels.includeMassFraction(massRatio);
	m_kernels.kernelWrite2File(fileName); 
}


Eigen::VectorXcd VolterraSolver::timeIntegration(const int timeIndex, const double includeSelfConsistent) const 
{
	Eigen::VectorXcd integral = 0.5 * m_timeStep * m_kernels(timeIndex) * (m_perturbationCoef(0) + includeSelfConsistent*m_responseCoef(0));
	for (int i = 1; i < (timeIndex); ++i){
		integral += m_timeStep * m_kernels(timeIndex - i) * (m_perturbationCoef(i) + includeSelfConsistent*m_responseCoef(i)); 
	}
	integral += 0.5 * m_timeStep * m_kernels(0) * m_perturbationCoef(timeIndex);	
	return integral; 
}



void VolterraSolver::solveVolterraEquation(const bool isSelfConsistent, const int integrateFromIndex) {
	Eigen::MatrixXcd identity{Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1)};
	double includeSelfConsistent{selfConsistentDouble(isSelfConsistent)};
	for (int timeIndex = integrateFromIndex; timeIndex < m_numbTimeSteps; ++timeIndex){
		printTimeIndex(timeIndex);
		m_responseCoef(timeIndex) = m_responseCoef(integrateFromIndex-1) + ((identity - includeSelfConsistent*0.5*m_kernels(timeIndex)).inverse()) 
									* timeIntegration(timeIndex, includeSelfConsistent);
	}
} 

void VolterraSolver::setInitialCondition(const std::vector<Eigen::VectorXcd> & ic) {
	for (int time = 0; time < ic.size(); ++time) {m_responseCoef(time) = ic[time];}
}

/* Perturbation Evolution */ 
/* ---------------------- */ 

void VolterraSolver::solveVolterraEquationPerturbation(const std::string & perturbationFilename, const bool isSelfConsistent)
{
	m_perturbationCoef.coefficentReadIn(perturbationFilename);
	solveVolterraEquation(isSelfConsistent); 
}



void VolterraSolver::coefficentEvolution(const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent) {
	solveVolterraEquationPerturbation(perturbationFilename, isSelfConsistent);
	m_responseCoef.write2File(outFilename, m_skip);
}

template <class Tbf>
void VolterraSolver::densityEvolution(const Tbf & bf, const std::string & outFilename, const std::string & perturbationFilename, const bool isSelfConsistent) {
	solveVolterraEquationPerturbation(perturbationFilename, isSelfConsistent);
	m_responseCoef.writeDensity2File(outFilename, bf, m_skip);
} 


template <class Tbf>
std::vector<double> VolterraSolver::energyEvolution(const Tbf & bf, const std::string & perturbationFilename, const bool isSelfConsistent) {
	solveVolterraEquationPerturbation(perturbationFilename, isSelfConsistent);
	return m_responseCoef.energyEvolution(bf.scriptE());
}

void VolterraSolver::deltaPerturbationTest() {
	for (int t = 0; t < m_numbTimeSteps; ++t) {m_responseCoef(t) = m_kernels(t) * m_perturbationCoef(0);}
}

void VolterraSolver::deltaPerturbationConsistent() {
	deltaPerturbationTest(); 
	Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(m_maxRadialIndex+1, m_maxRadialIndex+1);
	for (int timeIndex = 1; timeIndex < m_numbTimeSteps; ++timeIndex) {
		printTimeIndex(timeIndex); 
		Eigen::VectorXcd integral = 0.5 * m_timeStep * m_kernels(timeIndex) *  m_responseCoef(0);
		for (int i = 1; i < (timeIndex); ++i){ integral += m_timeStep * m_kernels(timeIndex - i) * m_responseCoef(i); }
		integral += m_responseCoef(timeIndex);	
		m_responseCoef(timeIndex) =  ((identity - 0.5*m_kernels(timeIndex)).inverse()) * integral;  
		std::cout << m_responseCoef(timeIndex).dot(m_responseCoef(timeIndex)) <<'\n';
	}
}


/* Bar Evolution */
/* ------------- */

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

void VolterraSolver::barRotationUnsaving(Bar2D & bar, const bool isSelfConsistent, const bool isFreelyRotating, const bool isEvolving)
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
}

void VolterraSolver::barRotation(Bar2D & bar, const bool isSelfConsistent)
{
	for (int timeIndex = 0; timeIndex < m_numbTimeSteps; ++timeIndex) {
		m_perturbationCoef(timeIndex) = bar.barCoeff(-1); 
		bar.drift(m_timeStep, 0); 
	}
	solveVolterraEquation(isSelfConsistent);
}


/* Spiral Function */
/* --------------- */ 

template <class T>
void VolterraSolver::spiralEvolution(T & spiral) {
	spiral(0) = m_xi * spiral(0);
	m_responseCoef(0) = spiral(0); 
	solveVolterraEquation(true);
	
	spiral.resizeVector(m_numbTimeSteps);
	transferCoeff(spiral);
}


/* Saving Functions */ 
/* ---------------- */ 

void saveVector(std::ofstream & out, const std::vector<double> & vec) {
	for (auto i = vec.begin(); i != vec.end()-1; ++i){out << *i<< ','; }
	out << vec.back() << '\n';
}

template <class Tbf>
void VolterraSolver::density1dCorotating(const std::string & outFilename, const Tbf & bf, const double omegaP, const double innerR, const double outerR, const int skip) const 
{
	std::ofstream out(outFilename);
	double spacing{(outerR-innerR)*0.005}; 
	std::vector<double> radii; 
	
	for (double r = innerR; r <= outerR; r += spacing) {radii.emplace_back(r); }
	saveVector(out, radii);
	std::complex<double> unitComplex(0,1);

	for (int t=0; t<m_numbTimeSteps; t += skip) 
	{ 
		std::vector<double> den; 
		den.reserve(radii.size());
		for (auto r : radii) {
			double sum{0};
			for (int n = 0; n<= bf.maxRadialIndex(); ++n) {sum += 2 * bf.density(r, n)
				* (exp(unitComplex*(m_fourierHarmonic * omegaP * m_timeStep*t))* m_responseCoef(t)(n)).real();}
			den.emplace_back(sum); 
		}
		saveVector(out, den);
	}
	out.close(); 

} 



#endif



