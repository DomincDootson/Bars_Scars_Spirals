#ifndef EVOLUTIONKERNELS
#define EVOLUTIONKERNELS 

#include <vector>
#include <complex>

#include <thread>
#include <functional>

#include <Eigen/Dense>

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"

class EvolutionKernels
{
public:
	EvolutionKernels(int numbTimeStep) : m_kernels(numbTimeStep) {}
	EvolutionKernels(std::string kernelFilename, int numbTimeStep) : m_kernels(numbTimeStep) {kernelReadIn(kernelFilename);}


	EvolutionKernels(int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, double timeStep) : m_kernels(numbTimeSteps),
	 m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_numbTimeSteps{numbTimeSteps}, m_timeStep{timeStep} {}

	~EvolutionKernels() {}

	Eigen::MatrixXcd operator()(int timeIndex) const {return m_kernels[timeIndex];} 
	Eigen::MatrixXcd& operator()(int timeIndex) {return m_kernels[timeIndex];} 
	
	template <class Tdf>
	void kernelCreation(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc, const std::complex<double> & omega = 0); 
	void getVolterraParams(const int maxRadialIndex, const int fourierHarmonic, const int numbTimeStep, const double timeStep);

	void kernelWrite2File(const std::string & kernelFilename) const; 
	void includeMassFraction(const double massRatio) {for (auto kernel : m_kernels) {kernel *= massRatio;}}
	void activeFraction(const double xi) {m_activeFraction = xi;for (int time = 0; time < m_numbTimeSteps; ++time) { m_kernels[time] *= m_activeFraction;}}
	void resetActiveFraction() {activeFraction(1/m_activeFraction); m_activeFraction = 1;}

	int numbTimeSteps() const {return m_kernels.size();}
	int numbBF() const {return m_maxRadialIndex;}
	double timeStep() const {return m_timeStep;}

	template <class Tdf> 
	Eigen::MatrixXcd responseMatrix(const std::complex<double> & omega, const Tdf & df, const ActionAngleBasisContainer & basisFunc); 

	Eigen::MatrixXcd responseMatrix(const std::complex<double> & omega) const ; 

private:
	std::vector<Eigen::MatrixXcd> m_kernels;

	// These are all member variables that we will need if generating the kernel. 
	Eigen::MatrixXd m_om1Grid, m_om2Grid, m_dfdEGrid, m_dfdJGrid, m_elJacobian;
	int m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_maxFourierHarmonic;
	double m_timeStep, m_spacing, m_activeFraction;

	template <class Tdf>
	void gridSetUp(const Tdf & df, const ActionAngleBasisContainer & basisFunc);

	void kernelReadIn(const std::string & kernelFilename);
	void kernelSquareReadIn(const std::string & kernelFilename); 
	void kernelDiagReadIn(const std::string & kernelFilename);
	void evolutionParams(int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, int maxFourierHarmonic, double timeStep, double spacing);
	
	
	template <class Tdf>
	void generateKernel(const Tdf & df, const ActionAngleBasisContainer & basisFunc, const std::complex<double> & omega = 0); 

	std::complex<double> kernelElement(const int npRow, const int npCol, const ActionAngleBasisContainer & basisFunc, const double time, const std::complex<double> & omega) const;
	void kernelAtTime(const ActionAngleBasisContainer & basisFunc, const int timeIndex, const std::complex<double> & omega); 
	void multipleKernelAtTime(const ActionAngleBasisContainer & basisFunc, const int timeIndex, const int nCores, const std::complex<double> & omega);

	std::complex<double> integration2d(const Eigen::MatrixXcd & grid2Integrate) const;
};

void EvolutionKernels::getVolterraParams(const int maxRadialIndex, const int fourierHarmonic, const int numbTimeStep, const double timeStep){
	m_maxRadialIndex = maxRadialIndex;
	m_fourierHarmonic = fourierHarmonic;
	m_numbTimeSteps = numbTimeStep;
	m_timeStep = timeStep;
}

template <class Tdf>
void EvolutionKernels::gridSetUp(const Tdf & df, const ActionAngleBasisContainer & basisFunc){
	
	m_spacing = basisFunc.spacing();
	m_om1Grid = df.omega1Grid(basisFunc.size(0), m_spacing); 
	m_om2Grid = df.omega2Grid(basisFunc.size(0), m_spacing);
	
	m_maxFourierHarmonic = basisFunc.maxFourierHarmonic();

	m_dfdEGrid = df.dFdEgrid(basisFunc.spacing(), m_om1Grid); 
	m_dfdJGrid = df.dFdJgrid(basisFunc.spacing(), m_om2Grid); 	

	m_elJacobian = df.energyAngMomJacobain(basisFunc.size(0), basisFunc.spacing());	
}

template <class Tdf>
void EvolutionKernels::generateKernel(const Tdf & df, const ActionAngleBasisContainer & basisFunc,  const std::complex<double> & omega)
{
	Eigen::MatrixXd inverseScriptE{basisFunc.inverseScriptE()};
	gridSetUp(df, basisFunc);
	int nCores{2};

	for (int timeIndex = 0; timeIndex < m_numbTimeSteps; timeIndex += nCores)
	{
		multipleKernelAtTime(basisFunc, timeIndex, nCores, omega);
		if (timeIndex % 10 == 0){
			std::cout << "Fraction of kernels completed: " << round(100*timeIndex/((double) m_numbTimeSteps))<< '%' <<  '\n';
		}
	}	
	//for (int timeIndex = 0; timeIndex < m_numbTimeSteps; ++timeIndex) {m_kernels[timeIndex] = inverseScriptE * m_kernels[timeIndex];} //PUT THIS BACK IN
}

template <class Tdf>
void EvolutionKernels::kernelCreation(const std::string fileName, const Tdf & df, const ActionAngleBasisContainer & basisFunc,  const std::complex<double> & omega)
{
	generateKernel(df, basisFunc, omega); 
	kernelWrite2File(fileName);
}


void EvolutionKernels::multipleKernelAtTime(const ActionAngleBasisContainer & basisFunc, const int timeIndex, const int nCores, const std::complex<double> & omega)
{
	std::vector<std::thread> threads;
	for (int i = 0; i < nCores; ++i){
		threads.push_back(std::thread( [&, this, basisFunc, timeIndex, i, omega] {kernelAtTime(basisFunc, timeIndex +i, omega);} ));
	}
	for (auto &th : threads) {
	th.join();
	}
}

void EvolutionKernels::kernelAtTime(const ActionAngleBasisContainer & basisFunc, const int timeIndex, const std::complex<double> & omega)
{
	m_kernels[timeIndex] = (Eigen::MatrixXcd::Zero(m_maxRadialIndex+1, m_maxRadialIndex+1));
	for (int i = 0; i <= m_maxRadialIndex; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			(m_kernels[timeIndex])(i,j) = kernelElement(i, j, basisFunc, timeIndex * m_timeStep, omega);
			(m_kernels[timeIndex])(j,i) = (m_kernels[timeIndex])(i,j);
		}
	}
}


std::complex<double> EvolutionKernels::kernelElement(const int npRow, const int npCol, const ActionAngleBasisContainer & basisFunc, const double time, const std::complex<double> & omega) const
{
	Eigen::MatrixXcd integrand(basisFunc.size(0), basisFunc.size(1));
	std::complex<double> unitComplex(0,1);
	for (int i = 0; i < integrand.rows(); ++i) 
	{
		for (int j = 1; j <= i; ++j)
		{
			integrand(i,j) = 0;
			for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++m1) 
			{
				integrand(i,j) += basisFunc(npRow, m1, i, j) * basisFunc(npCol, m1, i, j)
				* (m_fourierHarmonic*m_dfdJGrid(i,j) +m1*m_dfdEGrid(i,j)) 
				* exp(-unitComplex*time * ((m1 * m_om1Grid(i,j) + m_fourierHarmonic * m_om2Grid(i,j))))
				* exp(unitComplex * time * omega);	
			}
		}
	}
	return -unitComplex * (2*M_PI)*(2*M_PI) * integration2d(integrand); 
}


std::complex<double> EvolutionKernels::integration2d(const Eigen::MatrixXcd & grid2Integrate) const 
{	
	std::complex<double> integral{0}, rowIntegral{0};
	for (int i = 0; i < grid2Integrate.rows(); ++i)
	{
		for (int j = 1; j < i; ++j)
		{
			rowIntegral +=  m_spacing * grid2Integrate(i,j)*m_elJacobian(i, j)*(1/m_om1Grid(i,j));
		}
		integral +=  m_spacing * rowIntegral;
		rowIntegral = 0;
	}
	return integral;
}

void EvolutionKernels::kernelWrite2File(const std::string & kernelFilename) const
{
	std::ofstream out(kernelFilename);
	out << m_maxRadialIndex << " " << m_fourierHarmonic << " " << m_numbTimeSteps << " " << m_maxFourierHarmonic << " " << m_timeStep << " " << m_spacing << '\n';
	
	for (int time = 0; time < m_kernels.size(); ++ time){
		for (int i = 0; i<=m_maxRadialIndex; ++i){
			for (int j = 0; j <= m_maxRadialIndex; ++j){
				if (j == m_maxRadialIndex && i ==m_maxRadialIndex) {out << real(m_kernels[time](i,j)) << " " << imag(m_kernels[time](i,j)) <<'\n';}
				else {out << real(m_kernels[time](i,j)) << " " << imag(m_kernels[time](i,j)) << " ";}
			}
		}
	}
	out.close();
	std::cout << "Kernel saved to: " << kernelFilename << '\n';
}

void EvolutionKernels::evolutionParams(int maxRadialIndex, int fourierHarmonic, int numbTimeSteps, int maxFourierHarmonic, double timeStep, double spacing){
	m_maxRadialIndex = maxRadialIndex; 
	m_fourierHarmonic = fourierHarmonic;
	m_numbTimeSteps = numbTimeSteps;
	m_maxFourierHarmonic = maxFourierHarmonic;
	m_timeStep = timeStep;
	m_spacing = spacing;
	std::cout << numbTimeSteps <<'\n';
	assert(numbTimeSteps == m_kernels.size() && "The read in kernel is not the same length as the Volterra Solver.");
	assert(numbTimeSteps == m_kernels.size() && "The read in kernel is not the same length as the Volterra Solver.");
}

void EvolutionKernels::kernelReadIn(const std::string & kernelFilename) // Needs checking, but should be okay
{
	std::cout << "Reading in kernel from: " << kernelFilename <<'\n';
	std::ifstream kernelIn(kernelFilename);
	int maxRadialIndex, fourierHarmonic, numbTimeSteps, maxFourierHarmonic;
	double timeStep, spacing; 
	kernelIn >> maxRadialIndex >> fourierHarmonic >> numbTimeSteps >> maxFourierHarmonic >> timeStep >> spacing;
	evolutionParams(maxRadialIndex, fourierHarmonic, numbTimeSteps, maxFourierHarmonic, timeStep, spacing);
	kernelIn.close();

	if (maxFourierHarmonic ==0 && spacing == 0) {kernelDiagReadIn(kernelFilename);} // Diag kernel won't come from AA rep of BF 
	else {kernelSquareReadIn(kernelFilename);}
}


void EvolutionKernels::kernelDiagReadIn(const std::string & kernelFilename) {
	std::cout << "Reading in Diagonal Kernel\n";
	std::ifstream kernelIn(kernelFilename);
	int maxRadialIndex, fourierHarmonic, numbTimeSteps, maxFourierHarmonic;
	double timeStep, spacing, re, im; 
	kernelIn >> maxRadialIndex >> fourierHarmonic >> numbTimeSteps >> maxFourierHarmonic >> timeStep >> spacing;
	std::complex<double> unitComplex(0,1);

	for (int time = 0; time < numbTimeSteps; ++time){
		m_kernels[time].resize(maxRadialIndex+1, maxRadialIndex+1);
		m_kernels[time].setZero(); 
		for (int i = 0; i <= maxRadialIndex; ++ i){
			kernelIn >> re >> im;
			m_kernels[time](i,i) = re + unitComplex*im;
		}
	} 
	kernelIn.close();
}

void EvolutionKernels::kernelSquareReadIn(const std::string & kernelFilename) {
	std::ifstream kernelIn(kernelFilename);
	int maxRadialIndex, fourierHarmonic, numbTimeSteps, maxFourierHarmonic;
	double timeStep, spacing, re, im; 
	kernelIn >> maxRadialIndex >> fourierHarmonic >> numbTimeSteps >> maxFourierHarmonic >> timeStep >> spacing;
	std::complex<double> unitComplex(0,1);

	for (int time = 0; time < numbTimeSteps; ++time){
		m_kernels[time].resize(maxRadialIndex+1, maxRadialIndex+1);
		for (int i = 0; i <= maxRadialIndex; ++ i){
			for (int j = 0; j <= maxRadialIndex; ++j){
				kernelIn >> re >> im;
				m_kernels[time](i,j) = re + unitComplex*im;
			}
		}
	} 
	kernelIn.close();
}


/* Response Matrix */
/* --------------- */

Eigen::MatrixXcd EvolutionKernels::responseMatrix(const std::complex<double> & omega) const {
	Eigen::MatrixXcd rm{0.5 * m_kernels.front()}; bool isStable{omega.imag() < 0};
	
	const std::complex<double> unitComplex(0,1);
	if (isStable)
	{
		const std::complex<double> iOmega{conj(omega) * unitComplex};
		for (int time = 1; time < m_numbTimeSteps-1; ++time) {rm += m_kernels[time].conjugate() * exp(iOmega * (time * m_timeStep));}
		rm += 0.5 * m_kernels.back().conjugate() * exp(iOmega * ((m_numbTimeSteps-1) * m_timeStep));
		return - m_timeStep * rm; 
	}

	else {
		const std::complex<double> iOmega{omega * unitComplex};
		for (int time = 1; time < m_numbTimeSteps-1; ++time) {rm += m_kernels[time] * exp(iOmega * (time * m_timeStep));}
		rm += 0.5 * m_kernels.back() * exp(iOmega * ((m_numbTimeSteps-1) * m_timeStep));
		return m_timeStep*rm; 
	} 
}

template <class Tdf> 
Eigen::MatrixXcd EvolutionKernels::responseMatrix(const std::complex<double> & omega, const Tdf & df, const ActionAngleBasisContainer & basisFunc) {
	bool isStable{omega.imag() > 0}; std::complex<double> omegaConj = std::conj(omega);
	
	if (isStable) {generateKernel(df, basisFunc, omega);}
	else {generateKernel(df, basisFunc, omegaConj);}
	Eigen::MatrixXcd rm{0.5 * m_kernels.front()}; 

	for (int time = 0; time < m_numbTimeSteps-1; ++time) {rm += m_kernels[time];}
	rm += 0.5 * m_kernels.back();

	return m_timeStep * rm.conjugate(); 
} 


#endif




