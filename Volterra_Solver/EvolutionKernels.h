#ifndef EVOLUTIONKERNELS
#define EVOLUTIONKERNELS 

#include <vector>
#include <complex>

#include <Eigen/Dense>

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "VolterraSolver.h"

class EvolutionKernels
{
public:
	EvolutionKernels(int numbTimeStep) : m_kernels(numbTimeStep) {}
	EvolutionKernels(std::string kernelFilename, int numbTimeStep) : m_kernels(numbTimeStep) {kernelReadIn(kernelFilename);}

	~EvolutionKernels() {}

	Eigen::MatrixXcd operator()(int timeIndex) const {return m_kernels[timeIndex];} 
	
	template <class Tdf>
	void kernelCreation(const Tdf & df, const ActionAngleBasisContainer & basisFunc, const Eigen::MatrixXcd & scriptE); // SHould we pass it a script E here? and also a reference to a DF. --> Put sctio
	void getVolterraParams(const int maxRadialIndex, const int fourierHarmonic, const int numbTimeStep, const double timeStep);

private:
	std::vector<Eigen::MatrixXcd> m_kernels;

	// These are all member variables that we will need if generating the kernel. 
	Eigen::MatrixXd m_om1Grid, m_om2Grid, m_dfdEGrid, m_dfdJGrid, m_elJacobian;
	int m_maxRadialIndex, m_fourierHarmonic, m_numbTimeSteps, m_maxFourierHarmonic;
	double m_timeStep, m_spacing;

	template <class Tdf>
	void gridSetUp(const Tdf & df, const ActionAngleBasisContainer & basisFunc);

	
	void kernelReadIn(const std::string & kernelFilename);
	void kernelWrite2File() const; //const std::string & kernelFilename) WE NEED TO ADD IN FILENAMES
	
	std::complex<double> kernelElement(const int npRow, const int npCol, const ActionAngleBasisContainer & basisFunc, const double time) const;
	void kernelAtTime(const ActionAngleBasisContainer & basisFunc, const int timeIndex); 

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
	m_om1Grid = basisFunc.omega1Grid(df); 
	m_om2Grid = basisFunc.omega2Grid(df); 
	m_spacing = basisFunc.spacing();
	m_maxFourierHarmonic = basisFunc.maxFourierHarmonic();

	m_dfdEGrid = df.dfdEGrid(basisFunc.spacing(), m_om1Grid); 
	m_dfdJGrid = df.dfdJGrid(basisFunc.spacing(), m_om2Grid); 
	m_elJacobian = df.energyAngMomJacobain(basisFunc.size(0), basisFunc.spacing());	
}

template <class Tdf>
void EvolutionKernels::kernelCreation(const Tdf & df, const ActionAngleBasisContainer & basisFunc, const Eigen::MatrixXcd & scriptE)
{
	// Might it be worth making sure that we've read in the construction parameters? 
	gridSetUp(df, basisFunc);
	for (int timeIndex = 0; timeIndex < m_numbTimeSteps; ++timeIndex)
	{
		kernelAtTime(basisFunc, timeIndex);
	}

	kernelWrite2File();
}

void EvolutionKernels::kernelAtTime(const ActionAngleBasisContainer & basisFunc, const int timeIndex)
{
	m_kernels[timeIndex] = Eigen::MatrixXcd::Zero(m_maxRadialIndex, m_maxRadialIndex);
	for (int i = 0; i < m_maxRadialIndex; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			m_kernels[timeIndex](i,j) = kernelElement(i, j, basisFunc, timeIndex * m_timeStep);
			m_kernels[timeIndex](j,i) = m_kernels[timeIndex](i,j);
		}
	}
}


std::complex<double> EvolutionKernels::kernelElement(const int npRow, const int npCol, const ActionAngleBasisContainer & basisFunc, const double time) const
{
	Eigen::MatrixXcd integrand(basisFunc.size(0), basisFunc.size(1));
	std::complex<double> unitComplex(0,1);	
	for (int i = 0; i < integrand.rows(); ++i) 
	{
		for (int j = 1; j <= i; ++j)
		{
			integrand(i,j) = 0;
			for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++m1) // We do the sum over n1 index
			{
				integrand(i,j) += basisFunc(npRow, m1, i, j) * basisFunc(npCol, m1, i, j) * (m_fourierHarmonic*m_dfdJGrid(i,j) + m1*m_dfdEGrid(i,j))
				* exp(unitComplex*time * ((m1 * m_om1Grid(i,j) + m_fourierHarmonic * m_om2Grid(i,j))));	
			}
		}
	}

	 return integration2d(integrand); 
}


std::complex<double> EvolutionKernels::integration2d(const Eigen::MatrixXcd & grid2Integrate) const 
{	
	std::complex<double> integral{0}, rowIntegral{0};
	for (int i = 0; i < grid2Integrate.rows(); ++i)
	{
		for (int j = 1; j < i; ++j)
		{
			rowIntegral += m_spacing * grid2Integrate(i,j) * (1/m_om1Grid(i,j))*m_elJacobian(i * m_spacing, j * m_spacing);
		}
		integral +=  m_spacing * rowIntegral;
		rowIntegral = 0;
	}
	return integral;
}

void EvolutionKernels::kernelWrite2File() const
{
	//can we write some funciton to save the kernel please
}



#endif




