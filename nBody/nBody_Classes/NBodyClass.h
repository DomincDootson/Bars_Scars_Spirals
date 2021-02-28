#ifndef NBODYCLASS
#define NBODYCLASS

#include <Eigen/Dense>

#include "Bodies.h"
#include "Box.h"
#include "PerturbationGrid.h"

#include "../../DF_Class/Mestel.h"
#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

template <class Tbf>
class NBodyClass
{
public:
	NBodyClass(const int nParticles, const int numbTimeSteps, const int fourierHamonic, const double timesStep, const Tbf & bf) : 
	 m_foreground(nParticles), m_background(nParticles),
	 m_foregroundBox(120, 26.0, 0.18,1), m_backgroundBox(120, 26.0, 0.18,1),
	 m_pertGrid(20,400),
	 m_DF(),
	 m_basisFunction(bf),
	 m_numbTimeSteps{numbTimeSteps}, m_fourierHarmonic{fourierHamonic}, m_timeStep{timesStep}
	 {}

	~NBodyClass() {}

private:
	Bodies m_foreground, m_background;
	Box m_foregroundBox, m_backgroundBox;

	PerturbationGrid m_pertGrid;
	const Mestel m_DF;

	Tbf m_basisFunction;

	// some params 

	const int m_numbTimeSteps, m_fourierHarmonic;
	const double m_timeStep;
};


#endif