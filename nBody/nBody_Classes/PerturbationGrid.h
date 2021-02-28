#ifndef PERTURBATIONGRID
#define PERTURBATIONGRID

#include <valarray>

#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../../Volterra_Solver/ExpansionCoeff.h"

#include "Bodies.h"

class PerturbationGrid
{
public:
	PerturbationGrid(const double rMaxGrid, const int nGrid) : m_potentialArray(nGrid, nGrid),
	m_spacing{2*rMaxGrid/ ((double) nGrid-1)}
	// Expansion coeff will go here
	{} // I want this to take the BF and creat the indidual grids

	~PerturbationGrid() {}

	bool updateGridNow(const double nBodyTime);
	
	template <class Tbf>
	void updateGrid(const Tbf & bf);

	std::valarray<double> accels(const Bodies & ptle) const;

private:
		
	Eigen::ArrayXXd m_potentialArray;
	const double m_spacing;
	//Expansion coeff
	


};

#endif