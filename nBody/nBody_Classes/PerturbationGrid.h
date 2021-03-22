#ifndef PERTURBATIONGRID
#define PERTURBATIONGRID

#include <valarray>

#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../../Volterra_Solver/ExpansionCoeff.h"

#include "Bodies.h"

// I think a few checks wouldn't go a miss, maybe checking to see if l =0 perturbation conserves AM?? 

class PerturbationGrid
{
public:
	PerturbationGrid(const double rMaxGrid, const int nGrid, const std::string & filename, const int numbTimeStep) : m_potentialArray(nGrid, nGrid),
	m_spacing{2*rMaxGrid/ ((double) nGrid-1)}, m_rMaxGrid{rMaxGrid}, m_centre{0.5 * (nGrid-1)},
	m_perturbationCoeff(filename, numbTimeStep),
	m_perturbationIndex{0}
	// Expansion coeff will go here
	{} // I want this to take the BF and creat the indidual grids

	~PerturbationGrid() {}

	bool updateGridNow(const int nBodyIndex, const int totalNsteps) const;
	
	template <class Tbf>
	void updateGrid(const Tbf & bf);

	std::valarray<double> perturbationAccels(const Bodies & ptle) const;

private:
		
	Eigen::ArrayXXd m_potentialArray;
	const double m_spacing, m_rMaxGrid, m_centre;
	
	ExpansionCoeff m_perturbationCoeff; 
	double m_perturbationIndex;
	// Some internal clock parameters

	bool checkOnGrid(const double xPosition, const double yPosition) const;

	double xPos2Index(const double xPosition) const {return xPosition/m_spacing + m_centre;} // 0 < i,j < m_numbgridspace-1 !! Note not strictly
	double yPos2Index(const double yPosition) const {return yPosition/m_spacing + m_centre;} 

	double xAccel(const int i, const int j) const {return -(1/(2*m_spacing)) * (m_potentialArray(i,j+1) - m_potentialArray(i,j-1));} 
	double xAccel(const double xPosition, const double yPosition) const;

	double yAccel(const int i, const int j) const {return -(1/(2*m_spacing)) * (m_potentialArray(i+1,j) - m_potentialArray(i+1,j));} 
	double yAccel(const double xPosition, const double yPosition) const;

	void takeTimeStep() {m_perturbationIndex += 1;}


};


// Functions that get accels //

int floorInt(const double index) {return index/1;} // Overload to make sure that an integer is returned
int ceilInt(const double index) {return (index+1)/1;}

bool PerturbationGrid::checkOnGrid(const double xPosition, const double yPosition) const{
	if ((xPosition*xPosition + yPosition*yPosition) < m_rMaxGrid*m_rMaxGrid){
		return true;}
	return false;
}


double PerturbationGrid::xAccel(const double xPosition, const double yPosition) const 
{
	if (!checkOnGrid(xPosition, yPosition)) {return 0;}
	double i{xPos2Index(xPosition)}, j{yPos2Index(yPosition)};
	int i0{floorInt(i)}, j0{floorInt(j)}, i1{ceilInt(i)}, j1{ceilInt(j)};
	double accel_i0_j0{xAccel(i0, j0)}, accel_i1_j0{xAccel(i1, j0)}, accel_i0_j1{xAccel(i0, j1)}, accel_i1_j1{xAccel(i1, j1)};

	return accel_i0_j0*(1-i)*(1-j) + accel_i1_j0*i*(1-j0) + accel_i0_j1*(1-i)*j + accel_i1_j1*i*j; //linear interpolation 
}

double PerturbationGrid::yAccel(const double xPosition, const double yPosition) const 
{
	if (!checkOnGrid(xPosition, yPosition)) {return 0;}
	double i{xPos2Index(xPosition)}, j{yPos2Index(yPosition)};
	int i0{floorInt(i)}, j0{floorInt(j)}, i1{ceilInt(i)}, j1{ceilInt(j)};
	double accel_i0_j0{yAccel(i0, j0)}, accel_i1_j0{yAccel(i1, j0)}, accel_i0_j1{yAccel(i0, j1)}, accel_i1_j1{yAccel(i1, j1)};

	return accel_i0_j0*(1-i)*(1-j) + accel_i1_j0*i*(1-j0) + accel_i0_j1*(1-i)*j + accel_i1_j1*i*j; //linear interpolation 
}


std::valarray<double> PerturbationGrid::perturbationAccels(const Bodies & ptle) const 
{
	std::valarray<double> accels(2*ptle.n);
	for (int nParticle = 0; nParticle < ptle.n; nParticle +=2){
		accels[nParticle]   = xAccel(ptle.xy[nParticle], ptle.xy[nParticle+1]);
		accels[nParticle+1] = yAccel(ptle.xy[nParticle], ptle.xy[nParticle+1]);
	}
	return accels;
}

// Functions that do the updating of grid  --> maybe add an 'internal clock' so it knows when to update itself 

bool PerturbationGrid::updateGridNow(const int nBodyIndex, const int totalNsteps) const {
	if (nBodyIndex % (totalNsteps/m_perturbationCoeff.nTimeStep()) == 0){
		return true;
	}
	return false;
}

template <class Tbf>
void PerturbationGrid::updateGrid(const Tbf & bf) {
	m_potentialArray = bf.potentialArrayReal(m_perturbationCoeff(m_perturbationIndex), m_potentialArray.cols(), m_rMaxGrid);
	takeTimeStep();
}


#endif