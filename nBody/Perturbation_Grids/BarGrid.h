#ifndef BARGRID
#define BARGRID

#include <Eigen/Dense>

#include "PerturbationGrid.h"
#include "../../Bar2D/Bar2D.h"

class BarGrid : public PerturbationGrid
{
public:
	template <class T>
	BarGrid(const T & bf, const double rMaxGrid, const int nGrid, const Eigen::VectorXcd &expansionCoeff, const double omega0) : 
	PerturbationGrid(bf, rMaxGrid, nGrid),
	m_bar(expansionCoeff, bf, omega0)
	{updateGrid();}
	~BarGrid() {}

	bool updateGridNow(const int nBodyIndex, const int totalNsteps) const; 
	void updateGrid(); 

	void driftBar(const double timeStep) {m_bar.drift(timeStep);}
	void kickBar(const double timeStep, const Eigen::VectorXcd &diskCoeff, const double freelyRotating)
	{m_bar.kick(timeStep, diskCoeff, freelyRotating);}

	void saveBarEvolution(const std::string & barFile) {m_bar.saveBarEvolution(barFile);}


private: 
	Bar2D m_bar;
};


bool BarGrid::updateGridNow(const int nBodyIndex, const int totalNsteps) const {  // update every totalNsteps
	if ((nBodyIndex % totalNsteps == 0) && nBodyIndex != 0){
		return true; 
	}
	return false;
}

void BarGrid::updateGrid() {
	Eigen::ArrayXXcd potential = Eigen::ArrayXXcd::Zero(m_potentialArray.rows(), m_potentialArray.cols());
	for (int i = 0; i < v_individualPotentials.size(); ++i){
		potential += m_bar.barCoeff()[i] * v_individualPotentials[i];
	}
	if (m_fourierHarmonic == 0){m_potentialArray = potential.real();}
	else {m_potentialArray = 2 * potential.real();}
	takeTimeStep();
}


#endif