#ifndef BARGRID
#define BARGRID

#include <Eigen/Dense>

#include "PerturbationGrid.h"
#include "../../Bar2D/Bar2D.h"

class BarGrid : public PerturbationGrid
{
public:
	template <class T>
	BarGrid(const T & bf, const double rMaxGrid, const int nGrid, const Bar2D & bar) : 
	PerturbationGrid(bf, rMaxGrid, nGrid),
	m_bar{bar}
	{}
	~BarGrid() {}

	bool updateGridNow(const int nBodyIndex, const int totalNsteps) const; 
	void updateGrid(); 
	void updateGrid(const double time); 

	void driftBar(const double timeStep, const double freelyRotating) {m_bar.drift(timeStep, freelyRotating);}
	void kickBar(const double timeStep, const Eigen::VectorXcd &diskCoeff, const double freelyRotating, const double time)
	{m_bar.kick(timeStep, diskCoeff, freelyRotating, time);}

	void saveBarEvolution(const std::string & barFile) {m_bar.saveBarEvolution(barFile);}
	Bar2D m_bar;

private: 

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

void BarGrid::updateGrid(const double time) {
	Eigen::ArrayXXcd potential = Eigen::ArrayXXcd::Zero(m_potentialArray.rows(), m_potentialArray.cols());
	for (int i = 0; i < v_individualPotentials.size(); ++i){
		potential += m_bar.barCoeff(time)[i] * v_individualPotentials[i];
	}
	if (m_fourierHarmonic == 0){m_potentialArray = potential.real();}
	else {m_potentialArray = 2 * potential.real();}
	takeTimeStep();
}


#endif