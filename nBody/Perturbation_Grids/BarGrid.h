#ifndef BARGRID
#define BARGRID

#include <Eigen/Dense>

#include "PerturbationGrid.h"
#include "../Bodies/Bodies.h"

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

	void  sectionPlotterSampler(Bodies & ptle, const double rUpper) const;
	double angle() const {return m_bar.angle();}
	double patternSpeed() const {return m_bar.patternSpeed();}

	double corotatingPx(const double xdot, const double y) const {return xdot + m_bar.patternSpeed() * y;}
	double corotatingPy(const double ydot, const double x) const {return ydot - m_bar.patternSpeed() * x;}

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




void BarGrid::sectionPlotterSampler(Bodies & ptle, const double rUpper) const
{
	double barAngle{m_bar.angle()}, barSpeed{m_bar.patternSpeed()}, rLower{1.5};
    double xUpper{rUpper * cos(barAngle)}, yUpper{rUpper * sin(barAngle)};
    
    double energy = log(rUpper) + potential(xUpper, yUpper); // Sample along the axis of the bar
    std::cout << "Jacobi " << energy << " " << " " << barSpeed << " " << .5*pow(rUpper*barSpeed, 2)<< '\n';
    
    for (int i =0; i < ptle.n;++i) {  
       ptle.xy[2*i]   = rLower + (rUpper - rLower) * (i/((double) ptle.n)) ; ptle.xy[2*i+1] = barAngle;
       // calculate the PE (background + perturbation)
       double PE = log(ptle.xy[2*i]) + potential(ptle.xy[2*i] * cos(barAngle), ptle.xy[2*i] * sin(barAngle));
       ptle.vxvy[2*i] = 0; ptle.vxvy[2*i+1] = sqrt(2 * (energy - PE));
    }
    ; // Do sampling in polar coords. 
    ptle.convert2Cartesian(); 
}

#endif