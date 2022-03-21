#ifndef NBODYBAR
#define NBODYBAR

#include <Eigen/Dense>

#include "NBody.h"

#include "../Perturbation_Grids/BarGrid.h"
#include "../Orbit_Sections/OrbitSections.h"

template <class Tbf>
class NBodyBar : public NBody<Tbf>
{
public:
	NBodyBar(const int nParticles, const int numbTimeSteps, const double timesStep, const Tbf & bf, const Bar2D & bar, const double sigma = 0.35) : 
	NBody<Tbf>(nParticles, numbTimeSteps, timesStep, bf, 0.4, sigma),
	m_barGrid(bf, 20, 300, bar)
	{}
	~NBodyBar() {}

	void testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating); 
	void nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating);

	void orbitSectionPerturbation(const std::string & filename, const bool isSelfConsistent) {this->orbitSections(filename, isSelfConsistent, m_barGrid);}
	void angularMomentumSections(const std::string & filename, const bool isSelfConsistent); 
	void countTrappedOrbits(double runTime = 750); 

private:
	BarGrid m_barGrid; 

	void barUpdate(const double time, const double freelyRotating, const int updateBarEvery);
};

template <class Tbf>
void NBodyBar<Tbf>::testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100}; m_barGrid.updateGrid(0);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		this->outputInfo(time, out);
	    this->backgroundParticleEvolution(false);
	    
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(time*this->m_timeStep, freelyRotating, updateBarEvery); }
		this->foregroundParticleEvolution(false, m_barGrid);
	}
	m_barGrid.saveBarEvolution(barFile);
	out.close();
}

template <class Tbf>
void NBodyBar<Tbf>::nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100}; m_barGrid.updateGrid(0);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		this->outputInfo(time, out);
	    this->backgroundParticleEvolution(true);
	    
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(time*this->m_timeStep, freelyRotating, updateBarEvery); }
		this->foregroundParticleEvolution(true, m_barGrid);

	}
	out.close();
	m_barGrid.saveBarEvolution(barFile);
}

template <class Tbf>
void NBodyBar<Tbf>::barUpdate(const double time, const double freelyRotating, const int updateBarEvery){
	m_barGrid.driftBar(updateBarEvery * this->m_timeStep, freelyRotating);
	m_barGrid.updateGrid(time);

	Eigen::VectorXcd coef = this->m_foreground.responseCoefficents(this->m_basisFunction) - this->m_background.responseCoefficents(this->m_basisFunction);
	m_barGrid.kickBar(updateBarEvery * this->m_timeStep, coef, freelyRotating, time);	
}


template <class Tbf>
void NBodyBar<Tbf>::angularMomentumSections(const std::string & filename, const bool isSelfConsistent) 
{
	OrbitSections sectionsClass(20); int minIndex{0}, index{0}, skip{100}; m_barGrid.updateGrid(); 
	do 
	{	
		m_barGrid.driftBar(this->m_timeStep, 0);
		if (m_barGrid.updateGridNow(index, skip)) {m_barGrid.updateGrid();} 

		sectionsClass.driftStep(this->m_timeStep);
		std::valarray<double> accels = (this->accelsFromBackground(sectionsClass.m_ptle)+m_barGrid.perturbationAccels(sectionsClass.m_ptle));
		if (isSelfConsistent) {accels += this->accelsFromDisk(sectionsClass.m_ptle);}

		sectionsClass.m_ptle.vxvy += accels * this->m_timeStep;
		sectionsClass.driftStep(this->m_timeStep);

		sectionsClass.angularMomentumSections(m_barGrid.angle());

		if (minIndex != sectionsClass.minIndex()) {minIndex = sectionsClass.minIndex(); std::cout << "Min Index: " << minIndex << '\n';}  	

 		index +=1;

	} while (sectionsClass.continueSections()); 
	std::cout << index * this->m_timeStep;
	sectionsClass.outputSections(filename);
}

template <class Tbf>
void NBodyBar<Tbf>::countTrappedOrbits(double runTime) {
	OrbitSections sectionsClass(22); int minIndex{0}, index{0}, skip{100}; m_barGrid.updateGrid(); 
	sectionsClass.setUpforCounting(); 

	for (double time  = 0; time < runTime; time += this->m_timeStep) {
		m_barGrid.driftBar(this->m_timeStep, 0);
		if (m_barGrid.updateGridNow(index, skip)) {m_barGrid.updateGrid();} 

		sectionsClass.driftStep(this->m_timeStep);
		std::valarray<double> accels = (this->accelsFromBackground(sectionsClass.m_ptle)+m_barGrid.perturbationAccels(sectionsClass.m_ptle));

		sectionsClass.m_ptle.vxvy += accels * this->m_timeStep;
		sectionsClass.driftStep(this->m_timeStep);

		if (time > 0.25 * runTime) {sectionsClass.countingSections(m_barGrid.angle());}  	
 		index +=1;
 		if (index % 100000 == 0) {std::cout << "Fraction of time completed: " << time / runTime << '\n';}
	} 	
	std::cout << "Fraction of trapped orbits: " << sectionsClass.fractionOfOrbitsTrapped() << '\n'; 
}


#endif 


