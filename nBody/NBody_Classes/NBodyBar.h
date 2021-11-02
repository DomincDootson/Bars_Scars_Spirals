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
	NBodyBar(const int nParticles, const int numbTimeSteps, const double timesStep, const Tbf & bf, const Bar2D & bar) : 
	NBody<Tbf>(nParticles, numbTimeSteps, timesStep, bf),
	m_barGrid(bf, 20, 400, bar)
	{}
	~NBodyBar() {}

	void testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating); 
	void nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating);

	void orbitSectionPerturbation(const std::string & filename, const bool isSelfConsistent) {this->orbitSections(filename, isSelfConsistent, m_barGrid);}
	void barOrbitSections(const std::string & filename, const bool isSelfConsistent); 
	void angularMomentumSections(const std::string & filename, const bool isSelfConsistent); 

	void individualOrbits(int numbStep, double x, double y, double vx, double vy); 

private:
	BarGrid m_barGrid; 

	void barUpdate(const double time, const double freelyRotating, const int updateBarEvery);
};

template <class Tbf>
void NBodyBar<Tbf>::testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100}; m_barGrid.updateGrid(0); std::ofstream barCoef("barCoeffN.csv");
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		this->outputInfo(time, out);
	    this->backgroundParticleEvolution(false);
	    
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(time*this->m_timeStep, freelyRotating, updateBarEvery); 
		barCoef << (m_barGrid.m_bar).barCoeff(time*this->m_timeStep)(0).real() << ',' << (m_barGrid.m_bar).barCoeff(time*this->m_timeStep)(0).imag() << '\n';
		}
		this->foregroundParticleEvolution(false, m_barGrid);
	}
	out.close();
	m_barGrid.saveBarEvolution(barFile);
}

template <class Tbf>
void NBodyBar<Tbf>::nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100}; m_barGrid.updateGrid(0);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (time % this->m_skip == 0) {
			this->outputCoefficents(out); std::cout << "Fraction of NBody: " << time/((double) this->m_numbTimeSteps) << '\n';
		}
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

	Eigen::VectorXcd coef = this->m_foreground.responseCoefficents(this->m_basisFunction) 
	- this->m_background.responseCoefficents(this->m_basisFunction);
	m_barGrid.kickBar(updateBarEvery * this->m_timeStep, coef, freelyRotating, time);	
}




template <class Tbf>
void NBodyBar<Tbf>::barOrbitSections(const std::string & filename, const bool isSelfConsistent) // put in drift step for the bar
{
	OrbitSections sectionsClass(m_barGrid, 15); int minIndex{0}, index{0}, skip{100}; m_barGrid.updateGrid(); 
	do 
	{	
		m_barGrid.driftBar(this->m_timeStep, 0);
		if (m_barGrid.updateGridNow(index, skip)) {m_barGrid.updateGrid();}

		sectionsClass.driftStep(this->m_timeStep);
		std::valarray<double> accels = (this->accelsFromBackground(sectionsClass.m_ptle)+m_barGrid.perturbationAccels(sectionsClass.m_ptle));
		if (isSelfConsistent) {accels += this->accelsFromDisk(sectionsClass.m_ptle);}

		sectionsClass.m_ptle.vxvy += accels * this->m_timeStep;
		sectionsClass.driftStep(this->m_timeStep);


		if (minIndex != sectionsClass.minIndex()) {minIndex = sectionsClass.minIndex(); std::cout << "Min Index: " << minIndex << '\n';}  	
		 	
	 	index +=1;	


	} while (sectionsClass.continueSections(m_barGrid));
	
	sectionsClass.outputSections(filename); 
}

template <class Tbf>
void NBodyBar<Tbf>::angularMomentumSections(const std::string & filename, const bool isSelfConsistent) 
{
	OrbitSections sectionsClass(m_barGrid, 15); int minIndex{0}, index{0}, skip{100}; m_barGrid.updateGrid(); 
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

	} while (sectionsClass.continueSections()); //sectionsClass.continueSections()
	sectionsClass.outputSections(filename); 
}


template <class Tbf> 
void NBodyBar<Tbf>::individualOrbits(int numbStep, double x, double y, double vx, double vy) 
{
	std::ofstream out("particleEvolution.csv"); int skip{100}; m_barGrid.updateGrid();
	Bodies ptle0(x, y, vx, vy), ptle1(x, y, vx, vy);
	out << ptle0.xy[0] << ',' << ptle0.xy[1] << ',' << ptle1.xy[0] << ','<< ptle1.xy[1] << '\n';
	
	

	// Sample two particles
	for (int step = 0; step < numbStep; ++step){
		m_barGrid.driftBar(this->m_timeStep, 0);
		if (m_barGrid.updateGridNow(step,skip)) {m_barGrid.updateGrid();}
        
        ptle0.xy += 0.5*ptle0.vxvy*this->m_timeStep; ptle1.xy += 0.5*ptle1.vxvy*this->m_timeStep;
		ptle0.vxvy += (this->accelsFromBackground(ptle0))*this->m_timeStep; ptle1.vxvy += (this->accelsFromBackground(ptle1)+m_barGrid.perturbationAccels(ptle1))*this->m_timeStep;
		ptle0.xy += 0.5*ptle0.vxvy*this->m_timeStep; ptle1.xy += 0.5*ptle1.vxvy*this->m_timeStep;
 		

		//std::cout << m_barGrid.perturbationAccels(ptle1)[0] << '\n';
 		//std::cout << std::abs(m_barGrid.m_bar.barCoeff()[29]) << '\n';
		
			
	}
	out.close(); 
} 


#endif 


