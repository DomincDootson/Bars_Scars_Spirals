#ifndef NBODYBAR
#define NBODYBAR

#include <Eigen/Dense>

#include "NBody.h"
#include "../Perturbation_Grids/BarGrid.h"

template <class Tbf>
class NBodyBar : public NBody<Tbf>
{
public:
	NBodyBar(const int nParticles, const int numbTimeSteps, const double timesStep, const Tbf & bf, const Eigen::VectorXcd &expansionCoeff, const double omega0) : 
	NBody<Tbf>(nParticles, numbTimeSteps, timesStep, bf),
	m_barGrid(bf, 20, 400, expansionCoeff, omega0)
	{}
	~NBodyBar() {}

	void testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating); 
	void nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating);

private:
	BarGrid m_barGrid; 

	void barUpdate(const double freelyRotating, const int updateBarEvery);
};

template <class Tbf>
void NBodyBar<Tbf>::testParticleEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100};
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (time % this->m_skip == 0) {
			this->outputCoefficents(out); std::cout << "Fraction of Test particle: " << time/((double) this->m_numbTimeSteps) << '\n';
		}
	    this->backgroundParticleEvolution(false);
	    
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(freelyRotating, updateBarEvery); }
		this->foregroundParticleEvolution(false, m_barGrid);
	}
	out.close();
	m_barGrid.saveBarEvolution(barFile);
}

template <class Tbf>
void NBodyBar<Tbf>::nBodyEvolution(const std::string & diskFile, const std::string & barFile, const double freelyRotating){ 
	std::ofstream out(diskFile); int updateBarEvery{100};
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (time % this->m_skip == 0) {
			this->outputCoefficents(out); std::cout << "Fraction of NBody: " << time/((double) this->m_numbTimeSteps) << '\n';
		}
	    this->backgroundParticleEvolution(true);
	    
		if (m_barGrid.updateGridNow(time,updateBarEvery)) { barUpdate(freelyRotating, updateBarEvery); }
		this->foregroundParticleEvolution(true, m_barGrid);

	}
	out.close();
	m_barGrid.saveBarEvolution(barFile);
}

template <class Tbf>
void NBodyBar<Tbf>::barUpdate(const double freelyRotating, const int updateBarEvery){
	m_barGrid.driftBar(updateBarEvery * this->m_timeStep);
	m_barGrid.updateGrid();

	Eigen::VectorXcd coef = this->m_foreground.responseCoefficents(this->m_basisFunction) 
	- this->m_background.responseCoefficents(this->m_basisFunction);
	m_barGrid.kickBar(updateBarEvery * this->m_timeStep, coef, freelyRotating);	
}
#endif 