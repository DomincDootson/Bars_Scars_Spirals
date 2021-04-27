#ifndef NBODYPERTURBATION
#define NBODYPERTURBATION

#include "NBody.h"
#include "../Perturbation_Grids/ChoosenPerturbationGrid.h"

template <class Tbf>
class NBodyPerturbation : public NBody<Tbf>
{
public:
	NBodyPerturbation(const int nParticles, const int numbTimeSteps, const double timesStep, const Tbf & bf) : 
	NBody<Tbf>(nParticles, numbTimeSteps, timesStep, bf),
	m_pertGrid(bf, 20,400, "/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/nBody/Perturbation_0.out", 2000)
	{}

	~NBodyPerturbation() {}

	void testParticleEvolution(const std::string & filename); 
	void nBodyEvolution(const std::string & filename);

private: 
	ChoosenPerturbationGrid m_pertGrid;	
};

template <class Tbf>
void NBodyPerturbation<Tbf>::testParticleEvolution(const std::string & filename){ 
	std::ofstream out(filename);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (time % this->m_skip == 0) {this->outputCoefficents(out); std::cout << "Fraction of test particle: " << time/((double) this->m_numbTimeSteps) << '\n';}
	    this->backgroundParticleEvolution(false);
		if (m_pertGrid.updateGridNow(time, this->m_numbTimeSteps)) {
			m_pertGrid.updateGrid();}
		this->foregroundParticleEvolution(false, m_pertGrid);
	}
	out.close();
	m_pertGrid.saveArray("pertGrid.csv");
}

template <class Tbf>
void NBodyPerturbation<Tbf>::nBodyEvolution(const std::string & filename){ 
	std::ofstream out(filename);
	for (int time = 0; time < this->m_numbTimeSteps; ++time){
		if (time % this->m_skip == 0) {this->outputCoefficents(out); std::cout << "Fraction of NBody: " << time/((double) this->m_numbTimeSteps) << '\n';}
	    this->backgroundParticleEvolution(true);
	    
		if (m_pertGrid.updateGridNow(time,this->m_numbTimeSteps)) {m_pertGrid.updateGrid();}
		this->foregroundParticleEvolution(true, m_pertGrid);	
	}
	out.close();
}

// Because we are inheriting from a templated class, use this-> to tell the complier to inherit it. 

#endif 