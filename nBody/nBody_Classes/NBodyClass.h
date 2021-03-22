#ifndef NBODYCLASS
#define NBODYCLASS

#include <Eigen/Dense>
#include <complex>
#include <string>

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
	 m_DF(),
	 m_foreground("/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/nBody/cumulativeDensity.csv", m_DF, nParticles), 
	 m_background("/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/nBody/cumulativeDensity.csv", m_DF, nParticles), // We dont neet to sample
	 m_diskBox(120, 26.0, 0.18,1), m_m0Box(120, 26.0, 0.18,1),
	 m_pertGrid(20,400, "/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/nBody/Perturbation_5_0.out", 2000),
	 m_basisFunction(bf),
	 m_numbTimeSteps{numbTimeSteps}, m_fourierHarmonic{fourierHamonic}, m_skip{100}, m_timeStep{timesStep}, m_xi{.5}
	 {}

	~NBodyClass() {}

	void cumulativeDistribution(const std::string & filename) const {m_DF.cumulativeDensity(filename);}

	void nBodyEvolution(const std::string & filename);

	void particleSampling() {m_foreground.samplingDF();}

private:
	Bodies m_foreground, m_background;
	Box m_diskBox, m_m0Box;

	PerturbationGrid m_pertGrid;
	const Mestel m_DF;

	Tbf m_basisFunction;

	const int m_numbTimeSteps, m_fourierHarmonic, m_skip;
	const double m_timeStep, m_xi;

	void outputCoefficents(std::ofstream & out);

	void backgroundParticleEvolution(); 
	void foregroundParticleEvolution(); 
	valarray<double>  accelsFromDisk(const Bodies & ptle);

	void m0Grid();
};






template <class Tbf>
void NBodyClass<Tbf>::nBodyEvolution(const std::string & filename){
	std::ofstream out(filename);
	if (m_fourierHarmonic == 0){m0Grid();}
	
	for (int time = 0; time < m_numbTimeSteps; ++time){
	if (time % m_skip == 0) {outputCoefficents(out); std::cout << "Fraction of NBody: " << time/((double) m_numbTimeSteps) << '\n';}
	
	backgroundParticleEvolution();
	
	if (m_pertGrid.updateGridNow(time, m_numbTimeSteps)) {m_pertGrid.updateGrid(m_basisFunction);}
	foregroundParticleEvolution();
	}
	out.close();
}


// General N-Body Functions //
// ------------------------ //  

std::string outComplexNumber(const std::complex<double> number){
	std::string outString = std::to_string(real(number));
	if (imag(number) < 0){return outString + '-' +std::to_string(abs(imag(number)));} 
	else {return outString + '+' + std::to_string(imag(number));}
}

template <class Tbf>
void NBodyClass<Tbf>::outputCoefficents(std::ofstream & out){
	Eigen::VectorXcd coef = m_foreground.responseCoefficents(m_basisFunction, m_xi) - m_background.responseCoefficents(m_basisFunction, m_xi);
	for (int i =0; i < coef.size()-1; ++i){
		out << outComplexNumber(coef(i)) <<',';
	}
	out << outComplexNumber(coef(coef.size()-1)) << '\n';
}

template <class Tbf>
valarray<double>  NBodyClass<Tbf>::accelsFromDisk(const Bodies & ptle){
	m_diskBox.zero();
    m_diskBox.bodies2density_m2(ptle, m_fourierHarmonic, 2400, 720); // Put in values for nPhi and nRing
    m_diskBox.density2pot();
    valarray<double> accels = m_diskBox.pot2accels(ptle);
    
    if (m_fourierHarmonic ==0){ accels -= m_m0Box.pot2accels(ptle);}
    
    for(int i=0;i<ptle.n;i++) 
    {
    	double x{ptle.xy[2*i]}, y{ptle.xy[2*i+1]}, R2{x*x+y*y}, ax{accels[2*i]}, ay{accels[2*i+1]}; 
    	ax += m_DF.xAccel(x, y);    	
     	ay += m_DF.yAccel(x, y);

    	accels[2*i] = ax;
    	accels[2*i+1] = ay;
    }
    return accels; 
}
template <class Tbf>
void NBodyClass<Tbf>::backgroundParticleEvolution() {
	m_background.xy   += m_background.vxvy * m_timeStep * 0.5;
	m_background.vxvy += accelsFromDisk(m_background) * m_timeStep;
	m_background.xy   += m_background.vxvy * m_timeStep * 0.5;	
}
template <class Tbf>
void NBodyClass<Tbf>::foregroundParticleEvolution() {
	m_foreground.xy   += m_foreground.vxvy * m_timeStep * 0.5;
	m_foreground.vxvy += (accelsFromDisk(m_foreground) + perturbationAccels(m_foreground)) * m_timeStep;
	m_foreground.xy   += m_foreground.vxvy * m_timeStep * 0.5;	
}
template <class Tbf>
void NBodyClass<Tbf>::m0Grid() {
	m_m0Box.zero();
	m_m0Box.bodies2density_m2(m_background, 0, 2400, 720);
	m_m0Box.density2pot();
}


#endif