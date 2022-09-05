#ifndef SCARREDMESTEL
#define SCARREDMESTEL

#include <vector>

#include "DFClass.h"
#include "Scar.h"
#include "AngularMomentumScar.h"

template <class T> 
class ScarredMestel : public DFClass
{
public:
	ScarredMestel(const std::vector<T> & scars, double vc = 1, double r0 = 1, double litteSigma = 0.35, double xi = 1, 
	double rInner = 1, double rOuter = 11.5, double nuTaper = 4, double muTaper = 5) : 
	m_vc{vc}, m_r0{r0}, m_littleSigma{litteSigma}, m_q{pow(m_vc/m_littleSigma,2) - 1}, m_xi{xi}, 
	m_rInner{rInner}, m_rOuter{rOuter}, m_nuTaper{nuTaper}, m_muTaper{muTaper}, m_massScaling{1} {
		addScars(scars);
		readInPotentialVector();
	}

	ScarredMestel(double vc = 1, double r0 = 1, double litteSigma = 0.35, double xi = 1, 
	double rInner = 1, double rOuter = 11.5, double nuTaper = 4, double muTaper = 5) : 
	m_vc{vc}, m_r0{r0}, m_littleSigma{litteSigma}, m_q{pow(m_vc/m_littleSigma,2) - 1}, m_xi{xi}, 
	m_rInner{rInner}, m_rOuter{rOuter}, m_nuTaper{nuTaper}, m_muTaper{muTaper}, m_massScaling{1} {
		readInPotentialVector();
	}
	
	
	~ScarredMestel() {}

	double potentialBackground(double radius) const {return m_vc*m_vc*log(radius/m_r0);}
	double potential(double radius) const {return m_vc*m_vc*log(radius/m_r0);}//(1-m_xi) * potentialBackground(radius) + m_xi * potentialFromVector(radius);}//(1-m_xi) * potential(radius) + 
	double distFunc(double E, double J) const;

	void setDiskMass(const double diskMassValue) {m_massScaling = diskMassValue/diskMass();}

	template <class Ts> 
	void addScar(const Ts & scarToAdd) {v_scars.push_back(scarToAdd);}
	template <class Ts> 
	void addScars(const std::vector<Ts> & scars) {for (auto scar : scars) {v_scars.push_back(scar);}}

	double vRSampling() const;  
	double potentialFromVector(const double radius) const; 
	
private:
	std::vector<T> v_scars; 
	std::vector<double> v_potential; 

	const double m_vc, m_r0, m_littleSigma, m_q, m_xi, m_rInner, m_rOuter, m_nuTaper, m_muTaper;

	double m_massScaling, m_potentialSpacing; 

	double innerTaper(double J) const;
	double outerTaper(double J) const;
	
	double radialAction(const double E, const double L) const {return (exp(E/(m_vc*m_vc)) * m_vc * m_r0 - sqrt(M_E) * L) / sqrt(2 * M_PI);}
	double equivalentLconstJr(const double E, const double L, const double pattern, const double jacobi) const;
	double toRootFind (double L, double jR, double pattern, double jacobi) const {return jacobi - pattern * L - m_vc*m_vc * log((sqrt(M_E) * L + sqrt(2*M_PI)*jR)/(m_r0*m_vc));}

	double jMax(const double radius) const;
	double dfMax(const double radius, const double vR) const;
	

	double vRScale() const {return 3*m_littleSigma;}; // Give the upper limit for marganlising over
	double vPhiScale() const {return m_vc + 3*m_littleSigma;}

	void readInPotentialVector();
};

/* DF Functions */
/* ------------ */  

template <class T>
double ScarredMestel<T>::innerTaper(double J) const{
	return pow(J, m_nuTaper) / (pow(m_rInner*m_vc, m_nuTaper) + pow(J, m_nuTaper));
}

template <class T>
double ScarredMestel<T>::outerTaper(double J) const{
	return 1/(1+pow(J/(m_rOuter*m_vc), m_muTaper));
}

template <class T>
double ScarredMestel<T>::distFunc(double E, double J) const // Note that we have set G = 1
{
	double scarEffect{1};
	for (auto scar : v_scars) {scarEffect *= scar.groove(E, J, 1);} // We need the other way for the De Rijke Scar scar! 

	//for (auto scar : v_scars) {scarEffect *= scar.groove(E, J, equivalentLconstJr(E, J, scar.patternSpeed(), scar.jacobi()));}
	return  scarEffect*innerTaper(J)*outerTaper(J) *
	(m_xi)*pow((J/(m_r0*m_vc)), m_q) * exp(-E/pow(m_littleSigma,2)) * 
	(((pow(m_vc,2)/(2*M_PI*m_r0)) * pow(m_vc,m_q)) * pow(pow(2, 0.5*m_q) * 
	sqrt(M_PI)*tgamma(0.5*m_q+0.5)*pow(m_littleSigma,2+m_q), -1)); 
}

/* Scarring Functions */ 
/* ------------------ */

template <class T>
double ScarredMestel<T>::equivalentLconstJr(const double E, const double L, const double pattern, const double jacobi) const {
	double jR{radialAction(E,L)}, epsilon{0.000001}; int nStep{20};
	double lowerL{epsilon}, upperL{50 * m_vc}, holding{0.5*(lowerL + upperL)}; 

	for (int i = 0; i < nStep; ++i) {
		if (toRootFind(lowerL, jR, pattern, jacobi) * toRootFind(holding, jR, pattern, jacobi) <= 0) {upperL = holding; holding = 0.5 * (lowerL + upperL);}
		else {lowerL = holding; holding = 0.5 * (lowerL + upperL);}
	}
	return holding; 
}


/* Potential Functions */
/* ------------------- */ 

template <class T>
void ScarredMestel<T>::readInPotentialVector() {
	//savingDensity("DF_Class/Densities/scarredMestelDensity.csv", 0.03, 1000);	
	//std::system("/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/nBody/potentialFromDensity"); 
	std::ifstream inFile("DF_Class/Densities/scarredMestelPotential.csv");
	
	int nPoint; double holding; inFile >> nPoint >> m_potentialSpacing; 

	v_potential.resize(nPoint);
	for (int i = 0; i < v_potential.size(); ++i) {
		inFile >> holding; 
		v_potential[i] = holding;}
	inFile.close(); 
}

template <class T>
double ScarredMestel<T>::potentialFromVector(const double radius) const {
	if (v_potential.size() ==0) {return 0;}
	if (radius > (v_potential.size()-2)*m_potentialSpacing) {std::cout << "Off potential grid\n"; exit(0);}
	double index{radius/m_potentialSpacing}; int lowerIndex{(int) floor(index)};
	return v_potential[lowerIndex] + (index - lowerIndex) * (v_potential[lowerIndex+1] - v_potential[lowerIndex]); 
}


/* Sampling Functions */
/* ------------------ */ 

template <class T>
double ScarredMestel<T>::vRSampling() const {
	std::random_device generator; std::normal_distribution<double> vrDF(0, m_littleSigma);
	return vrDF(generator);
}


template <class T>
double ScarredMestel<T>::jMax(const double radius) const{

	if (radius > m_rOuter) {return 5 * m_rOuter * m_vc;}
	else {return sqrt(2) * 5 * m_littleSigma * radius;}
}

template <class T>
double ScarredMestel<T>::dfMax(const double radius, const double vR) const {
	return distFunc(0.5*vR*vR +0.5*m_q*pow(m_littleSigma,2)+potential(radius), sqrt(m_q)*m_littleSigma*radius);
}


#endif
