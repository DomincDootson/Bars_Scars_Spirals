#ifndef WAVE 
#define WAVE

#include <cmath>
#include <Eigen/Dense>

#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"

class Wave
{
public:
	Wave(const double omega, const double centre, const double k) :
	m_fourierHarmonic{2}, m_omega{omega}, m_centre{centre}, m_k{k}, m_sigma{0.5}, m_amplitude{0.01}
	{}
	~Wave() {}

	// We need some funciton that can return the initial coefficents. 


	
	const int m_fourierHarmonic;
	const double m_omega, m_centre, m_k, m_sigma, m_amplitude; 
	
	double envelope(const double radius) const;
	std::complex<double> density(const double radius, const double phi) const;

	std::complex<double> resolveCoefficent(const int n) const;
	Eigen::VectorXcd resolveIC() const; 
};



double Wave::envelope(const double radius) const {
	double phi = 0.5 * (M_PI/m_sigma) * (radius - m_centre);
	if (abs(phi) < 0.5*M_PI) {return m_amplitude * cos(phi)/(8.0*m_sigma * m_centre);}  
	else {return 0;}
}

std::complex<double> Wave::density(const double radius, const double phi) const {
	std::complex<double> unitComplex(0,1.0);
	return envelope(radius) * exp(unitComplex * (m_fourierHarmonic*phi + m_k*radius)); 
}


#endif



/*
Things to do with this class:

	1. Some way to give it an inital frequency and radius and get k, maybe for the time being just pass both? 
	2. We need a constructor to make the initial conditions , i.e. to get the coefficents

*/ 