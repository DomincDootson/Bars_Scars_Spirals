#ifndef SPIRALBASIS
#define SPIRALBASIS

#include <complex>
#include <boost/math/special_functions/bessel.hpp>
#include "PotentialDensityPair.h"

class SpiralBasis : public PotentialDensityPair 
{
public:
	SpiralBasis(std::vector<double> params, int m_radialIndex, int m_fourierHarmonic)
		: PotentialDensityPair(m_radialIndex, m_fourierHarmonic), 
		m_r0{params[0] * (m_radialIndex+1)}, m_sigma{params[1]}, m_kR{params[2]}
		{} 

	~SpiralBasis() {} 


	double density(double r) const; 
	double potential(double r) const; 

	std::complex<double> analyticKernel(const double timeDelta, const double sigmaR) const;
	double scriptWElementAnalytic(const int m1, const double jA, const double rG) const; 
	

private:
	const double m_r0, m_sigma, m_kR;

	double circularFreq() const {return 1/m_r0;} // We've assumed v_{c} = 1 so J = m_r0 
	double epicycleFreq() const {return circularFreq() * sqrt(2);}
	double surfaceDensity() const {return (1/(2*M_PI*m_r0));}

	double chi(const double sigmaR) const {return pow((sigmaR * m_kR)/epicycleFreq(), 2);}
	std::complex<double> analyticKernel(const double mR, const double timeDelta, const double sigmaR) const;
};


double SpiralBasis::potential(double r) const { 
	return sqrt(2/(m_kR*m_r0)) * cos(m_kR * r) * pow(M_PI * m_sigma*m_sigma, -0.25)  // We include sqrt(2) so that orthonormal
	* exp(- 0.5 * pow((r-m_r0)/m_sigma,2)); 
}

double SpiralBasis::density(double r) const { 
	return - m_kR * (1/(2*M_PI)) * potential(r); 
}


/* Calculating the Analytic Kernel */
/* ------------------------------- */ 

std::complex<double> SpiralBasis::analyticKernel(const double timeDelta, const double sigmaR) const {
	int upperMr{8}; std::complex<double> sum{0};
	for (int mr = 1; mr <upperMr; ++mr) {sum += analyticKernel(mr, timeDelta, sigmaR);}// std::cout << mr << " " << sum << '\n';
	
	return sum; 
}

std::complex<double> SpiralBasis::analyticKernel(const double mR, const double timeDelta, const double sigmaR) const {
	double x{chi(sigmaR)}, expTimeBessel; std::complex<double> unitComplex(0,1);
	if (x > 500) {expTimeBessel = 1/sqrt(2*M_PI*x) * (1 - (4*mR*mR-1)/(8*x) * (1 - (4*mR*mR-9)/(16*x)));} // Use asymptotic form for large chi
	else {expTimeBessel = exp(-x) * boost::math::cyl_bessel_i(mR, x);}

	return surfaceDensity() * ((2 * M_PI * m_kR)/(epicycleFreq() * x)) * 
	(-mR * expTimeBessel * exp(-unitComplex * (mR * epicycleFreq() + m_fourierHarmonic * circularFreq()))*timeDelta); // We need a way to caluclate the bessel function including exp(-chi) so it doesn't diverge
} 

/* Calculating the Analytic BF */
/* --------------------------- */ 

double SpiralBasis::scriptWElementAnalytic(const int m1, const double jA, const double rG) const {
	std::complex<double> unitComplex(0,1);
	double value  = (exp(unitComplex * (m_kR * rG - m1 * 0.5 * M_PI)) + exp(-unitComplex * (m_kR * rG - m1*0.5 * M_PI))).real();
 	return boost::math::cyl_bessel_j(m_fourierHarmonic, jA * m_kR) * (potential(rG)/(cos(m_kR * rG))) * value; 
}


#endif 