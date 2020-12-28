#ifndef GAUSSIANLOGBASIS
#define GAUSSIANLOGBASIS

#include "PotentialDensityPair.h"

#include <iostream>

class GaussianLogBasis : public PotentialDensityPair
{
public:
	GaussianLogBasis(int m_radialIndex, int m_fourierHarmonic)
		: PotentialDensityPair(m_radialIndex, m_fourierHarmonic), m_maxIndex{24}, m_innerRadius{.5}, m_outerRadius{15},
		m_r0{r0()}, m_sigma{sigma()}
		
		{
			std::cout << m_r0 << " " << m_sigma << '\n';
		} 
	
	virtual ~GaussianLogBasis() {}

	double density(double r) const; 
	double potential(double r) const; 

private:
	const int m_maxIndex;
	const double m_innerRadius, m_outerRadius, m_r0, m_sigma;
	
	double r0() const;
	double sigma() const;
};

double GaussianLogBasis::r0() const
{
	double t = pow(m_outerRadius/m_innerRadius, 1/((double) m_maxIndex));
	return m_innerRadius * pow(t, m_radialIndex);
}

double GaussianLogBasis::sigma() const
{
	double t = pow(m_outerRadius/m_innerRadius, 1/((double) m_maxIndex));
	return m_innerRadius * pow(t, m_radialIndex+1) - m_r0;
}

double GaussianLogBasis::potential(double radius) const
{
	return m_sigma*m_sigma*exp(-pow(((radius-m_r0)/m_sigma),2)); 
}

double GaussianLogBasis::density(double radius) const
{
	return potential(radius) *
	(4*pow(radius,4) - 4*pow(radius*m_sigma,2) - pow(m_fourierHarmonic*m_sigma*m_sigma,2) +2*radius*m_r0*(-4*radius*radius +m_sigma*m_sigma+2*radius*m_r0))
	/(4*3.14159 * pow(radius*m_sigma*m_sigma,2));
}


#endif