#ifndef GAUSSIANLOGBASIS
#define GAUSSIANLOGBASIS

#include "PotentialDensityPair.h"

class GaussianLogBasis : public PotentialDensityPair
{
public:
	GaussianLogBasis(int m_radialIndex, int m_fourierHarmonic)
		: PotentialDensityPair(m_radialIndex, m_fourierHarmonic), m_maxIndex{24}, m_r0{r0()}, m_sigma{sigma()}, 
		m_innerRadius(.5), m_outerRadius(15)
		{} 
	
	virtual ~GaussianLogBasis() {}

	double density(double r) const; 
	double potential(double r) const; 

private:
	const int m_maxIndex;
	const double m_r0, m_sigma, m_innerRadius, m_outerRadius;
	
	double r0() const;
	double sigma() const;
};

#endif