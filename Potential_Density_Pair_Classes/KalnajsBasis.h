#ifndef KALNAJSBASIS
#define KALNAJSBASIS

#include <vector>

#include "PotentialDensityPair.h"

class KalnajsBasis : public PotentialDensityPair
{

public:
	KalnajsBasis(int m_radialIndex, int m_fourierHarmonic)
		: PotentialDensityPair(m_radialIndex, m_fourierHarmonic),
		m_kalnajsIndex{4}, m_kalnajsScale{20},
		m_alphaPrime{alphaPrime()},
		m_betaPrime{betaPrime()} 
		{} 
	
	virtual ~KalnajsBasis() {}

	double density(double r) const; 
	double potential(double r) const; 

private:
	const int m_kalnajsIndex;
	const double m_kalnajsScale;

	// Into these we will absorb script P and script S
	const std::vector<double> m_alphaPrime;
	const std::vector<double>  m_betaPrime;


	double scriptP() const;
	double scriptS() const;

	double alphaElement(int i, int j) const;
	double  betaElement(int j) const;

	std::vector<double> alphaPrime() const;
	std::vector<double>  betaPrime() const;


};

#endif