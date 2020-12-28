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

double pochhammer(double a, int i) // Pochhammer function
{
	if (i==0)
	{
		return 1;
	}
	
	double output{a};

	for (int j = 1; j < i; ++j)
	{
		output *= (j+a);
	}
	return output;
}


double KalnajsBasis::scriptP() const  
	{
		int k{m_kalnajsIndex}, l{m_fourierHarmonic}, n{m_radialIndex};
		return pow(((2*k + l + 2*n + 0.5) * tgamma(2*k + l + n + 0.5) * tgamma(l + n + 0.5)) / 
			   (tgamma(2*k + n + 1) * pow(tgamma(l+1),2) * tgamma(n+1)),0.5);
	}


double KalnajsBasis::scriptS() const 
	{
		int k{m_kalnajsIndex}, l{m_fourierHarmonic}, n{m_radialIndex};
		return (tgamma(k+1) / (M_PI * tgamma(2*k+1) * tgamma(k+0.5)))*pow(((2*k+l+2*n+0.5) * tgamma(2*k+n+1) * tgamma(2*k+l+n+0.5))/
			(tgamma(l+n+0.5) * tgamma(n+1)) , 0.5);
	} 


double KalnajsBasis::alphaElement(int i, int j) const 
{
	int k{m_kalnajsIndex}, l{m_fourierHarmonic}, n{m_radialIndex};
	return (pochhammer(-k, i) * pochhammer(l+0.5, i) * pochhammer(2*k+l+n+0.5,j) * pochhammer(i+l+0.5,j) * pochhammer(-n,j)) / 
	(pochhammer(l+1,i) * pochhammer(1,i) * pochhammer(l+i+1,j) * pochhammer(l+0.5,j) * pochhammer(1,j));
}

double KalnajsBasis::betaElement(int j) const 
{
	int k{m_kalnajsIndex}, l{m_fourierHarmonic}, n{m_radialIndex};
	return (pochhammer(2*k+l+n+0.5,j) * pochhammer(k+1,j) * pochhammer(-n,j)) / 
	(pochhammer(2*k+1,j) * pochhammer(k+0.5, j) * pochhammer(1,j));
}


std::vector<double> KalnajsBasis::alphaPrime() const 
{
	std::vector<double> alphaPrime;
	double scriptP = KalnajsBasis::scriptP();

	for (int j = 0; j <= m_radialIndex; ++j)
	{
		for (int i = 0; i <= m_kalnajsIndex; ++i)
		{
			alphaPrime.push_back(alphaElement(i,j) * scriptP);
		}
	}

	return alphaPrime;
}
	
std::vector<double>  KalnajsBasis::betaPrime() const 
{
	std::vector<double> betaPrime;
	double scriptS = KalnajsBasis::scriptS();
	for (int j = 0; j<=m_radialIndex; ++j)
	{
		betaPrime.push_back(betaElement(j) * scriptS);
	}

	return betaPrime;
}


double KalnajsBasis::potential(double r) const 
{
	double density{0}, rScaled{r/m_kalnajsScale};
	if (rScaled>1){return 0;}

	for (int j = 0; j <= m_radialIndex; ++j)
	{
		for (int i = 0; i <= m_kalnajsIndex; ++i)
		{
			density += m_alphaPrime[i + j*(m_kalnajsIndex+1)] * pow(rScaled, 2*i+2*j);
		}
	}

	return -pow(m_G/m_kalnajsScale, 0.5)*pow(rScaled, m_fourierHarmonic) * density;
}


double KalnajsBasis::density(double r) const 
{
	double potential{0}, rScaled{r/m_kalnajsScale};
	if (rScaled>1){return 0;}
	for (int j = 0; j <= m_radialIndex; ++j)
	{
		potential += m_betaPrime[j] * pow(1-rScaled*rScaled, j);
	} 

	return (pow(-1,m_radialIndex)*pow(m_G,-0.5)*pow(m_kalnajsScale,-1.5)) 
	* pow(1-rScaled*rScaled, m_kalnajsIndex-0.5) *pow(rScaled, m_fourierHarmonic) * potential;
} 

#endif