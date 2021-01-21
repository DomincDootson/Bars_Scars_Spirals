#ifndef MESTEL
#define MESTEL

#include <cmath>

#include "DFClass.h"

class Mestel : public DFClass
{
public:
	Mestel() : m_vc{1}, m_r0{1}, m_littleSigma{.25}, m_q{pow(m_vc/m_littleSigma,2) - 1}, m_xi{1}, m_rInner{1}, m_rOuter{10}, m_nuTaper{4}, m_muTaper{5}  {}
	~Mestel() {}

	virtual double potential(double radius) const {return m_vc*m_vc*log(radius/m_r0);}
	virtual double distFunc(double E, double J) const;

	
private:
	const double m_vc, m_r0, m_littleSigma, m_q, m_xi, m_rInner, m_rOuter, m_nuTaper, m_muTaper;
	double innerTaper(double J) const;
	double outerTaper(double J) const;
};

double Mestel::innerTaper(double J) const{
	return pow(J, m_nuTaper) / (pow(m_rInner*m_vc, m_nuTaper) + pow(J, m_nuTaper));
}

double Mestel::outerTaper(double J) const{
	return 1/(1+pow(J/(m_rOuter*m_vc), m_muTaper));
}


double Mestel::distFunc(double E, double J) const // Note that we have set G = 1
{
	return  innerTaper(J)*outerTaper(J) *(m_xi)*pow((J/(m_r0*m_vc)), m_q) * exp(-E/pow(m_littleSigma,2)) * 
	(((pow(m_vc,2)/(2*M_PI*m_r0)) * pow(m_vc,m_q)) * pow(pow(2, 0.5*m_q) * 
	sqrt(M_PI)*tgamma(0.5*m_q+0.5)*pow(m_littleSigma,2+m_q), -1)); 

}

#endif