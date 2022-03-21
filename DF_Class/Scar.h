#ifndef SCAR 
#define SCAR

#include <cmath>

class Scar
{
public:
	Scar(double Lcirc, double alpha, double w = 0, double a = 0, double b = 0, double c = 0, double d = 0) : 
	m_Lcirc{Lcirc}, m_pattern{(1-sqrt(2)/2)/m_Lcirc}, m_jacobi{0.5 + log(Lcirc) + m_pattern*m_Lcirc}, // We have assumed a Mestel background with vc = 1
	m_alpha{alpha}, m_a{a}, m_b{b}, m_c{c}, m_d{d}, m_w{w}
	{}
	~Scar() {}

	double groove(const double E, const double L, const double equivalentL) const;
	double patternSpeed() const {return m_pattern;}
	double jacobi() const {return m_jacobi;}

	void outParams() const {std::cout << m_Lcirc << " " << m_pattern << " " << m_jacobi << '\n';}
//private:
	const double m_Lcirc, m_pattern, m_jacobi, m_alpha, m_a, m_b, m_c, m_d, m_w; 

	double getX(const double equivalentL) const;

	bool isInWidth(const double E, const double L, double x) const;
	double widthFrac(const double E, const double L, double x) const; 
	double jacobi(const double E, const double L) const {return E + L * m_pattern;}
};


double Scar::getX(const double equivalentL) const { 
	return 1- equivalentL/m_Lcirc; 
}

double Scar::groove(const double E, const double L, const double equivalentL) const { 
	double x{getX(equivalentL)}; 
	if (x < 0) {x = 0;} // We can get rounding errors

	return 1 + widthFrac(E,L,x) * (m_a * pow(1-x, -m_alpha) + m_b + m_c * x + m_d * x * x - 1);	
}

bool Scar::isInWidth(const double E, const double L, double x) const {
	if (x > .65) {return false;}
	if (abs(jacobi(E,L) - m_jacobi) < m_w)  {return true;} 
	else {return false;}
}

double Scar::widthFrac(const double E, const double L, double x) const { 
	if (!isInWidth(E,L, x)) {return 0;}
	return cos(0.5 * M_PI * ((jacobi(E,L) - m_jacobi)/(m_w))); 
}

#endif