#ifndef ANGULARMOMENTUMSCAR
#define ANGULARMOMENTUMSCAR

class AngularMomentumScar
{
public:
	AngularMomentumScar(double position, double width, double depth) : 
	m_position{position}, m_width{width}, m_depth{depth}
	{if (m_depth < -1) {std::cout << "Depth of AM scar must be > -1\n"; exit(0);}}
	
	~AngularMomentumScar() {}

	double groove(const double E, const double L, const double equivalentL = 1) const {
		if (abs(L-m_position) < m_width) {return grooveFunction(L);}
		else {return 1;}
	} 

private:
	const double m_position, m_width, m_depth;

	double grooveFunction(const double L) const {return 1 + m_depth * pow(cos((M_PI/2) * (L-m_position)/m_width),10);}
};

#endif