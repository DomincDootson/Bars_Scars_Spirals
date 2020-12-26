#ifndef POTENTIALDENSITYPAIR
#define POTENTIALDENSITYPAIR

#include <vector>
#include <cmath>

class PotentialDensityPair
{
public:
	PotentialDensityPair(int n, int l) 
	:
		m_fourierHarmonic{l}, m_radialIndex{n}, m_G(1){}
	
	virtual ~PotentialDensityPair() {}

	virtual double density(double r) const = 0; 
	virtual double potential(double r) const = 0; 


	double scriptWElement(int const m1, std::vector<double> const &radii, std::vector<double> const &theta1, std::vector<double> const &theta2, std::vector<double> const &theta1Deriv);

protected:
	const int m_fourierHarmonic, m_radialIndex;
	const int m_G;
	
};


double PotentialDensityPair::scriptWElement(int const m1, std::vector<double> const &radii, std::vector<double> const &theta1, std::vector<double> const &theta2, std::vector<double> const &theta1Deriv)
{
	int nstep{static_cast<int> (radii.size())};
	double integral{0}, upperU{0.5 * M_PI}, stepSize{upperU/nstep};
	for (int i = 1; i< nstep-1; ++i){
		integral += (theta1Deriv[i]) * potential(radii[i]) * sin(2*stepSize*i) *stepSize * cos(m1 * theta1[i] + m_fourierHarmonic*theta2[i]);
	}
	return (radii.back() - radii.front()) * (integral/M_PI);
}
#endif