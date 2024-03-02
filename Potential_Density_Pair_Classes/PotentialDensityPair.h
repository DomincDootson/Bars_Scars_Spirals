#ifndef POTENTIALDENSITYPAIR
#define POTENTIALDENSITYPAIR

#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Dense>

class PotentialDensityPair
{
public:
	PotentialDensityPair(int n, int l) : m_fourierHarmonic{l}, m_radialIndex{n}, m_G(1), m_rMax{0} 
	{}
	
	virtual ~PotentialDensityPair() {}

	virtual double density(double r) const = 0; 
	virtual double potential(double r) const = 0; 


	double scriptWElement(int const m1, std::vector<double> const &radii, std::vector<double> const &theta1, std::vector<double> const &theta2, std::vector<double> const &theta1Deriv);


	Eigen::ArrayXXcd   densityGrid(const int nGrid, const double rMax) const;
	Eigen::ArrayXXcd potentialGrid(const int nGrid, const double rMax) const;

protected:

	const int m_fourierHarmonic, m_radialIndex;
	const int m_G;
	
	double m_rMax;
};


double PotentialDensityPair::scriptWElement(int const m1, std::vector<double> const &radii, std::vector<double> const &theta1, std::vector<double> const &theta2, std::vector<double> const &theta1Deriv)
{
	int nstep{static_cast<int> (radii.size())};
	double integral{0}, upperU{0.5 * M_PI}, stepSize{upperU/(nstep-1)};
	for (int i = 1; i< nstep-1; ++i){
		integral += (theta1Deriv[i]) * potential(radii[i]) * sin(2*stepSize*i) * stepSize * cos(m1*theta1[i] + m_fourierHarmonic*theta2[i]);
	}
	
	if (integral != integral){std::cout << radii.front() << " " << radii.back() << " " << m_radialIndex << " " << m1 <<'\n'; std::cout << integral <<'\n';exit(0) ;return 0;}
	return (radii.back() - radii.front()) * (integral/M_PI);
}



double angle(double x, double y) { return atan2(y, x); }

Eigen::ArrayXXcd PotentialDensityPair::densityGrid(const int nGrid, const double rMax) const
{	
	Eigen::ArrayXXcd grid = Eigen::ArrayXXcd(nGrid, nGrid);
	double x{}, y{}, theta{}, r{}, spacing{2*rMax/((double) nGrid-1)}, centre{0.5*(nGrid-1)};
	std::complex<double> unitComplex(0,1);
	for (int i = 0; i < grid.rows(); ++i)
	{
		for (int j = 0; j < grid.cols(); ++j)
		{
			x = spacing * (i - centre); y = spacing * (j - centre);
			r = sqrt(x*x+y*y); theta = angle(x,y);
			
			grid(i,j) = exp(-unitComplex * (theta * m_fourierHarmonic)) * density(r); // I am 99.8 % sure that the minus sign shouldn't be hear, but as it will only affect plots, I'll leave it
		}
	}
	return grid;
}

Eigen::ArrayXXcd PotentialDensityPair::potentialGrid(const int nGrid, const double rMax) const
{
	Eigen::ArrayXXcd grid = Eigen::ArrayXXcd::Zero(nGrid, nGrid);
	double x{}, y{}, theta{}, r{}, spacing{2*rMax/((double) nGrid-1)}, centre{0.5*(nGrid-1)};
	std::complex<double> unitComplex(0,1);
	for (int i = 0; i < grid.rows(); ++i)
	{
		for (int j = 0; j < grid.cols(); ++j)
		{
			x = spacing * (i - centre); y = spacing * (j - centre);
			r = sqrt(x*x+y*y); theta = angle(x,y);
			
			grid(i,j) = exp(unitComplex * (theta * m_fourierHarmonic)) * potential(r); // we have remove a minus sign in the expo, NO SIGN SHOULD BE CORRECT
		}
	}
	return grid;
}

#endif