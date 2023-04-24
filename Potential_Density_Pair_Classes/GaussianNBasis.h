#ifndef GAUSSIANNBASIS
#define GAUSSIANNBASIS 

#include <fstream>
#include <vector>

#include "PotentialDensityPair.h"

class GaussianNBasis : public PotentialDensityPair
{
public:
	GaussianNBasis(std::ifstream & inFile, const double step, const double rMax, const int n, const int l) : PotentialDensityPair(n, l),
	m_step{step}, m_rMax{rMax}, m_r0{r0(n)}, m_sigma{sigma(n)},
	v_potential(((int) 1/m_step)+1), v_density(((int) 1/m_step)+1) // It includes the endpoint. 
	{
		for (int i = 0; i < v_potential.size(); i++) {inFile >> v_potential[i];}
		for (int i = 0; i < v_density.size(); i++) {inFile >> v_density[i];}
	}
	
	~GaussianNBasis() {}

	double density(double r) const;
	double potential(double r) const; 

	double rMax() const {return m_rMax;} 

	double maxRadius() const {return m_rMax;}
	

private:
 const double m_step, m_rMax, m_r0, m_sigma;
 std::vector<double> v_potential, v_density; 	

 double ind(double r) const {return r/(m_step*m_rMax);}
 int floorIndex(double r) const {return (int) floor(r/(m_step*m_rMax));}
 int ceilIndex(double r) const {return (int) ceil(r/(m_step*m_rMax));} 

 double r0(const int index, const int maxIndex = 48, const double innerRadius=0.5, const double outerRadius=15) const;
 double sigma(const int index, const int maxIndex = 48, const double innerRadius=0.5, const double outerRadius=15) const;
};

double GaussianNBasis::r0(const int index, const int maxIndex, const double innerRadius, const double outerRadius) const
{
	double t = pow(outerRadius/innerRadius, 1/((double) maxIndex));
	return innerRadius * pow(t, index);
}

double GaussianNBasis::sigma(const int index, const int maxIndex, const double innerRadius, const double outerRadius) const
{
	double t = pow(outerRadius/innerRadius, 1/((double) maxIndex));
	return innerRadius * pow(t, index+1) - m_r0;
}


double GaussianNBasis::potential(double r) const {
	if (r>m_rMax) {return 0;}
	double index{ind(r)};
	int lowerIndex{floorIndex(r)}, upperIndex{ceilIndex(r)};
	
	return v_potential[lowerIndex] + ((index-lowerIndex)) * (v_potential[upperIndex]-v_potential[lowerIndex]);
}

double GaussianNBasis::density(double r) const {return 1/(m_sigma*sqrt(2*M_PI)) * exp(-0.5 *pow(r-m_r0,2)/(m_sigma*m_sigma));}
#endif

// Things to check

// Compling them
// Returning the correct potential and density 
