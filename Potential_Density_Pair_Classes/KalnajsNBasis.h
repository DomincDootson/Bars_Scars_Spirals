#ifndef KALNAJSNBASIS
#define KALNAJSNBASIS 

#include <fstream>
#include <vector>

#include "PotentialDensityPair.h"

class KalnajsNBasis : public PotentialDensityPair
{
public:
	KalnajsNBasis(std::ifstream & inFile, const double step, const double rMax, const int n, const int l) : PotentialDensityPair(n, l),
	m_step{step}, m_rMax{rMax},
	v_potential(((int) 1/m_step)+1), v_density(((int) 1/m_step)+1) // It includes the endpoint. 
	{
		for (int i = 0; i < v_potential.size(); i++) {inFile >> v_potential[i];}
		for (int i = 0; i < v_density.size(); i++) {inFile >> v_density[i];}
	}
	
	~KalnajsNBasis() {}

	double density(double r) const;
	double potential(double r) const; 
	

private:
 const double m_step, m_rMax;
 std::vector<double> v_potential, v_density; 	

 double ind(double r) const {return r/(m_step*m_rMax);}
 int floorIndex(double r) const {return (int) floor(r/(m_step*m_rMax));}
 int ceilIndex(double r) const {return (int) ceil(r/(m_step*m_rMax));} 
};

double KalnajsNBasis::density(double r) const {
	double index{ind(r)}; if (r>m_rMax) {return 0;} 
	int lowerIndex{floorIndex(r)}, upperIndex{ceilIndex(r)};
	return v_density[lowerIndex] + ((index-lowerIndex)) * (v_density[upperIndex]-v_density[lowerIndex]);
}

double KalnajsNBasis::potential(double r) const {
	if (r>m_rMax) {return 0;}
	double index{ind(r)};
	int lowerIndex{floorIndex(r)}, upperIndex{ceilIndex(r)};
	
	return v_potential[lowerIndex] + ((index-lowerIndex)) * (v_potential[upperIndex]-v_potential[lowerIndex]);
}
#endif

// Things to check

// Compling them
// Returning the correct potential and density 
