#ifndef PERTURBEDFCONTAINER
#define PERTURBEDFCONTAINER 

#include <vector>
#include <string>
#include <iostream>

#include "PerturbedDF.h"

class PerturbedDFContainer
{
public:
	PerturbedDFContainer(); // Assume that the first two are zero (ICs)
	~PerturbedDFContainer();

	void write2file(const std::string & filename) const; // write function that save a DF to file 


private:
	const int m_m1; 


	std::vector<PerturbedDF> v_perturbedDF; 
	
};

#endif