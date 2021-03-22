#include <iostream>
#include <Eigen/Dense>

#include "nBody_Classes/NBodyClass.h"
#include "../DF_Class/Mestel.h"
#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"

#include <cmath>

int main(){
	/*Mestel DF;
	DF.cumulativeDensity("cumulativeDensity.csv");*/


	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd(params, 10,1);

	NBodyClass nbody(100000, 20000, 1, 0.0025, pd); 


	return 0;
} 
