#include <iostream>
#include <Eigen/Dense>

#include "NBody_Classes/NBodyPerturbation.h"
#include "NBody_Classes/NBodyBar.h"
#include "../DF_Class/Mestel.h"
#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"

#include <cmath>

void perturbationFile(std::string filename)
{
	std::ofstream out(filename); out << 10 << '\n';
	for (int time =0; time < 2000; ++time){
		for (int i =0; i < 10; ++i){
			if (i==0 || i==1 || i==2) {out << 0.01*sin(M_PI*(time/2000.0)) << " " << 0 << " ";} // 0.01*sin(M_PI*(time/2000.0))
			else {out << 0 << " " << 0 <<" ";}
		}
		out << 0 << " " << 0 <<'\n';
	}
}

int main(){
	/*Mestel DF;
	DF.cumulativeDensity("cumulativeDensity.csv");*/
	std::vector<double> params{4, 20};

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(11);
	coeff(0) = 0.01;

	PotentialDensityPairContainer<KalnajsBasis> pd(params, 10,2);
	NBodyBar nbody(100000, 20000, 0.001, pd, coeff, 0.5);  //100000

	nbody.testParticleEvolution("coeffEvolution.csv", "barEvolution.csv", 0);







	//nbody.testParticleEvolution("diskCoef.csv", "barEvolution.csv", 0);


	return 0;
} 


// Things to do
//	1) 