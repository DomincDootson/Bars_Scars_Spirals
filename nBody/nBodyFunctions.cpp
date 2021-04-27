#include <iostream>
#include <Eigen/Dense>

#include "NBody_Classes/NBodyPerturbation.h"
#include "NBody_Classes/NBodyBar.h"
#include "../DF_Class/Mestel.h"
#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "../Bar2D/Bar2D.h"

#include <cmath>

const int NUMBPARTICLES{100000};
const int NSTEPS{20000};
const double TIMESTEP{0.001};


Eigen::VectorXcd gaussianBar(const PotentialDensityPairContainer<GaussianLogBasis> & pd) {
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(pd.maxRadialIndex() + 1);
	for (int i = 0; i < pd.maxRadialIndex() + 1; ++i){
		coeff(i) = 0.01*pd.potential(2.06271, i);
	}
	return - (pd.scriptE()).inverse() * coeff;
}

std::string evolutionFilename(int runNumber) {
	return "../Plotting/GaussianTorque/Evolution_" + std::to_string(runNumber) + ".csv";
} 

std::string coefficentFilename(int runNumber) {
	return "../Plotting/GaussianTorque/Coefficent_" + std::to_string(runNumber) + ".csv";
} 


void barEvolution() 
{
	std::vector<double> params{24, .5, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> pd(params, 24, 2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(24+1);
	std::complex<double> unitComplex(0,1);
	

	for (int i =0; i < 5; ++i){
		std::cout << '\n' << '\n' << "Realisation: " << i << '\n';

		coeff = gaussianBar(pd);	
		Bar2D bar(coeff, 0.1, "../Bar2D/barSize.out");

		NBodyBar nbodyBar(NUMBPARTICLES, NSTEPS, TIMESTEP, pd, bar);  
		nbodyBar.testParticleEvolution(coefficentFilename(i), evolutionFilename(i), 0);
	}
}
