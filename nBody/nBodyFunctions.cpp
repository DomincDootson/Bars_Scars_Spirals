#include <iostream>
#include <Eigen/Dense>

#include "NBody_Classes/NBodyPerturbation.h"
#include "NBody_Classes/NBodyBar.h"
#include "NBody_Classes/NBodySpiral.h"


#include "../DF_Class/Mestel.h"
#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"

#include "../Bar2D/Bar2D.h"
#include "../Spiral2D/Spiral2D.h"

#include <cmath>


const int NUMBPARTICLES{500000};
const int NSTEPS{30000};   //300000
const double TIMESTEP{0.01};


Eigen::VectorXcd gaussianBar(const PotentialDensityPairContainer<GaussianLogBasis> & pd) {
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(pd.maxRadialIndex() + 1);
	for (int i = 0; i < pd.maxRadialIndex() + 1; ++i){
		coeff(i) = 0.001*pd.potential(2.06271, i);
	}
	return - (pd.scriptE()).inverse() * coeff;
}

std::string evolutionFilenameGaussian(int runNumber) {
	return "../Plotting/GaussianTorque/Evolution_" + std::to_string(runNumber) + ".csv";
} 

std::string coefficentFilenameGaussian(int runNumber) {
	return "../Plotting/GaussianTorque/Coefficent_" + std::to_string(runNumber) + ".csv";
} 


std::string evolutionFilenameKalnajs(int runNumber) {
	return "../Plotting/KalnajsTorque/Evolution_" + std::to_string(runNumber) + ".csv";
} 

std::string coefficentFilenameKalnajs(int runNumber) {
	return "../Plotting/KalnajsTorque/Coefficent_" + std::to_string(runNumber) + ".csv";
} 

std::string stemFilename(const std::string & stem, const int runNumber) {
	return "../Plotting/KalnajsTorque/" + stem + "_" + std::to_string(runNumber) + ".csv";
}


void barEvolutionGaussian() 
{
	std::vector<double> params{50, .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> pd(params, 50, 2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(50+1);
	std::complex<double> unitComplex(0,1);
	

	for (int i =0; i < 5; ++i){
		std::cout << '\n' << '\n' << "Realisation: " << i << '\n';

		coeff = gaussianBar(pd);	
		Bar2D bar(coeff, 0.01, "../Bar2D/barSize.out");

		NBodyBar nbodyBar(1, NSTEPS, TIMESTEP, pd, bar);  
		nbodyBar.testParticleEvolution(coefficentFilenameGaussian(i), evolutionFilenameGaussian(i), 0);
		exit(0); 
		//nbodyBar.nBodyEvolution(coefficentFilenameGaussian(i), evolutionFilenameGaussian(i), 0);
	}
}

void barEvolutionKalnajs(const std::string & stem, const bool isSelfConsistent, const double littleSigma)
{
	std::vector<double> params{4, 20};
 	PotentialDensityPairContainer<KalnajsNBasis> pd("../Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_20_2.dat");

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(48+1);
	coeff(0) = 0.01;
	

	for (int i =0; i < 5; ++i){
		std::cout << '\n' << '\n' << "Realisation: " << i << '\n';

		Bar2D bar(coeff, 0.05, "../Bar2D/barSize.out");

		NBodyBar nbodyBar(NUMBPARTICLES, NSTEPS, TIMESTEP, pd, bar, littleSigma);  
		if (isSelfConsistent) {nbodyBar.nBodyEvolution(coefficentFilenameKalnajs(i), stemFilename(stem, i), 0);}
		else {nbodyBar.testParticleEvolution(coefficentFilenameKalnajs(i), stemFilename(stem, i), 0);} 

	}
}

void orbitSection() 
{
	std::vector<double> params{24, .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> pd(params, 24, 2);	

	std::vector<double> omega{0.5, .1, .5}, strength{1, 1, 0};
	std::vector<std::string> filenames{"../Plotting/Orbit_Sections/ResonantSections.csv", "../Plotting/Orbit_Sections/NonResonantSections.csv", "../Plotting/Orbit_Sections/NoPerturbationSections.csv"};

	for (int i = 0; i < 1; ++i) {
		Eigen::VectorXcd coeff = strength[i]*gaussianBar(pd);	
		Bar2D bar(coeff, omega[i]);
		NBodyBar nbodyBar(1, 200, 0.0005, pd, bar);  
		//nbodyBar.barOrbitSections(filenames[i], false);


		nbodyBar.angularMomentumSections(filenames[i], false);
		//nbodyBar.countTrappedOrbits(); 
		exit(0);
		std::cout << "Finished the funciton: " << i << '\n';
	}
}


void kalanajTest()
{
	std::vector<double> params{4, 20};
 	PotentialDensityPairContainer<KalnajsBasis> pd(params, 10, 2);


 	 NBodyPerturbation nbody(NUMBPARTICLES, NSTEPS, TIMESTEP, pd);
 	 nbody.nBodyEvolution("../Plotting/Evolution.csv");

}


void spiralTesting() {
 	PotentialDensityPairContainer<KalnajsNBasis> pd("../Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

 	Spiral2D spiral(pd, -1, 5, 0.001);

 	NBodySpiral nbody(500000, 100000*0.5, TIMESTEP, pd, spiral);
 	nbody.nBodyEvolution("../Plotting/Spiral_Data/nBodyEvolution.csv");
}
