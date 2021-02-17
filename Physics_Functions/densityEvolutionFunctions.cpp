#include <iostream>
#include <typeinfo>
#include <Eigen/Dense>
#include <vector>
#include <complex>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/Mestel.h"

#include "../Volterra_Solver/VolterraSolver.h"
#include "../Volterra_Solver/ExpansionCoeff.h"

#include "../Bar2D/Bar2D.h"


#include "densityEvolutionFunctions.h"

void kalnajsKernelsVaryingSigma(int l)
{
	std::vector<double> littleSigma{.25, .35, .50};
	
	std::string file = "Kalnajs/Kalnajs_4_20";
	ActionAngleBasisContainer test(file, 10, l, 5, 101, 20);
	
	VolterraSolver solver(10, l, 2000, 0.025);
	for (int i = 0; i < littleSigma.size(); ++i){
		
		std::cout << "Calculating kernel for littleSigma: " << littleSigma[i] << '\n';
		Mestel DF(1, 1, littleSigma[i]);
		
		std::string kernel = "Disk_Kicking/littleSigma_" + std::to_string((int) round(littleSigma[i]*100)) 
		+ "/Kalnajs"+ "_" + std::to_string(l) +".out";
		
		solver.generateKernel(kernel, DF, test);
	}
}


std::string perturbationFilename(double littleSigma, double radius, int fourierHarmonic)
{
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Perturbation_" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + ".out";
}

std::string kickingDensityFilename(double littleSigma, double radius, int fourierHarmonic)
{
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Density" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + ".out";
}

template <typename T>
void outPutPerturbation(const T & pd, double littleSigma, double radius, int numbTimeSteps = 2000)
{
	Eigen::VectorXcd t0Coeff(pd.maxRadialIndex()+1);
	for (int i = 0; i <= pd.maxRadialIndex(); ++i){
		t0Coeff[i] = pd.potential(radius, i);
	}
	// TODO - Include scipt E multplication
	ExpansionCoeff holding(t0Coeff, numbTimeSteps, pd.maxRadialIndex());
	holding.writePerturbation2File(perturbationFilename(littleSigma, radius, pd.fourierHarmonic()));
}


void diskKickingPerturbations()
{
	std::vector<double> radii {.1,.5,1, 2, 3, 5};

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0(params, 10,0), pd1(params, 10, 1), pd2(params, 10, 2);

	for (auto i = radii.begin(); i != radii.end(); ++i){
		outPutPerturbation(pd0, .25, *i); outPutPerturbation(pd1, .25, *i); outPutPerturbation(pd2, .25, *i);
		outPutPerturbation(pd0, .35, *i); outPutPerturbation(pd1, .35, *i); outPutPerturbation(pd2, .35, *i);
		outPutPerturbation(pd0, .50, *i); outPutPerturbation(pd1, .50, *i); outPutPerturbation(pd2, .50, *i);
	}
}

template <typename T>
void individualDiskKicking(const T & pd, double littleSigma, double radius)
{
	int angHarmonic{pd.fourierHarmonic()};
	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	std::string densityFile{kickingDensityFilename(littleSigma, radius, angHarmonic)};
	std::string kernel = "Disk_Kicking/littleSigma_" + std::to_string((int) round(littleSigma*100)) + "/Kalnajs"+ "_" + std::to_string(angHarmonic) +".out";

	VolterraSolver solver(kernel, pd.maxRadialIndex(), angHarmonic, 2000, 0.025);
	solver.activeFraction(.25);
	solver.volterraSolver(pd, densityFile, perturbationFile, true);
}




void diskKicking()
{
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0(params, 10,0), pd1(params, 10, 1), pd2(params, 10, 2);

	std::vector<double> littleSigma{.25, .35, .50}, radii{.1, .5, 1, 2, 5};

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto radius = radii.begin(); radius != radii.end(); ++radius){
			individualDiskKicking(pd0, *s, *radius);
			individualDiskKicking(pd1, *s, *radius);
			individualDiskKicking(pd2, *s, *radius);
		}
	}
}



