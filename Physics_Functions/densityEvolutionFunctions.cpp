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

const int NUMBTIMESTEPS{2000}; // some global paramters for the Volterra Solver
const double DELTAT{0.025};
const double XI{.5};

std::string perturbationFilename(double littleSigma, double radius, int fourierHarmonic)
{
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Perturbation_" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + ".out";
}

std::string kickingDensityFilename(double littleSigma, double radius, int fourierHarmonic)
{
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Density" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + ".csv";
}

std::string kernelFilename(double littleSigma, int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Kalnajs" + "_" +std::to_string(fourierHarmonic) + ".out";
}



std::vector<double> littleSigmas(){
	return {0.25, 0.35, 0.50};
}

std::vector<double> perturbationRadii() {
	return {.1,.5,1, 2, 3, 5,15};
}

void kalnajsKernelsVaryingSigma(int l)
{
	std::vector<double> littleSigma{littleSigmas()};
	
	std::string file = "Kalnajs/Kalnajs_4_20";
	ActionAngleBasisContainer test(file, 10, l, 5, 101, 20);
	
	VolterraSolver solver(10, l, NUMBTIMESTEPS, DELTAT);
	for (int i = 0; i < littleSigma.size(); ++i){
		
		std::cout << "Calculating kernel for littleSigma: " << littleSigma[i] << '\n';
		Mestel DF(1, 1, littleSigma[i]);
		
		std::string kernel{kernelFilename(littleSigma[i], l)};
		
		solver.generateKernel(kernel, DF, test);
	}
}

template <typename T>
void outPutPerturbation(const T & pd, double littleSigma, double radius, int numbTimeSteps = NUMBTIMESTEPS)
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
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0(params, 10,0), pd1(params, 10, 1), pd2(params, 10, 2);

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto r = radii.begin(); r != radii.end(); ++r){
			outPutPerturbation(pd0, *s, *r);	
			outPutPerturbation(pd1, *s, *r);	
			outPutPerturbation(pd2, *s, *r);	
		}
	}
}

template <typename T>
void individualDiskKicking(const T & pd, double littleSigma, double radius)
{
	int angHarmonic{pd.fourierHarmonic()};
	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	std::string densityFile{kickingDensityFilename(littleSigma, radius, angHarmonic)};
	std::string kernel{kernelFilename(littleSigma, angHarmonic)};
	std::cout << "Solving: " << littleSigma << " " << radius << " " << angHarmonic <<'\n';
	

	VolterraSolver solver(kernel, pd.maxRadialIndex(), angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(XI);
	solver.densityEvolution(pd, densityFile, perturbationFile, true);
}
 
void diskKicking()
{
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0(params, 10,0), pd1(params, 10, 1), pd2(params, 10, 2);

	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto radius = radii.begin(); radius != radii.end(); ++radius){
			individualDiskKicking(pd0, *s, *radius);
			individualDiskKicking(pd1, *s, *radius);
			individualDiskKicking(pd2, *s, *radius);
		}
	}
}


void outPutEnergy(std::ofstream & out, const std::vector<double> & energies, double littleSigma, int angHarmonic, double radius){
	int skip{10};
	out << littleSigma << ',' << angHarmonic << ',' << radius << ',';
	for (int i = 0; i<energies.size()-skip; i += skip){
		out << energies[i] <<',';
	}
	out << energies[energies.size() - skip] << '\n';
}

template <class Tbf>
std::vector<double> individualEnergyEvolution(const Tbf & pd, double littleSigma, double radius)
{
	int angHarmonic{pd.fourierHarmonic()};
	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	std::string kernel{kernelFilename(littleSigma, angHarmonic)};
	std::cout << "Solving: " << littleSigma << " " << radius << " " << angHarmonic <<'\n';

	VolterraSolver solver(kernel, pd.maxRadialIndex(), angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(XI);
	return solver.energyEvolution(pd, perturbationFile, true);
}


void energyEvolution()
{
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0{params, 10, 0}, pd1{params, 10, 1}, pd2{params, 10, 2};

	std::vector<PotentialDensityPairContainer<KalnajsBasis>> potentialDensityPairs{pd0, pd1, pd2};
	
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};

	std::ofstream out("Disk_Kicking/Energy_Evolution/KalnajsEnergyEvolution.csv");

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto pd = potentialDensityPairs.begin(); pd != potentialDensityPairs.end(); ++pd){
			for (auto r = radii.begin(); r != radii.end(); ++r){
				std::vector<double> energies = individualEnergyEvolution(*pd, *s, *r);
				outPutEnergy(out, energies, *s, pd -> fourierHarmonic(), *r);
			}
		}
	} 	
	out.close();
}


// 5) ADD SOME GLOBAL VARIABLES ?? 