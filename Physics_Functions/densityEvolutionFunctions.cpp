#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <cstring>

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
double XI{.5};

bool SELF{true}; 

void makeSelfConsistent() {SELF = true;}
void makeTestParticle()   {SELF = false;}


std::string isTestParticle()
{
	if (!SELF)
	{
		return "_test";
	}
	return "";
}

std::string perturbationFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Perturbation_" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() + ".out";
}

std::string kickingCoefFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Coeff" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() +".csv";
}

std::string kickingDensityFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Density" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() +".csv";
}

std::string density2dFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Density2D" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() + ".csv";
}

std::string kernelFilename(const double littleSigma, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Kalnajs" + "_" +std::to_string(fourierHarmonic) + ".out";
}



std::vector<double> littleSigmas(){
	return {0.25, 0.35, 0.50};
}

std::vector<double> perturbationRadii() {
	std::vector<double> holding;
	for (double radius = .5; radius <= 7; radius +=.1){
		holding.push_back(radius);
	}
	return holding;
}

void kalnajBF()
{
	Mestel DF;
	double Rka{15};
	std::vector<double> params{4, Rka};
	for (int i =0; i<=2; ++i){
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, i);

		ActionAngleBasisContainer test(10, i, 5, 101, 10);
		std::string file = "Kalnajs/Kalnajs_4_" + std::to_string((int) Rka); // Use file function name here
		test.scriptW(PD, DF, file);
	} 	
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

void maxDensityRadii(){
	std::ofstream out("../maxDensity.csv");
	std::vector littleSigma{littleSigmas()};
	for (int i = 0; i < littleSigma.size()-1; ++i) {out << littleSigma[i] << ',';}
	out << littleSigma.back() << '\n';

	for (double rInner = .1; rInner < 5; rInner += .1){
		out << rInner << ',';
		std::cout << rInner << '\n';
		for (auto i = 0; i < littleSigma.size(); ++i){
			Mestel DF(1, 1, littleSigma[i], 1, rInner);
			if (i < littleSigma.size() -1){ out << DF.maxDensityRadius() << ',';}
			else {out << DF.maxDensityRadius() << '\n';}
		}
	}
	out.close(); 
}

// Perturbation Functions // 
// ---------------------- //

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
 
void cleaningPerturbations()
{
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()}, angHarmonic{0,1,2};
	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto r = radii.begin(); r != radii.end(); ++r){
			for (auto l = angHarmonic.begin(); l != angHarmonic.end(); ++l){
				if (!remove(perturbationFilename(*s, *r, *l).c_str()))
				{
					std::cout << "Sucessfully deleted: " << perturbationFilename(*s, *r, *l) << '\n';
				}
			}
		}
	}
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

void vectorLine2File(std::ofstream & out, const std::vector<double> & vector)
{
	for (auto it = vector.begin(); it != vector.end()-1; ++it){
		out << *it << ',';
	}
	out << vector.back() << '\n';
}

std::vector<double> radii4Line(int rMax = 10, int nStep = 100){
	double step{rMax/((double) nStep-1)};
	std::vector<double> radii;

	for (int i = 0; i < nStep; ++i){
		radii.push_back(i*step);
	}
	return radii;
}

void plottingPerturbations()
{
	Eigen::VectorXcd coeff(11);
	std::vector<double> params{4, 20}, radii{perturbationRadii()}, radii4Output{radii4Line()};
	PotentialDensityPairContainer<KalnajsBasis> pd1{params, 10, 1};

	std::ofstream out("Disk_Kicking/PlottingPerturbation.csv");

	for (auto r = radii.begin(); r != radii.end(); r += 1){
		for (int n = 0; n < coeff.size(); ++n){
			coeff[n] = pd1.potential(*r, n);
		}
		out << *r << ','; 
		vectorLine2File(out, pd1.oneDpotential(radii4Output, coeff));
		std::cout << *r << '\n';
	}
	out.close();
}


// Coefficent Evolution //
// -------------------- //

void indivdualCoeffEvolution(const int nRadialHarmonics, const int angHarmonic, const double littleSigma, const double radius){
	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	std::string coeff{kickingCoefFilename(littleSigma, radius, angHarmonic)};
	std::string kernel{kernelFilename(littleSigma, angHarmonic)};
	std::cout << "Solving: " << littleSigma << " " << radius << " " << angHarmonic <<'\n';
	

	VolterraSolver solver(kernel, nRadialHarmonics , angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(XI);
	solver.coefficentEvolution(coeff, perturbationFile, SELF);
}

void coefficentEvolution()
{
	diskKickingPerturbations();
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto radius = radii.begin(); radius != radii.end(); ++radius){
			indivdualCoeffEvolution(10, 0, *s, *radius);
			indivdualCoeffEvolution(10, 1, *s, *radius);
			indivdualCoeffEvolution(10, 2, *s, *radius);
		}
	}
	cleaningPerturbations();
}


// 1D Density Evolution //
// -------------------- //


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
	solver.densityEvolution(pd, densityFile, perturbationFile, SELF);
}
 
void diskKicking()
{
	diskKickingPerturbations();
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
	cleaningPerturbations();
}

// 2D Density Evolution //
// -------------------- //

void density2D(double littleSigma, double radius, int angHarmonic){
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd{params, 10, angHarmonic};	

	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	outPutPerturbation(pd, littleSigma, radius);

	std::string kernel{kernelFilename(littleSigma, angHarmonic)};
	VolterraSolver solver(kernel, pd.maxRadialIndex(), angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(XI);

	solver.solveVolterraEquation(perturbationFile, SELF);

	std::string outFilename{density2dFilename(littleSigma, radius, angHarmonic)};
	solver.density2dEvolution(outFilename, pd);
}


// Energy Evolution //
// ---------------- // 

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
	return solver.energyEvolution(pd, perturbationFile, SELF);
}


void energyEvolution(const std::string & energyFilename)
{
	diskKickingPerturbations();
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0{params, 10, 0}, pd1{params, 10, 1}, pd2{params, 10, 2};
	std::vector<PotentialDensityPairContainer<KalnajsBasis>> potentialDensityPairs{pd0, pd1, pd2};
	
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};
	std::ofstream out(energyFilename);

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto pd = potentialDensityPairs.begin(); pd != potentialDensityPairs.end(); ++pd){
			for (auto r = radii.begin(); r != radii.end(); ++r){
				std::vector<double> energies = individualEnergyEvolution(*pd, *s, *r);
				outPutEnergy(out, energies, *s, pd -> fourierHarmonic(), *r);
			}
		}
	} 	
	out.close();
	cleaningPerturbations();
}