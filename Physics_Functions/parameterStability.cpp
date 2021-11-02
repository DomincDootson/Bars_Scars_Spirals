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

#include "../DF_Function/DFfunction.h"

const int K_M1_STANDARD{7}, K_DELTAXSTEP_STANDARD{101}, K_NMAX_STANDARD{10};
const double K_TIME_STANDARD{0.2}, K_ACTIVE_FRACTION{0.43};
const std::string K_KERNEL{"Kalnajs/Testing/StabilityKernel.out"}, K_PERTURBATION_FILE{"Kalnajs/Testing/perturbation.out"};
const std::vector<double> K_PARAMS{4, 20};
const bool K_SELF_CONSISTENT{true};


/* General Functions */

int nStep(double timeStep, double timeEnd = 400) 
{
	return ((int) round(timeEnd/timeStep));
}

// would it be possible to add the value of our variation to the output? 

void energy2File(double value, std::vector<double> & energy, std::ofstream & out){
	out << value << ',';
	for (int i = 0; i < energy.size()-1; ++i) {out << energy[i] << ',';}
	out << energy.back() <<'\n';
}

// function that writes the peurbation to file 

void kalnajsPerturbation(int nMax, double timeStep, double timeEnd = 50) {
	std::ofstream out(K_PERTURBATION_FILE);

	int nSteps{nStep(timeStep, timeEnd)};
	out << nMax <<'\n'; 
	for (int i = 0; i < nSteps; ++i) {
		out << 0.01 * sin(M_PI * ((i*timeStep)/20)) << " " << 0 << " ";
		for (int n = 1; n < nMax; ++n) {out << 0 << " " << 0 <<" ";}
		out << 0 << " " << 0 << '\n';
	}
	
	out.close(); 
}

void kalnajsDelta(int nMax, double timeStep, double timeEnd = 50, int nMode = 0) {
	std::ofstream out(K_PERTURBATION_FILE);

	int nSteps{nStep(timeStep, timeEnd)};
	out << nMax << '\n';
	for (int n = 0; n < nMax; ++n) {
		if (n != nMode) {out << 0 << " " << 0 <<" ";}
		else {out << 0.01 << " " << 0 <<" ";}
	}

	if (nMode == nMax) {out << 0.01 << " " << 0 << '\n';}
	else {out << 0 << " " << 0 << '\n';}
	
	for (int i = 1; i < nSteps; ++i) {
		for (int n = 0; n < nMax; ++n) {out << 0 << " " << 0 <<" ";}
		out << 0 << " " << 0 << '\n';
	}
	
	out.close(); 
}

void gaussianDelta(int nMax, double timeStep, double timeEnd = 50) {
	std::ofstream out(K_PERTURBATION_FILE);

	int nSteps{nStep(timeStep, timeEnd)};
	out << nMax <<'\n'; 

	out << 0.0 << " " << 0 << " ";
	//out << 0.01 << " " << 0 << " ";
	for (int n = 0; n < nMax; ++n) {
		if (n == 30) {out << 0.01 << " " << 0 <<" ";}
		else {out << 0 << " " << 0 <<" ";}
	} 
	

	out << 0.0 << " " << 0 << '\n';
	for (int i = 1; i < nSteps; ++i) {
		for (int n = 0; n < nMax; ++n) {out << 0 << " " << 0 <<" ";}
		out << 0 << " " << 0 << '\n';
	}
	
	out.close(); 
}


/* BF & Kernel Generation Kalnajs */


void generatingKalnajsBF(int m1, int deltaXstep, int nMax, int m2 = 2) {
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, nMax, m2);
	
	Mestel DF(1, 1, 0.61);

	ActionAngleBasisContainer test("Kalnajs", nMax, 2, m1, deltaXstep, 20);
	test.scriptW(PD, DF ,"Kalnajs/Testing");
}

void generatingKalnajsKernels(int m1, int deltaXstep, int nMax, double timeStep, int m2 = 2, const std::string & kernelFile = K_KERNEL)
{
	

	ActionAngleBasisContainer test("Kalnajs/Testing", "Kalnajs", nMax, m2, m1, deltaXstep, 20);
	Mestel DF(1, 1, 0.35);

	VolterraSolver solver2(nMax, m2, nStep(timeStep), timeStep);

	solver2.generateKernel(kernelFile, DF, test);
}


void generatingGaussianBF(int m1, int deltaXstep, int nMax, int m2 = 2) {
	std::vector<double> params{static_cast<double>(nMax), .15, 15};
	PotentialDensityPairContainer<GaussianLogBasis> PD(params, nMax, m2);
	
	Mestel DF(1, 1, 0.35);

	ActionAngleBasisContainer test("GaussianLog", nMax, 2, m1, deltaXstep, 20);
	test.scriptW(PD, DF ,"GaussianLog/Testing");
}

void generatingGaussianKernels(int m1, int deltaXstep, int nMax, double timeStep, int m2 = 2, const std::string & kernelFile = K_KERNEL)
{
	generatingGaussianBF(m1, deltaXstep, nMax, 2);
	ActionAngleBasisContainer test("GaussianLog/Testing", "GaussianLog", nMax, m2, m1, deltaXstep, 20);
	Mestel DF(1, 1, 0.35);

	VolterraSolver solver2(nMax, m2, nStep(timeStep), timeStep);

	solver2.generateKernel(kernelFile, DF, test);
}


void saveScriptE(const Eigen::MatrixXd & scriptE, int m2 = 2) {
	std::ofstream out("Kalnajs/Testing/scriptE_" +std::to_string(scriptE.rows()) + "_" +std::to_string(m2) +".out");

	for (int i = 0; i < scriptE.rows(); ++i)
	{
		for (int j = 0; j < scriptE.cols(); ++j) {
			if (i == (scriptE.cols()-1)) {out << scriptE(i,j) << '\n';}
			else if ((i == (scriptE.cols()-1)) && (j == (scriptE.rows()-1))) {out << scriptE(i,j);}
			else {out << scriptE(i,j);}
		}
	}
	out.close(); 
} 

/* Kalnajs Testing */ 


void stabilityM1Kalanajs(const std::string & filename) {
	std::cout << "Checking stability w.r.t. number of m1 summed over.\n\n";
	std::vector<int> m1Vector = {2, 5, 6, 7}; // Please put in ascending order
	generatingKalnajsBF(m1Vector.back(), K_DELTAXSTEP_STANDARD, K_NMAX_STANDARD);

	std::ofstream out(filename);

	for (auto m1 : m1Vector) {
		PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, K_NMAX_STANDARD, 2);

		generatingKalnajsKernels(m1, K_DELTAXSTEP_STANDARD, K_NMAX_STANDARD, K_TIME_STANDARD);

		VolterraSolver solver(K_KERNEL, K_NMAX_STANDARD, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		kalnajsPerturbation(K_NMAX_STANDARD, K_TIME_STANDARD);

		std::vector<double> energy = solver.energyEvolution(PD, K_PERTURBATION_FILE, K_SELF_CONSISTENT);
		energy2File(m1, energy, out);
	}
	out.close();
}

void stabilityTimeStep(const std::string & filename) {
	std::cout << "Checking stability w.r.t. time step.\n\n";
	std::vector<double> timeSteps = {0.1, 0.2, 0.5, 1}; // Please put in ascending order
	//generatingKalnajsBF(K_M1_STANDARD, K_DELTAXSTEP_STANDARD, K_NMAX_STANDARD);

	std::ofstream out(filename);

	for (auto timeStep : timeSteps) {
		PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, K_NMAX_STANDARD, 2);

		generatingKalnajsKernels(K_M1_STANDARD, K_DELTAXSTEP_STANDARD, K_NMAX_STANDARD, timeStep);

		VolterraSolver solver(K_KERNEL, K_NMAX_STANDARD, 2, nStep(timeStep), timeStep);
		solver.activeFraction(K_ACTIVE_FRACTION);

		kalnajsPerturbation(K_NMAX_STANDARD, timeStep);

		std::vector<double> energy = solver.energyEvolution(PD, K_PERTURBATION_FILE, K_SELF_CONSISTENT);
		energy2File(timeStep, energy, out);
	}
	out.close();
}

void stabilityDeltaXGridKalnajs(const std::string & filename) {
	std::cout << "Checking stability w.r.t. grid size.\n\n";
	std::vector<int> deltaXSteps = {81, 101, 151, 201}; // Please put in ascending order
	
	std::ofstream out(filename);

	for (auto deltaXstep : deltaXSteps) {
		generatingKalnajsBF(K_M1_STANDARD, deltaXstep, K_NMAX_STANDARD);

		PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, K_NMAX_STANDARD, 2);

		generatingKalnajsKernels(K_M1_STANDARD, deltaXstep, K_NMAX_STANDARD, K_TIME_STANDARD);

		VolterraSolver solver(K_KERNEL, K_NMAX_STANDARD, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		kalnajsPerturbation(K_NMAX_STANDARD, K_TIME_STANDARD);

		std::vector<double> energy = solver.energyEvolution(PD, K_PERTURBATION_FILE, K_SELF_CONSISTENT);
		energy2File(deltaXstep, energy, out);
	}
	out.close();
}


void stabilityNBasisFunctions(const std::string & filename){
	std::cout << "Checking stability w.r.t. number of basis functions.\n\n";
	std::vector<int> nMaxValues = {1,2,3,4,5,6,7, 8, 9, 10}; // Please put in ascending order
	//generatingKalnajsBF(K_M1_STANDARD, K_DELTAXSTEP_STANDARD, nMaxValues.back());

	std::ofstream out(filename);

	for (auto nMax : nMaxValues) {
		PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, nMax, 2);
		//saveScriptE(PD.scriptE()); 
		auto kernelFile = [nMax] () {return "Kalnajs/Testing/Kernel_" + std::to_string(nMax) + ".out"; };
		generatingKalnajsKernels(K_M1_STANDARD, K_DELTAXSTEP_STANDARD, nMax, K_TIME_STANDARD,2, kernelFile());

		VolterraSolver solver(kernelFile(), nMax, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		
		//kalnajsPerturbation(nMax, K_TIME_STANDARD);
		kalnajsDelta(nMax, K_TIME_STANDARD);

		std::vector<double> energy = solver.energyEvolution(PD, K_PERTURBATION_FILE, K_SELF_CONSISTENT);
		energy2File(nMax, energy, out);

		auto file = [nMax] () {return "Plotting/BF_Comparison/Coefficent_" + std::to_string(nMax)+".csv";}; 
		solver.saveResponseCoeff(file()); 
	}
	out.close();
}

void varyingModePoke(const std::string & filename){
	std::vector<int> nMaxValues = {0,1,2,3,4,5, 6, 7, 8, 9, 10}; // Please put in ascending order
	//generatingKalnajsKernels(K_M1_STANDARD, K_DELTAXSTEP_STANDARD, 10, K_TIME_STANDARD, 2, "Kalnajs/Testing/Kernel_10.out");

	std::ofstream out(filename);

	for (auto nMode : nMaxValues) {
 		int nMax{10}; PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, nMax, 2);
		auto kernelFile = [nMax] () {return "Kalnajs/Testing/Kernel_" + std::to_string(nMax) + ".out"; };

		VolterraSolver solver(kernelFile(), nMax, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		kalnajsDelta(nMax, K_TIME_STANDARD, 400, nMode);

		std::vector<double> energy = solver.energyEvolution(PD, K_PERTURBATION_FILE, K_SELF_CONSISTENT);
		energy2File(nMax, energy, out);

		auto file = [nMode] () {return "Plotting/BF_Comparison/Coefficent_" + std::to_string(nMode)+".csv";}; 
		solver.saveResponseCoeff(file()); 
	}
	out.close();
}


void stabilityNBasisFunctionsGaussian(const std::string & filename){
	std::cout << "Checking stability w.r.t. number of basis functions.\n\n";
	std::vector<int> nMaxValues = {20, 50, 60, 70, 80}; // Please put in ascending order
	//generatingKalnajsBF(K_M1_STANDARD, K_DELTAXSTEP_STANDARD, nMaxValues.back());

	std::ofstream out(filename);

	for (auto nMax : nMaxValues) {
		PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, nMax, 2);
		//saveScriptE(PD.scriptE()); 
		auto kernelFile = [nMax] () {return "GaussianLog/Testing/Kernel_" + std::to_string(nMax) + ".out"; };
		generatingGaussianKernels(K_M1_STANDARD, 251, nMax, K_TIME_STANDARD,2, kernelFile()); // THIS ONE

		VolterraSolver solver(kernelFile(), nMax, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		//solver.saveKernelForPython("Plotting/kernel.csv"); 
		//kalnajsPerturbation(nMax, K_TIME_STANDARD);
		gaussianDelta(nMax, K_TIME_STANDARD); // THIS ONE 

		std::vector<double> energy = solver.energyEvolution(PD, K_PERTURBATION_FILE, K_SELF_CONSISTENT);
		energy2File(nMax, energy, out);

		auto file = [nMax] () {return "Plotting/BF_Comparison/Coefficent_Gaussian_" + std::to_string(nMax)+".csv";}; 
		solver.saveResponseCoeff(file()); 
	}
	out.close();
}


int main() {
	//stabilityM1Kalanajs("Plotting/KalnajsTesting/kalnajsVaryingm1.csv");
	//stabilityTimeStep("Plotting/KalnajsTesting/kalnajsVaryingTimestep.csv");
	//stabilityDeltaXGridKalnajs("Plotting/KalnajsTesting/kalnajsVaryingGridSize.csv");
	stabilityNBasisFunctions("Plotting/KalnajsTesting/kalnajsVaryingNbasis.csv");


	//stabilityNBasisFunctions("Plotting/BF_Comparison/kalnajsVaryingNbasis.csv");
	//varyingModePoke("Plotting/BF_Comparison/kalnajsVaryingNpoke.csv");
	//stabilityNBasisFunctionsGaussian("Plotting/BF_Comparison/gaussianVaryingNbasis.csv");
	return 0;
}