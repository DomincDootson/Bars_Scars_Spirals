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
const double K_TIME_STANDARD{0.2}, K_ACTIVE_FRACTION{.50}, K_TIME_END{400};
const std::string K_KERNEL{"Kalnajs/Testing/StabilityKernel.out"}, K_PERTURBATION_FILE{"Kalnajs/Testing/perturbation.out"}, G_PERTURBATION_FILE{"GaussianLog/Testing/perturbation.out"};
const std::vector<double> K_PARAMS{4, 20};
const bool K_SELF_CONSISTENT{true};


/* General Functions */

int nStep(double timeStep, double timeEnd = K_TIME_END) 
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

void kalnajsStep(int nMax, double timeStep) { // This function put a value of 0.01 for  1t_{0}
	std::ofstream out(K_PERTURBATION_FILE);

	int nSteps{nStep(timeStep)};
	out << nMax << '\n';

	for (int time = 0; time < 1/timeStep; time += 1) {
		for (int n = 0; n < nMax; ++n) {
			if (n != 0) {out << 0 << " " << 0 <<" ";}
			else {out << 0.01 << " " << 0 <<" ";}
	}

	out << 0 << " " << 0 << '\n';
	}

	for (double time = ((int) 1/timeStep); time < K_TIME_END; time += timeStep){
		for (int n = 0; n < nMax; ++n) {out << 0 << " " << 0 <<" ";}
		out << 0 << " " << 0 << '\n';
	}


}

void saveVector(const std::string & filename, const std::vector<double> & vec) {
	std::ofstream out(filename);
	for (auto i : vec) {out << i << '\n';}
	out.close();
}
std::vector<double> radiiVector() {
	std::vector<double> radii;
 	for (int i =1; i < 2000; ++i) {radii.push_back(i*0.01);}
 	return radii;
 }



void savePerturbation(const std::string & filename, const Eigen::VectorXcd & coeff) {
	ExpansionCoeff holding(nStep(K_TIME_STANDARD), coeff.size()+1);
	holding(0) =  coeff;
	for (int time = 1; time < holding.nTimeStep(); ++time) {holding(time) = 0 * coeff;}
	holding.writePerturbation2File(filename); 
}

void gaussianDelta (int nMax, int index = 15) {
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> kalnajs(params, 10, 2);

	std::vector<double> params1{static_cast<double>(30), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params1, 30, 2);

 	Eigen::VectorXcd coefG = Eigen::VectorXcd::Zero(30+1);
 	coefG(15) = -2;

 	
 	
 
 	Eigen::ArrayXXcd potentialG = PD.potentialArray(coefG, 1500, 20);
 	


 	std::vector<double> params2{static_cast<double>(nMax), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD1(params2, nMax, 2); 

 	Eigen::VectorXcd coefG1 = Eigen::VectorXcd::Zero(nMax+1);

 	potentialG = PD.potentialArray(coefG, 1500, 20);
 	coefG1 = PD1.potentialFitting(potentialG, 20);
 	std::cout << "Saving perturbation to: " << G_PERTURBATION_FILE << '\n';
 	savePerturbation("GaussianLog/Testing/perturbation_" + std::to_string(nMax) +".out", coefG1);
 	
 	
}


/* BF & Kernel Generation Kalnajs */


void generatingKalnajsBF(int m1, int deltaXstep, int nMax, int m2 = 2) {
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, nMax, m2);
	
	Mestel DF;

	ActionAngleBasisContainer test("Kalnajs", nMax, 2, m1, deltaXstep, 20);
	test.scriptW(PD, DF ,"Kalnajs/Testing");
}

void generatingKalnajsKernels(int m1, int deltaXstep, int nMax, double timeStep, int m2 = 2, const std::string & kernelFile = K_KERNEL, const double maxRadius = 20)
{
	ActionAngleBasisContainer test("Kalnajs/Testing", "Kalnajs", nMax, m2, m1, deltaXstep, maxRadius);
	Mestel DF;

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
	//generatingGaussianBF(m1, deltaXstep, nMax, 2);
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

		kalnajsStep(K_NMAX_STANDARD, timeStep);
		

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
	std::vector<int> nMaxValues = {10}; // Please put in ascending order
	
	generatingKalnajsBF(K_M1_STANDARD, K_DELTAXSTEP_STANDARD, nMaxValues.back());

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
	std::vector<int> nMaxValues = {60, 30, 40, 50}; // Please put in ascending order 
	//

	std::ofstream out(filename);

	for (auto nMax : nMaxValues) {
		std::vector<double> params2{static_cast<double>(nMax), .15, 15};
 		PotentialDensityPairContainer<GaussianLogBasis> PD(params2, nMax, 2); 
		//generatingGaussianBF(K_M1_STANDARD, 251, nMax);	

		auto kernelFile = [nMax] () {return "GaussianLog/Testing/Kernel_" + std::to_string(nMax) + ".out"; };
		//generatingGaussianKernels(K_M1_STANDARD, 251, nMax, K_TIME_STANDARD,2, kernelFile()); // THIS ONE
		

		VolterraSolver solver(kernelFile(), nMax, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		gaussianDelta(nMax, K_TIME_STANDARD); // THIS ONE 
		
		std::vector<double> energy = solver.energyEvolution(PD, "GaussianLog/Testing/perturbation_" + std::to_string(nMax) +".out", true);
		energy2File(nMax, energy, out);

		auto file = [nMax] () {return "Plotting/BF_Comparison/Coefficent_Gaussian_" + std::to_string(nMax)+".csv";}; 
		solver.saveResponseCoeff(file()); 
	}
	out.close();
}

void stabilityVaryingCoupling(const std::string & filename) {
	std::vector index{1,2,3,4,5,6,7,8,9, 10, 11};
	std::ofstream out(filename);
	PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, 10, 2);

	kalnajsDelta(10, K_TIME_STANDARD);
	for (auto i : index) {
		std::cout << "Decoupled Index: " << i << '\n'; 
		VolterraSolver solver("Kalnajs/Testing/Kernel_10.out", 10, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		std::vector<double> energy = solver.energyEvolutionDecoupled(PD, K_PERTURBATION_FILE, i, K_SELF_CONSISTENT);
		energy2File(i, energy, out);
	}
	out.close(); 
}

void varyingActiveFraction(const std::string & filename) {
	std::vector<double> xi{0.3, 0.35, 0.4, 0.45, .5};

	std::ofstream out(filename);
	PotentialDensityPairContainer<KalnajsBasis> PD(K_PARAMS, 10, 2);

	kalnajsDelta(10, K_TIME_STANDARD);

		for (auto af : xi) {
		std::cout << "Active Fraction: " << af << '\n'; 
		VolterraSolver solver("Kalnajs/Testing/Kernel_10.out", 10, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(af);

		std::vector<double> energy = solver.energyEvolutionDecoupled(PD, K_PERTURBATION_FILE, K_NMAX_STANDARD, K_SELF_CONSISTENT);
		energy2File(af, energy, out);
	}
	out.close(); 
}


void varyingRka(const std::string & filename) { 
	
	int nMax{10}; 
	std::vector<double> rKa{nMax, 20};

	std::ofstream out(filename);
	kalnajsDelta(nMax, K_TIME_STANDARD);

	for (auto r : rKa) {
		std::vector<double> params{4, r}; Mestel DF;
		PotentialDensityPairContainer<KalnajsBasis> PD(params, nMax, 2);
	
		//ActionAngleBasisContainer test("Kalnajs", nMax, 2, K_M1_STANDARD, K_DELTAXSTEP_STANDARD, r);
		//test.scriptW(PD, DF ,"Kalnajs/Testing");
		
		//std::cout << r/((double) K_DELTAXSTEP_STANDARD-1) << '\n' << '\n';
		//ActionAngleBasisContainer readInTest("Kalnajs/Testing", "Kalnajs", 10, 2, K_M1_STANDARD, K_DELTAXSTEP_STANDARD, r); 

		auto kernelFile = [r] () {return "Kalnajs/Testing/Kernel_10_" + std::to_string(r) + ".out"; };

		//VolterraSolver solver2(nMax, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		//VolterraSolver solver2(nMax, K_M1_STANDARD, 6, K_TIME_STANDARD);

		//solver2.generateKernel(kernelFile(), DF, readInTest);


		VolterraSolver solver(kernelFile(), 10, 2, nStep(K_TIME_STANDARD), K_TIME_STANDARD);
		solver.activeFraction(K_ACTIVE_FRACTION);

		std::vector<double> energy = solver.energyEvolutionDecoupled(PD, K_PERTURBATION_FILE, 10, K_SELF_CONSISTENT);

		energy2File(r, energy, out);
	}
	out.close(); 
}

int main() {
	//stabilityM1Kalanajs("Plotting/KalnajsTesting/kalnajsVaryingm1.csv");
	//stabilityTimeStep("Plotting/KalnajsTesting/kalnajsVaryingTimestep.csv");
	//stabilityDeltaXGridKalnajs("Plotting/KalnajsTesting/kalnajsVaryingGridSize.csv");
	//stabilityNBasisFunctions("Plotting/KalnajsTesting/kalnajsVaryingNbasis.csv");
	//stabilityVaryingCoupling("Plotting/KalnajsTesting/kalnajsVaryingCoupling.csv");
	//stabilityTimeStep("Plotting/KalnajsTesting/VaryingTimestep.csv");
	//stabilityNBasisFunctionsGaussian("Plotting/GaussianTesting/gaussianVaryingNbasis.csv");
	//generatingKalnajsKernels(7, K_DELTAXSTEP_STANDARD, 0, K_TIME_STANDARD, 2, "Plotting/Kalnaj_0.out"); 
	//stabilityNBasisFunctions("Plotting/BF_Comparison/kalnajsVaryingNbasisSmallStep.csv");
	//varyingActiveFraction("Plotting/KalnajsTesting/VaryingXI.csv");

	varyingRka("Plotting/KalnajsTesting/varyingRka.csv");
	//varyingModePoke("Plotting/BF_Comparison/kalnajsVaryingNpoke.csv");
	//stabilityNBasisFunctionsGaussian("Plotting/BF_Comparison/gaussianVaryingNbasis.csv");
	return 0;
}