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

#include "barEvolutionFunctions.h"


std::string kernelName(std::string dir, std::string stem, int Kka, int Rka, int m2)
{
	return dir + "/" + stem + "_" + std::to_string(Kka) + "_" + std::to_string(Rka) + "_" + std::to_string(m2) + ".out";
}

std::string kernelName(std::string dir, std::string stem, int Kka, double Rka, int m2) 
{
	return dir + "/" + stem + "_" + std::to_string(Kka) + "_" + std::to_string((int) round(Rka)) + "_" + std::to_string(m2) + ".out";
}


void barGrowthRate(const std::string & filename, double timeScale) {
	std::ofstream out(filename);
	out << 100 << '\n';
	double step{timeScale/100};
	for (double time = 0; time < timeScale; time += step) {out << time << " " << pow(sin((time /timeScale)*0.5*M_PI),2) << '\n';}
	out.close();
}

void kalnajsBarTest() {


	VolterraSolver solver("test200.out", 10, 2, 200, 0.1);

	solver.activeFraction(.25);	
	
	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	coef[0] = .01; 
	
	Bar2D bar(coef, 0.5);
	solver.barRotation(bar, "linearCoeff.csv", "linearEvolution.csv", false, false); 
	
}

void kalnajBFVaryingK()
{
	Mestel DF;
	std::vector<double> Kka{4, 5, 6, 7};
	for (int i = 0; i < 4; ++i){
		
		std::cout << "Calculationg BF for Kka: " << Kka[i] << '\n';

		std::vector<double> params{Kka[i], 10};
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2);

		ActionAngleBasisContainer test("Kalnajs", 10, 2, 5, 101, 10);
		std::string file = "Kalnajs/Kalnajs_" + std::to_string((int) Kka[i]) + "_10"; // Use file function name here
		test.scriptW(PD, DF, file);
	}
}

void kalnajBFVaryingR()
{
	Mestel DF;
	std::vector<double> Rka{5, 10, 15, 20};
	for (int i = 0; i < 4; ++i){
		
		
		std::cout << "Calculationg BF for Rka: " << Rka[i] << '\n';
		std::vector<double> params{4, Rka[i]};
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2);

		ActionAngleBasisContainer test("Kalnajs", 10, 2, 5, 101, 10);
		std::string file = "Kalnajs/Kalnajs_4_" + std::to_string((int) Rka[i]); // Use file function name here
		test.scriptW(PD, DF, file);
	}
}


void kalnajBasisFunctionsVaryingK() 
{
	std::vector<double> params{4, 10};
	std::ofstream out("quad_Density_Functions.csv");

	for (int i = 0; i<999; ++i)
	{
		out << i*0.01 << ',';
	}
	out << 999*0.01 << '\n';

	for (int k = 3; k <9;++k)
	{
		params[0] = k;
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2); // Remake the plots
		for (int i = 0; i<999; ++i)
		{
			out << PD.density(i*0.01,0) << ',';
		}
		out << PD.density(999*0.01,0)  << '\n';

	}
	out.close();
}

std::string evolutionFileName(std::string dir, double omega0, bool isSelfConsistent = true)
{
	std::string filename = dir + "/evolution"; 

	if (!isSelfConsistent){filename += "Test";}

	filename += std::to_string((int) round(100*omega0)) + ".csv";

	return filename;
}

std::string coeffFileName(std::string dir, double omega0, bool isSelfConsistent = true)
{
	std::string filename = dir + "/coeff";
	if (!isSelfConsistent) {filename += "Test";}

	filename += std::to_string((int) round(100*omega0)) + ".csv";		
	return filename;
}
Eigen::VectorXcd kalnajsResolving(double radius = 0, int nMax = 10) 
{
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff(0) = 0.01; 
	return coeff;
}

void kalnajsTorque(int nMax) {
	double barRadius{2}, timeStep{0.25};
	int numbTimeSteps{800}; 

	std::string kernel1 = "Kernels/Kalnajs100.out";
	 
	ActionAngleBasisContainer test("Kalnajs/Kalnajs_4_20", "Kalnajs", nMax, 2, 5, 101, 20);
	Mestel DF(1,1, .35);

	VolterraSolver solver1(nMax, 2, numbTimeSteps, timeStep);
	solver1.generateKernel(kernel1, DF, test);
	solver1.saveKernelForPython("Plotting/kernel.csv");
	
	// Produce kernel
	std::vector<double> growthRate = {10};

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd(params, nMax, 2);
	
	bool selfConsistent{true};
	for (auto rate : growthRate) {
		barGrowthRate("Bar2D/barSizeDiff.out", rate);

		Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff =  kalnajsResolving(barRadius, nMax);
		Bar2D bar(coeff, .1, "Bar2D/barSizeDiff.out");

		VolterraSolver solver(kernel1, nMax, 2, numbTimeSteps, timeStep);
		solver.activeFraction(.45);
		solver.barRotation(bar, "Plotting/KalnajsTorque/VaryingN/Coeff_" + std::to_string(nMax)+ ".csv", "Plotting/KalnajsTorque/VaryingN/Evolution_"+ std::to_string(nMax)+".csv", selfConsistent, false, true);
		//solver.writeDensity2File("Plotting/GaussianTorque/densityEvolution.csv", pd);
		//exit(0);
	}

}

void barVaryingAngularSpeed()
{
	std::vector<double> angSpeed = {.04, .08, .12, .16, .20};

	std::string kernelFileName = kernelName("Kernels", "Kalnajs", 4, 10, 2);
	VolterraSolver solver(kernelFileName, 10, 2, 2000, 0.025);

	solver.activeFraction(.25);	
	
	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	coef[0] = .1; 
	
	for (int i =0; i < angSpeed.size(); ++i){
		// Construct the bar
		Bar2D bar(coef, angSpeed[i]);
		std::string evolutionFilename = evolutionFileName("BarEvolution/VaryingOmega", angSpeed[i]);
		std::string outFilename = coeffFileName("BarEvolution/VaryingOmega", angSpeed[i]);
		solver.barRotation(bar, outFilename, evolutionFilename, true, true); 
	}
}

void barVaryingKka()
{
	std::vector<int> Kka{4, 5, 6, 7};
	for (int i = 0; i < 4; i++)
	{
		std::cout << "Evolution for Kka: " << Kka[i] << '\n';
		std::string kernelFileName = kernelName("Kernels", "Kalnajs", Kka[i], 10, 2);
		VolterraSolver solver(kernelFileName, 10, 2, 2000, 0.025);
		solver.activeFraction(.25);	
		 
		Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
		coef[0] = .1;
		Bar2D bar(coef, 0);

		std::string outFilename = coeffFileName("BarEvolution/VaryingKka", 0.01*Kka[i]);
		std::string evolutionFilename = evolutionFileName("BarEvolution/VaryingKka", 0.01*Kka[i]);

		solver.barRotation(bar, outFilename, evolutionFilename); 
	}
}

void barVaryingRka()
{
	std::vector<double> Rka{5, 10, 15, 20};
	for (int i = 0; i < 4; i++)
	{
		std::cout << "Evolution for Rka: " << Rka[i] << '\n';
		std::string kernelFileName = kernelName("Kernels", "Kalnajs", 4, Rka[i], 2);
		VolterraSolver solver(kernelFileName, 10, 2, 2000, 0.025);
		solver.activeFraction(.25);	
		 
		Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
		coef[0] = .1;
		Bar2D bar(coef, 0);

		std::string outFilename = coeffFileName("BarEvolution/VaryingRka", 0.01*Rka[i]);
		std::string evolutionFilename = evolutionFileName("BarEvolution/VaryingRka", 0.01*Rka[i]);

		solver.barRotation(bar, outFilename, evolutionFilename); 
	}
}

void barVaryingTurnOn() 
{
	double barRadius{2}, timeStep{0.25};
	int nMax{10}, numbTimeSteps{800}; 

	std::string kernel1 = "Kernels/Kalnajs100.out";
	 
	ActionAngleBasisContainer test("Kalnajs/Kalnajs_4_20", "Kalnajs", nMax, 2, 5, 101, 20);
	Mestel DF(1,1, .35);

	VolterraSolver solver1(nMax, 2, numbTimeSteps, timeStep);
	//solver1.generateKernel(kernel1, DF, test);
	
	
	// Produce kernel
	std::vector<double> growthRate = {5, 10, 15, 20, 30};

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd(params, nMax, 2);
	
	bool selfConsistent{false};
	for (auto rate : growthRate) {
		barGrowthRate("Bar2D/barSizeDiff.out", rate);

		Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff =  kalnajsResolving(barRadius, nMax);
		Bar2D bar(coeff, .1, "Bar2D/barSizeDiff.out");

		VolterraSolver solver(kernel1, nMax, 2, numbTimeSteps, timeStep);
		solver.activeFraction(.4);
		solver.barRotation(bar, coeffFileName("Plotting/KalnajsTorque", rate, selfConsistent), evolutionFileName("Plotting/KalnajsTorque", rate, selfConsistent), selfConsistent, false, true);
		//solver.writeDensity2File("Plotting/GaussianTorque/densityEvolution.csv", pd);
		//exit(0);
	}
} 

void barVaryingActiveFraction() {
	double barRadius{2}, timeStep{0.25};
	int nMax{10}, numbTimeSteps{2*800}; 

	std::string kernel1 = "Kernels/Kalnajs100.out";
	 
	ActionAngleBasisContainer test("Kalnajs/Kalnajs_4_20", "Kalnajs", nMax, 2, 5, 101, 20);
	Mestel DF(1,1, .35);

	VolterraSolver solver1(nMax, 2, numbTimeSteps, timeStep);
	//solver1.generateKernel(kernel1, DF, test);
	
	
	// Produce kernel
	std::vector<double> activeFractions = {.40, .42, .44, .46, .48, .5};

	std::vector<double> params{static_cast<double>(nMax), 4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd(params, nMax, 2);
	
	bool selfConsistent{true};
	for (auto xi : activeFractions) {
		barGrowthRate("Bar2D/barSizeDiff.out", 5);

		Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff =  kalnajsResolving(barRadius, nMax);
		Bar2D bar(coeff, .1, "Bar2D/barSizeDiff.out");

		VolterraSolver solver(kernel1, nMax, 2, numbTimeSteps, timeStep);
		solver.activeFraction(xi);
		solver.barRotation(bar, coeffFileName("Plotting/KalnajsTorque", xi, selfConsistent), evolutionFileName("Plotting/KalnajsTorque", xi, selfConsistent), selfConsistent, false, true);
		//solver.writeDensity2File("Plotting/GaussianTorque/densityEvolution.csv", pd);
		//exit(0);
	}
}

void longTermEvolution() {
	double barRadius{2}, timeStep{0.25};
	int nMax{10}, numbTimeSteps{2*800}; 

	std::string kernel1 = "Kernels/Kalnajs100.out";

	VolterraSolver solver(kernel1, nMax, 2, numbTimeSteps, timeStep);
	solver.activeFraction(.4);

	Bar2D bar(kalnajsResolving(2), .1, "Bar2D/barSizeDiff.out");
	std::cout << bar.torque(solver.longTimeCoeff(bar.barCoeff(), false)) << '\n'; 
	std::cout << bar.torque(solver.longTimeCoeff(bar.barCoeff(), true)) << '\n'; 
}

void kalnajsKernelsDiffTemp()
{
	double timeStep{0.25};
	int nMax{10}, numbTimeSteps{800};

	std::vector sigmas = {.35, .45};
	ActionAngleBasisContainer test("Kalnajs/Kalnajs_4_20", "Kalnajs", nMax, 2, 5, 101, 20);

	for (auto it = sigmas.begin(); it != sigmas.end(); ++it) { 
		Mestel DF(1, 1, *it);
		
		VolterraSolver solver(nMax, 2, numbTimeSteps, timeStep);

		std::string kernel = "Kernels/Kalnajs_2_" + std::to_string((int) round(100*(*it))) + ".out";
		//solver.generateKernel(kernel, DF, test); 
	}

	bool selfConsistent{true};
	for (auto it: sigmas) {
		std::string kernel = "Kernels/Kalnajs_2_" + std::to_string((int) round(100*(it))) + ".out";
		VolterraSolver solver(kernel, nMax, 2, numbTimeSteps, timeStep);
		solver.activeFraction(.4);

		Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff =  kalnajsResolving();
		Bar2D bar(coeff, .1, "Bar2D/barSize.out");
 
		solver.barRotation(bar, coeffFileName("Plotting/KalnajsTorque/DiffTemp", it, selfConsistent), evolutionFileName("Plotting/KalnajsTorque/DiffTemp", it, selfConsistent), selfConsistent, false, true);
	}
}


// Gaussian Bar //

Eigen::VectorXcd gaussianResolving(double radius = 2, int nMax = 80)
{
	std::vector<double> params{static_cast<double>(nMax), .15, 15};
	PotentialDensityPairContainer<GaussianLogBasis> pd(params, nMax, 2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(pd.maxRadialIndex() + 1);
	for (int i = 0; i < pd.maxRadialIndex() + 1; ++i){
		coeff(i) = 0.01*pd.potential(radius, i);
	}
	return - (pd.scriptE()).inverse() * coeff;
}

Eigen::VectorXcd kalnajsResolving()
{
	int nMax{10};
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd(params, nMax, 2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(pd.maxRadialIndex() + 1);
	for (int i = 0; i < pd.maxRadialIndex() + 1; ++i){
		coeff(i) = 0.01*pd.potential(2.06271, i);
	}
	return - (pd.scriptE()).inverse() * coeff;
}


void gaussianBarEvolution(){
	double barRadius{1.977}, timeStep{0.25}, angSpeed{.1};
	int nMax{50}, numbTimeSteps{200};
	bool selfConsistent{false};
	
	Bar2D barCold(gaussianResolving(barRadius), angSpeed, "Bar2D/barSize.out");
	VolterraSolver solverCold( "GaussianLog_Cold.out", nMax, 2, numbTimeSteps, timeStep);
	solverCold.activeFraction(.25);
	solverCold.barRotation(barCold, "Plotting/barCoefficents.csv", "Plotting/barEvolutionCold_Test.csv", selfConsistent, false);

	Bar2D barMed(gaussianResolving(barRadius), angSpeed, "Bar2D/barSize.out");
	VolterraSolver solverMed( "GaussianLog_Med.out", nMax, 2, numbTimeSteps, timeStep);
	solverMed.activeFraction(.25);
	solverMed.barRotation(barMed, "Plotting/barCoefficents.csv", "Plotting/barEvolutionMed_Test.csv", selfConsistent, false);

	Bar2D barWarm(gaussianResolving(barRadius), angSpeed, "Bar2D/barSize.out");
	VolterraSolver solverWarm( "GaussianLog_Warm.out", nMax, 2, numbTimeSteps, timeStep);
	solverWarm.activeFraction(.25);
	solverWarm.barRotation(barWarm, "Plotting/barCoefficents.csv", "Plotting/barEvolutionWarm_Test.csv", selfConsistent, false);
}

void gaussianKernelsDiffTemp()
{
	double barRadius{2}, timeStep{0.1};
	int nMax{24}, numbTimeSteps{500};

	std::vector sigmas = {.35, .45};
	ActionAngleBasisContainer test("GaussianLog", "GaussianLog", 24, 2, 10, 301, 20);

	for (auto it = sigmas.begin(); it != sigmas.end(); ++it) { 
		Mestel DF(1, 1, *it);
		
		VolterraSolver solver(24, 2, numbTimeSteps, timeStep);

		std::string kernel = "Kernels/GaussianLog_2_" + std::to_string((int) round(100*(*it))) + ".out";
		solver.generateKernel(kernel, DF, test); 
	}
}

void gaussianEvolutionDiffTemp() {
	double barRadius{2}, timeStep{0.1};
	int nMax{24}, numbTimeSteps{500};

	std::vector sigmas = {.35, .45};

	for (auto it = sigmas.begin(); it != sigmas.end(); ++it) { 
		std::string kernel = "Kernels/GaussianLog_2_" + std::to_string((int) round(100*(*it))) + ".out";

		Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(24+1); coeff =  gaussianResolving();
		Bar2D bar(coeff, .1, "Bar2D/barSize.out");
 

		VolterraSolver solver(kernel, 24, 2, 500, 0.1);
		solver.activeFraction(.25);
		solver.barRotation(bar, coeffFileName("Plotting/GaussianTorque", *it), evolutionFileName("Plotting/GaussianTorque", *it), false, false, true);
	}
}

void gaussianBarDiffGrowthRate() {
	double barRadius{2}, timeStep{0.25};
	int nMax{80}, numbTimeSteps{400};

	//ActionAngleBasisContainer test("GaussianLog", "GaussianLog", nMax, 2, 7, 251, 20);
	//Mestel DF(1,1, .35);

	//VolterraSolver solver1(nMax, 2, numbTimeSteps, timeStep);
	std::string kernel1 = "Kernels/Gaussian100.out";
	//solver1.generateKernel(kernel1, DF, test);


	// Produce kernel
	std::vector<double> growthRate = {5, 10, 15, 20, 30};

	std::vector<double> params{static_cast<double>(nMax), .15, 15};
	PotentialDensityPairContainer<GaussianLogBasis> pd(params, nMax, 2);
	
	for (auto rate : growthRate) {
		barGrowthRate("Bar2D/barSizeDiff.out", rate);

		Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff =  gaussianResolving(barRadius, nMax);
		Bar2D bar(coeff, .1, "Bar2D/barSizeDiff.out");

		VolterraSolver solver(kernel1, nMax, 2, numbTimeSteps, timeStep);
		solver.activeFraction(.25);
		solver.barRotation(bar, coeffFileName("Plotting/GaussianTorque", rate), evolutionFileName("Plotting/GaussianTorque", rate), false, false, true);
		//solver.writeDensity2File("Plotting/GaussianTorque/densityEvolution.csv", pd);
		//exit(0);
	}
}


 