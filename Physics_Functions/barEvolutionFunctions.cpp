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

	filename += "_Sormani_" + std::to_string((int) round(100*omega0)) + ".csv";

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


void kalnajsKernelsDiffTemp()
{
	double timeStep{0.25};
	int nMax{40}, numbTimeSteps{1200};

	std::vector sigmas = {.35, .45};
	//ActionAngleBasisContainer test("KalnajsN", "KalnajsN", nMax, 2, 4, 251, 15);

	for (auto it = sigmas.begin(); it != sigmas.end(); ++it) { 
		Mestel DF(1, 1, *it);
		
		VolterraSolver solver(nMax, 2, numbTimeSteps, timeStep);

		std::string kernel = "Kernels/Kalnajs_2_" + std::to_string((int) round(100*(*it))) + ".out";
		//solver.generateKernel(kernel, DF, test); 
	}
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat",40);
	bool selfConsistent{true};
	for (auto it: sigmas) {
		std::string kernel = "Kernels/Kalnajs_2_" + std::to_string((int) round(100*(it))) + ".out";
		VolterraSolver solver(kernel, nMax, 2, numbTimeSteps, timeStep);
		solver.activeFraction(.40);

		Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff(0) = 0.01; 
		Bar2D bar(coeff, .09, "Bar2D/barSize.out");
		
		//Bar2D bar(0.05, "Bar2D/barSize.out");

		bar.sormaniBar(pd); 
 
		solver.barRotation(bar, coeffFileName("Plotting/KalnajsTorque/Monari_Bar", it, selfConsistent), evolutionFileName("Plotting/KalnajsTorque/Monari_Bar", it, selfConsistent), selfConsistent, false, true);
	}
}

/* Sormani Plotting Bar */
/* -------------------- */  

void differentTempKernels() {
	double timeStep{1};
	int nMax{48}, numbTimeSteps{100};

	std::vector sigmas = {.35};
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", nMax, 2, 10, 251, 15);

	for (auto it = sigmas.begin(); it != sigmas.end(); ++it) { 
		Mestel DF(1, 1, *it, 1, 1);
		
		VolterraSolver solver(nMax, 2, numbTimeSteps, timeStep);

		std::string kernel = "Kernels/Kalnajs_2_" + std::to_string((int) round(100*(*it))) + "_New.out";
		solver.generateKernel(kernel, DF, test); 
	}
} 

void saveFittedSormani(const std::string & filename, const std::string & barFile) {
	int nMax{72};


	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", nMax);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); coeff(0) = 0.01; 
	Bar2D bar(coeff, .09, "Bar2D/barSize.out");

	bar.sormaniBar(pd, 1, barFile); 
	double step{10/400.0};
	std::vector<double> radii;
	for (double r = 0; r < 15; r += step) {radii.push_back(r);}

	std::vector<double> pot_fit = pd.oneDpotential(radii, bar.barCoeff());
	std::ofstream out(filename); 
	for (int i = 0; i < radii.size(); ++i) {out << radii[i] <<',' << pot_fit[i] <<',' << bar.sormaniPotential(radii[i], 99)<< '\n'; std::cout << radii[i] <<'\n';}
	out.close();
}


void sormaniBarRun(const double temp, const bool selfConsistent, const std::string & barFile = "Bar2D/Bar_Potentials/Sormani_Medium.out") {
	double timeStep{1};
	int nMax{48}, numbTimeSteps{100};

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", nMax);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); 
	coeff(0) = 0.01; 

	Bar2D bar(coeff, 0.18, "Bar2D/barSize.out"); 

	bar.sormaniBar(pd, 0.05, barFile); 


	std::string kernel = "Kernels/Kalnajs_2_" + std::to_string((int) round(100*(temp))) + "_New.out";
	VolterraSolver solver(kernel, nMax, 2, numbTimeSteps, timeStep); 
	solver.activeFraction(.50);
	
 	
 	solver.barRotation(bar, coeffFileName("Plotting/KalnajsTorque/Sormani_Bar", temp, selfConsistent), evolutionFileName("Plotting/KalnajsTorque/Sormani_Bar", temp, selfConsistent), selfConsistent, false, true);
 	//solver.barRotation(bar, "Plotting/test.csv", evolutionFileName("Plotting/KalnajsTorque/Sormani_Bar", temp, selfConsistent), selfConsistent, false, true);
 	//solver.barRotationUnsaving(bar, selfConsistent, false, true);
 	//solver.density2dEvolution(evolutionFileName("Plotting/Bar_Data/Sormani_Diff_Temp", temp, selfConsistent), pd, 25, 5);
 	//solver.potential2dEvolution(evolutionFileName("Plotting/Bar_Data/Sormani_Diff_Temp", temp, selfConsistent), pd, 25, 5);
}

void sormaniBarEvolution(const std::string & barFile) {
	//sormaniBarRun(0.35, true, barFile);
	sormaniBarRun(0.35, false, barFile);

	//sormaniBarRun(0.45, true, barFile);
	//sormaniBarRun(0.45, false, barFile);
}

/* Testing different Soramni Effects */
/* --------------------------------- */ 


void sormaniGeneral(const std::string & evolutionFilename, const std::string & barFile, double patternSpeed, bool isSelfConsistent) {
	double timeStep{1};
	int nMax{72}, numbTimeSteps{300};

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", nMax);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); Bar2D bar(coeff, patternSpeed, "Bar2D/barSize.out"); 

	bar.sormaniBar(pd, 1, barFile); // We normalise the bar such that the amplifitude is 1, this mean output is torque/\eplsilon^2  

	std::string kernel = "Kernels/Kalnajs_2_35.out";
	VolterraSolver solver(kernel, nMax, 2, numbTimeSteps, timeStep); 
	solver.activeFraction(.50);
	
 	solver.barRotation(bar, "test.csv", evolutionFilename, isSelfConsistent, false, true);
 	std::cout << "Save evolution to: " << evolutionFilename <<'\n'; 
}

void sormaniDensity(const std::string & densityFilename, const std::string & barFile, int skip, double patternSpeed, bool isSelfConsistent, bool outputDensity) {
	double timeStep{1};
	int nMax{72}, numbTimeSteps{300};

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", nMax);
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1); 
		Bar2D bar(coeff, patternSpeed, "Bar2D/barSize.out"); 

	//bar.sormaniBar(pd, 1, barFile); // We normalise the bar such that the amplifitude is 1, this mean output is torque/\eplsilon^2  

	std::string kernel = "Kernels/Kalnajs_2_35.out";
	VolterraSolver solver(kernel, nMax, 2, numbTimeSteps, timeStep); 
	solver.activeFraction(.50);
	
 	solver.barRotationUnsaving(bar, isSelfConsistent, false, true);
 	if (outputDensity) {solver.density2dEvolution(densityFilename, pd, skip, 5);}
 	else {solver.potential2dEvolution(densityFilename, pd, skip, 5);}
 	
 	std::cout << "Save evolution to: " << densityFilename <<'\n'; 
}


void sormaniPatternSpeed(const std::string & dir) {
	std::vector<double> omegaP{0, 0.06, 0.12, 0.18, 0.24}; 
	auto filenameConsistent = [dir] (double patternSpeed) {return dir + "Kalanajs_evolution_consistent_" +std::to_string(int (patternSpeed *100)) + ".csv";};
	auto filenameTest = [dir] (double patternSpeed)  {return dir + "Kalanajs_evolution_test_" +std::to_string(int (patternSpeed * 100)) + ".csv";};
	auto densityFilenameConsistent = [dir] (double patternSpeed) {return dir + "Kalanajs_density_consistent_" +std::to_string(int (patternSpeed *100)) + ".csv";};
	auto densityFilenameTest = [dir] (double patternSpeed)  {return dir + "Kalanajs_density_test_" +std::to_string(int (patternSpeed * 100)) + ".csv";};

	for (auto op = omegaP.begin(); op != omegaP.end(); ++op){
		sormaniGeneral(filenameConsistent(*op), "Bar2D/Bar_Potentials/Sormani_Medium.out", *op, true);
		sormaniDensity(densityFilenameConsistent(*op), "Bar2D/Bar_Potentials/Sormani_Medium.out", 25, *op, true, true);

		sormaniGeneral(filenameTest(*op), "Bar2D/Bar_Potentials/Sormani_Medium.out", *op, false);
		sormaniDensity(densityFilenameTest(*op), "Bar2D/Bar_Potentials/Sormani_Medium.out", 25, *op, false, true);
	}
}


void sormaniBarShape(const std::string & dir) {
	sormaniGeneral(dir + "Kalnajs_consistent_Small.csv", "Bar2D/Bar_Potentials/Sormani_Small.out", 1/5.5, true);
	sormaniDensity(dir + "Kalnajs_consistent_density_Small.csv", "Bar2D/Bar_Potentials/Sormani_Small.out", 25, 1/5.5, true, true);
	sormaniGeneral(dir + "Kalnajs_test_Small.csv", "Bar2D/Bar_Potentials/Sormani_Small.out", 1/5.5, false);
	sormaniDensity(dir + "Kalnajs_test_density_Small.csv", "Bar2D/Bar_Potentials/Sormani_Small.out", 25, 1/5.5, false, true);

	sormaniGeneral(dir + "Kalnajs_consistent_Medium.csv", "Bar2D/Bar_Potentials/Sormani_Medium.out", 1/5.5, true);
	sormaniDensity(dir + "Kalnajs_consistent_density_Medium.csv", "Bar2D/Bar_Potentials/Sormani_Medium.out", 25, 1/5.5, true, true);
	sormaniGeneral(dir + "Kalnajs_test_Medium.csv", "Bar2D/Bar_Potentials/Sormani_Medium.out", 1/5.5, false);
	sormaniDensity(dir + "Kalnajs_test_density_Medium.csv", "Bar2D/Bar_Potentials/Sormani_Medium.out", 25, 1/5.5, false, true);

	sormaniGeneral(dir + "Kalnajs_consistent_Large.csv", "Bar2D/Bar_Potentials/Sormani_Large.out", 1/5.5, true);
	sormaniDensity(dir + "Kalnajs_consistent_density_Large.csv", "Bar2D/Bar_Potentials/Sormani_Large.out", 25, 1/5.5, true, true);
	sormaniGeneral(dir + "Kalnajs_test_Large.csv", "Bar2D/Bar_Potentials/Sormani_Large.out", 1/5.5, false);
	sormaniDensity(dir + "Kalnajs_test_density_Large.csv", "Bar2D/Bar_Potentials/Sormani_Large.out", 25, 1/5.5, false, true);
} 


void sormaniConstantRatio(const std::string & dir) {
	sormaniGeneral(dir + "Kalnajs_consistent_Small.csv", "Bar2D/Bar_Potentials/Sormani_Small.out", 0.27, true);
	sormaniDensity(dir + "Kalnajs_consistent_density_Small.csv", "Bar2D/Bar_Potentials/Sormani_Small.out", 25, 0.27, true, false);

	sormaniGeneral(dir + "Kalnajs_consistent_Medium.csv", "Bar2D/Bar_Potentials/Sormani_Medium.out", 0.18, true);
	sormaniDensity(dir + "Kalnajs_consistent_density_Medium.csv", "Bar2D/Bar_Potentials/Sormani_Medium.out", 25, 0.18, true, false);

	sormaniGeneral(dir + "Kalnajs_consistent_Large.csv", "Bar2D/Bar_Potentials/Sormani_Large.out", 0.135, true);
	sormaniDensity(dir + "Kalnajs_consistent_density_Large.csv", "Bar2D/Bar_Potentials/Sormani_Large.out", 25, 0.135, true, false);
}