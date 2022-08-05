#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../Volterra_Solver/VolterraSolver.h"


void generatingKalnajsKernels(int m2, int nMax)
{
	ActionAngleBasisContainer test("KalnajsN/", "KalnajsN", nMax, m2, 7, 251, 15);
	Mestel DF;

	VolterraSolver solver2(nMax, m2, 500, 0.2);

	std::string kernel2 = "Kernels/KalnajsN_Small_"+ std::to_string((int) m2) + ".out";
	solver2.generateKernel(kernel2, DF, test);
}

std::string coeffFileName(int nMode,  const std::string & prefix, const std::string & dir = "Plotting/Waves_Data") {return dir + '/' + prefix + "_" +std::to_string(nMode) + ".csv";}

void saveArray(const Eigen::ArrayXXd grid) {
	std::ofstream out("Plotting/inhomo.csv");
	for (int i = 0; i < grid.rows(); ++i) {
		for (int j = 0; j < grid.rows()-1; ++j) { 
			out << grid(i,j) << ',';
		}
		out << grid(i, grid.cols()-1) << '\n';
	}
	out.close();
}


Eigen::VectorXcd deltaFunctionDensity(PotentialDensityPairContainer<KalnajsNBasis> & bf,  double radius, int nMax) {
	int nGrid{200}; double rMax{15}, spacing{(2*rMax) / ((double) nGrid)};

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax);
	Eigen::ArrayXXcd den = bf.densityArray(coeff, nGrid, rMax);
	den((int) nGrid/2, (int) (nGrid/2 + radius/spacing)) = 1; 
	
	coeff = bf.densityFitting(den, rMax);

	return coeff;
}

Eigen::VectorXcd deltaFunctionPotential(PotentialDensityPairContainer<KalnajsNBasis> & bf,  double radius, int nMax) {
	int nGrid{200}; double rMax{15}, spacing{(2*rMax) / ((double) nGrid)};

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax);
	Eigen::ArrayXXcd den = bf.potentialArray(coeff, nGrid, rMax);
	den((int) nGrid/2, (int) (nGrid/2 + radius/spacing)) = 1; 
	
	coeff = bf.potentialFitting(den, rMax);

	return coeff;
}

/* Mode Evolution Functions */
/* ------------------------ */ 

void selfConsistentWaves(int nMode) {
	std::string outFilename{coeffFileName(nMode, "Self_Consistent_0_Small")};

	VolterraSolver solver("Kernels/kalnajsN_Small_0.out", 38, 0, 500, 0.50);
	solver.activeFraction(.5);

	Eigen::VectorXcd ic = Eigen::VectorXcd::Zero(solver.maxRadialIndex()); 
	ic[nMode] = 1;
	solver.setInitialCondition(ic);

	solver.solveVolterraEquation(true);
	solver.saveResponseCoeff(outFilename);
}

void perturbationWaves(int nMode) {
	std::string outFilename{coeffFileName(nMode, "Perturbation_0_Small")};

	VolterraSolver solver("Kernels/KalnajsN_Small_0.out", 38, 0, 500, 0.25);
	solver.activeFraction(.5);

	Eigen::VectorXcd ic = Eigen::VectorXcd::Zero(solver.maxRadialIndex()); 
	ic[nMode] = 1;

	solver.kernelTesting(outFilename, nMode, true);	
}




/* Density Evolution Functions l = 0*/ 
/* -------------------------------- */ 


void selfConsistentDensity(double radius) {
	//generatingKalnajsKernels(0);
	PotentialDensityPairContainer<KalnajsNBasis>  bf("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat", 48);

	std::string outFilename{coeffFileName((int) radius, "RadiusPing_Density_0")};

	VolterraSolver solver("Kernels/kalnajsN_0.out", 48, 0, 500, 0.20);
	solver.activeFraction(.5);
	solver.setInitialCondition(deltaFunctionDensity(bf, radius, solver.maxRadialIndex()));

	solver.solveVolterraEquation(true);
	solver.density1dEvolution(outFilename, bf); 
}

void pullingDensity(double radius) {
	
	PotentialDensityPairContainer<KalnajsNBasis>  bf("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat", 48);
	VolterraSolver solver("Kernels/kalnajsN_0.out", 48, 0, 500, 0.20);
	solver.activeFraction(.5);

	std::string outFilename{coeffFileName((int) radius, "RadiusPull_Density_0")};

	int nPull{20}; std::vector<Eigen::VectorXcd> deltas; deltas.reserve(nPull); 
	Eigen::VectorXcd delta{deltaFunctionDensity(bf, radius, solver.maxRadialIndex())};
	for (int time = 0; time < nPull; ++time) {deltas.emplace_back(delta * (time /((double) nPull)));}


	
	solver.setInitialCondition(deltas);

	solver.solveVolterraEquation(true, nPull);
	solver.density1dEvolution(outFilename, bf); 
}

/* Potential Evolution Functions l = 0 */
/* ----------------------------------- */ 


void selfConsistentPotential(double radius) {
	//generatingKalnajsKernels(0);
	PotentialDensityPairContainer<KalnajsNBasis>  bf("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat", 48);

	std::string outFilename{coeffFileName((int) radius, "RadiusPing_Potential_0")};

	VolterraSolver solver("Kernels/kalnajsN_0.out", 48, 0, 500, 0.20);
	solver.activeFraction(.5);
	solver.setInitialCondition(deltaFunctionDensity(bf, radius, solver.maxRadialIndex()));

	solver.solveVolterraEquation(true);
	solver.density1dEvolution(outFilename, bf); 
}

void pullingPotential(double radius) {
	
	PotentialDensityPairContainer<KalnajsNBasis>  bf("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat", 48);
	VolterraSolver solver("Kernels/kalnajsN_0.out", 48, 0, 500, 0.20);
	solver.activeFraction(.5);

	std::string outFilename{coeffFileName((int) radius, "RadiusPull_Potential_0")};

	int nPull{20}; std::vector<Eigen::VectorXcd> deltas; deltas.reserve(nPull); 
	Eigen::VectorXcd delta{deltaFunctionDensity(bf, radius, solver.maxRadialIndex())};
	for (int time = 0; time < nPull; ++time) {deltas.emplace_back(delta * (time /((double) nPull)));}


	
	solver.setInitialCondition(deltas);

	solver.solveVolterraEquation(true, nPull);
	solver.density1dEvolution(outFilename, bf); 
}

/* Density Evolution Functions l = 2*/ 
/* -------------------------------- */ 

void selfConsistentQuadDensity(double radius) {
	//generatingKalnajsKernels(0);
	PotentialDensityPairContainer<KalnajsNBasis>  bf("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat", 48);

	std::string outFilename{coeffFileName((int) radius, "RadiusPing_Density_2")};

	VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 500, 0.50);
	solver.activeFraction(.5);
	solver.setInitialCondition(deltaFunctionDensity(bf, radius, solver.maxRadialIndex()));

	solver.solveVolterraEquation(true);
	solver.density2dEvolution(outFilename, bf); 
}
 
