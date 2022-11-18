#include <iostream>
#include <typeinfo>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <chrono>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"
//#include "../Potential_Density_Pair_Classes/SpiralBasis.h"


#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../Potential_Density_Pair_Classes/TheoreticalPDContainer.h"

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/Mestel.h"

#include "../Volterra_Solver/VolterraSolver.h"
#include "../Volterra_Solver/ExpansionCoeff.h"

#include "../Bar2D/Bar2D.h"

#include "../DF_Function/DFfunction.h"

#include "../Wave/Wave.h"


#include "generalFunctions.h"


// Basis Function Generation //
// ------------------------- // 

// Put in a function to generate the names for the files that the BF are kept it

void generatingKalnajsBF(int m2)
{
	Mestel DF;
	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");

	//ActionAngleBasisContainer test("KalnajsN_Inner", 48, m2, 4, 251, 10); 
	ActionAngleBasisContainer test("KalnajsN", 48, m2, 4, 251, 15); 
	test.scriptW(PD, DF, "KalnajsN"); // Use file function name here
}

void getSpiralParam() {
	std::vector<double> params{24, .5, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params, 24, 2);
}

void generatingSpiralBF(int m2)
{
	Mestel DF(1,1, .33); int nMax{80};
	
	std::vector<double> params{static_cast<double>(nMax), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params, nMax,m2);

	ActionAngleBasisContainer test("GaussianLog", nMax, m2, 7, 251, 20); 
	test.scriptW(PD, DF, "GaussianLog");
}

void savingKalnajsFunctions(const std::string & filename) {
	std::ofstream out(filename);
	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	for (double radius = 0.01; radius < 15; radius += 0.05) {
		out << radius << ',';
	}
	out << 15  << '\n';
	for (int n = 0; n <48; ++n) {
		for (double radius = 0.01; radius < 15; radius += 0.05) {
			out << PD.density(radius, n) << ',';
		}
		out << 0 << '\n';
	}
	out.close();
}

// Kernel Generation //
// ----------------- // 

void generatingKalnajsKernels(const std::string & filename, int m2, int nMax, double rInner)
{
	ActionAngleBasisContainer test("KalnajsN_Small", "KalnajsN", 48, 2, 7, 100, 10);
	Mestel DF(1, 1, 0.377, 1, 1, 11.5, 4, 5);

	VolterraSolver solver2(48, 2, 400, 0.1);

	//std::string kernel2 = "Kernels/kalnajsComparison_" +std::to_string((int) nMax) + "_" + std::to_string((int) rInner) + ".out";
	
	solver2.generateKernel(filename, DF, test);
}

void generatingKalnajsKernelsAxisymmetric(const std::string & filename) {
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, 0, 4, 251, 15);
	Mestel DF(1, 1, 0.238, 1, 1, 11.5, 4, 5);

	VolterraSolver solver2(48, 0, 100, 0.35543);
	solver2.generateKernel(filename, DF, test); 
}

void generatingGaussianKernels(int m2)
{
	ActionAngleBasisContainer test("GaussianLog", "GaussianLog", 40, m2, 7, 251, 20);
	Mestel DF;

	VolterraSolver solver1(48, m2, 200, 0.25);
	std::string kernel1 = "GaussianLogLong.out";
	solver1.generateKernel(kernel1, DF, test);

	/*VolterraSolver solver(40, m2, 200, 0.1);
	std::string kernel = "GaussianLogShort.out";
	solver.generateKernel(kernel, DF, test);*/ 

}

void generatingSpiralBFDiffTemp(int m2) {
	ActionAngleBasisContainer test("GaussianLog", "GaussianLog", 50, m2, 7, 251, 20);
	Mestel DFMed(1,1, .35);

	VolterraSolver solver(50, m2, 200, 0.25);
	std::string kernel = "GaussianLog_Med.out";
	solver.generateKernel(kernel, DFMed, test);	

	Mestel DFcold(1,1, .25);
	VolterraSolver solver1(50, m2, 200, 0.25);
	std::string kernel1 = "GaussianLog_Cold.out";
	solver1.generateKernel(kernel1, DFcold, test);

	Mestel DFwarm(1,1, .45);
	VolterraSolver solver2(50, m2, 200, 0.25);
	std::string kernel2 = "GaussianLog_Warm.out";
	solver2.generateKernel(kernel2, DFwarm, test);
}


void testingBarTorque() {
	std::ofstream out("perturbation.csv");
	out << 10 <<'\n';
	for (int time = 0; time < 2000; ++time){
		for (int i = 0; i < 11; ++i) {
			if (i ==0) {out << 0.01*sin(M_PI * (time /400.0)) << " " << 0 << " ";} 
			else if (i==10) {out << 0 << " " << 0 << '\n';}
			else {out << 0 << " " << 0 << " ";}
		}
	}
	out.close();
	//generatingKalnajsKernels(2);

	bool selfConsistent{true};

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(10+1); coeff(0) = 0.01;
	

	VolterraSolver solver("test200.out", 10, 2, 2000, 0.1);
	solver.activeFraction(.50);
	solver.coefficentEvolution("Plotting/Coeff.csv", "perturbation.csv", true);
	//solver.barRotation(bar, "Plotting/Coeff.csv", "Plotting/evolution.csv", selfConsistent, false, true);
} 

void saveVector(const std::string & filename, const std::vector<double> & vec) {
	std::ofstream out(filename);
	for (auto phi: vec) {
		out << phi << '\n';
	}
	std::cout << "Saving to: " << filename << '\n';
	out.close();
}

void saveKalnajs() {
	std::ofstream outD("Plotting/KalnajsTesting/kalnajsDensity.csv");
	std::ofstream outP("Plotting/KalnajsTesting/kalnajsPotential.csv");

	std::vector<double> params{4, 1};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2);	

	for (double r = 0; r <=1; r+=0.01) {
		outP << PD.potential(r, 0) << ',' << PD.potential(r, 1) << ',' << PD.potential(r, 2) << '\n';
		outD << PD.density(r, 0) << ',' << PD.density(r, 1) << ','  << PD.density(r, 2) << '\n';
	}
	outP.close();
	outD.close(); 
}

void kalnajTest() {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	

	Mestel DF; 

	ActionAngleBasisContainer test("KalnajsN", 48, 2, 7, 251, 20); 
	test.scriptW(pd, DF, "KalnajsN");
	//ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 10, 2, 7, 101, 20);


	VolterraSolver solver1(48, 2, 200, 0.25);
	std::string kernel1 = "Kernels/kalnajsComparison.out";
	solver1.generateKernel(kernel1, DF, test);
}


void energyTapping(int nMax, int rInner) {
	std::string kernel2 = "Kernels/kalnajsComparison_" +std::to_string((int) nMax) + "_" + std::to_string((int) rInner) + ".out"; 
	VolterraSolver solver2(kernel2, nMax, 2, 200, 0.25);
	solver2.activeFraction(.5);
	
	std::string filename = "Plotting/General_Data/kalnajsComparison_" +std::to_string((int) nMax) + "_" + std::to_string((int) rInner) + ".out"; 
	solver2.kernelTesting(filename, 0, true);
}



/* Misc */
/* ---- */

void savingDensity(const std::string & filename) {
	Mestel DF(1, 1, sqrt(1/12.4), 1, 1, 11.5, 4, 5);
	std::cout << DF.diskMass() << '\n';
	std::ofstream out(filename);
	int nStep{100}; double stepSize{0.05};
	out << nStep << ',' << stepSize <<'\n';
	out << 0 << ',' << 0 <<'\n';

	for (double radius = stepSize; radius < nStep * stepSize; radius += stepSize) {out << radius << ',' <<  DF.density(radius) << '\n';}

	out.close(); 
}