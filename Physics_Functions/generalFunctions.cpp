#include <iostream>
#include <typeinfo>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <chrono>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "../Potential_Density_Pair_Classes/SpiralBasis.h"


#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../Potential_Density_Pair_Classes/TheoreticalPDContainer.h"

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/Mestel.h"

#include "../Volterra_Solver/VolterraSolver.h"
#include "../Volterra_Solver/ExpansionCoeff.h"

#include "../Bar2D/Bar2D.h"

#include "../DF_Function/DFfunction.h"

// Basis Function Generation //
// ------------------------- // 

// Put in a function to generate the names for the files that the BF are kept it

void generatingKalnajsBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, m2);

	ActionAngleBasisContainer test("Kalnajs", 10, m2, 7, 101, 20); 
	test.scriptW(PD, DF, "Kalnajs/Kalnajs_4_20"); // Use file function name here
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

// Kernel Generation //
// ----------------- // 

void generatingKalnajsKernels(int m2)
{
	ActionAngleBasisContainer test("Kalnajs/Kalnajs_4_15", "Kalnajs", 10, m2, 7, 101, 20);
	Mestel DF;

	VolterraSolver solver2(10, m2, 200, 0.1);

	std::string kernel2 = "Kernels/kalnajsComparison.out";//"Kernels/Kalnajs" +std::to_string(m2) +".out";
	solver2.generateKernel(kernel2, DF, test);
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
	generatingKalnajsKernels(2);

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



void spiralBasisComparions() {
	std::vector<double> params{1, 0.1, 50};
	TheoreticalPDContainer<SpiralBasis> PD(params, 0, 2); 

	std::cout << PD.scriptE() << '\n' << '\n';

	Mestel DF; 
	ActionAngleBasisContainer testOld("OldSpiral", 0, 2, 0, 100, 2);

	for (int n = 5; n < 55; ++n) {
	std::cout << n << " ";
	
	std::cout <<testOld.scriptW(PD, DF, 0, 0, 56, n)/testOld.analyticScriptW(PD, DF, 0, 0,  56, n) -1 << '\n'; }

}


void spiralTestAnalytic() {
	std::vector<double> params{2, 0.05, 50};
	TheoreticalPDContainer<SpiralBasis> PD(params, 10, 2); 
 	PD.analyticKernel("Kernels/SpiralAnalytic.out", 100, 0.25, 0.1);
	VolterraSolver solver1("Kernels/SpiralAnalytic.out", 10, 2, 100, 0.25);
	

	std::ofstream out("Plotting/Spiral_Data/AnalyticKernel.csv");

	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			if (j==10) {out << abs(solver1(1)(i,j)) << '\n';}
			else {out << abs(solver1(1)(i,j)) << ',';}
		}
	}
	out.close();

} 

void spiralTestQuasi() {
	std::vector<double> params{2, 0.05, 50};
	TheoreticalPDContainer<SpiralBasis> PD(params, 10, 2); 
	Mestel DF(1,1,0.10); 


	ActionAngleBasisContainer test("Spiral", 10, 2, 7, 100, 5);
	test.analyticScriptW(PD, DF, "SpiralAABF");

	VolterraSolver solver1(10, 2, 100, 0.25);
	std::string kernel1 = "Kernels/SpiralQuasi.out";
	solver1.generateKernel(kernel1, DF, test);

	std::ofstream out("Plotting/Spiral_Data/QuasiKernel.csv");

	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			if (j==10) {out << abs(solver1(1)(i,j)) << '\n';}
			else {out << abs(solver1(1)(i,j)) << ',';}
		}
	}
	out.close();

}

void spiralTestTrue() {
	std::vector<double> params{0.2, 0.05, 50};
	TheoreticalPDContainer<SpiralBasis> PD(params, 10, 2); 
	Mestel DF(1,1,0.1); 

	ActionAngleBasisContainer testOld("OldSpiral", 10, 2, 7, 100, 5);
	testOld.scriptW(PD, DF, "SpiralAABF");


	VolterraSolver solver1(10, 2, 100, 0.25);
	std::string kernel1 = "Kernels/SpiralTrue.out";
	solver1.generateKernel(kernel1, DF, testOld);
}


void spiralEvolution() {
	std::vector<std::string> kernels = {"Kernels/SpiralAnalytic.out", "Kernels/SpiralQuasi.out","Kernels/SpiralTrue.out"};
	std::vector<std::string> files = {"Plotting/Spiral_Data/SpiralAnalyticTest.csv", "Plotting/Spiral_Data/SpiralQuasiTest.csv","Plotting/Spiral_Data/SpiralTrueTest.csv"};

	for (int i = 0; i < 3; ++i) {
		VolterraSolver solver1(kernels[i], 10, 2, 100, 0.25);
		solver1.kernelTesting(files[i], 3, false);
	}





}