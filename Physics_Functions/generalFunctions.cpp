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

// Basis Function Generation //
// ------------------------- // 

// Put in a function to generate the names for the files that the BF are kept it

void generatingKalnajsBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,m2);

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
	ActionAngleBasisContainer test("Kalnajs/Kalnajs_4_20", "Kalnajs", 10, m2, 7, 101, 20);
	Mestel DF(1, 1, 0.35);

	VolterraSolver solver2(10, m2, 2000, 0.1);

	std::string kernel2 = "test200.out";//"Kernels/Kalnajs" +std::to_string(m2) +".out";
	solver2.generateKernel(kernel2, DF, test);
}

void generatingGaussianKernels(int m2)
{
	ActionAngleBasisContainer test("GaussianLog", "GaussianLog", 40, m2, 7, 251, 20);
	Mestel DF(1,1, .45);

	VolterraSolver solver1(40, m2, 2000, 0.1);
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

void testingFitting() {
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> kalnajs(params, 10, 2);
	std::complex<double> i(0,1);

 	Eigen::VectorXcd coefG = Eigen::VectorXcd::Zero(10+1);
 
 	coefG(0) = 1;


 	Eigen::ArrayXXcd pot = kalnajs.potentialArray(coefG, 1500, 20);
 	
	
 	int nMax{24};
	std::vector<double> params1{static_cast<double>(nMax), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> pd(params1, nMax, 2);

 	/*Eigen::VectorXcd coefG = Eigen::VectorXcd::Zero(nMax+1);
 	coefG(10) = 1;
 	coefG(15) = 2; */ 

 	//Eigen::ArrayXXcd pot = pd.potentialArray(coefG, 1000, 20);
 	Eigen::VectorXcd fitcoef = pd.potentialFitting(pot, 20);

 	std::vector<double> radiiVector;
 	for (int i = 1; i < 2000; ++i) {radiiVector.push_back(i*0.01);}
 	saveVector("Plotting/fitGaussian.csv", pd.oneDpotential(radiiVector, fitcoef));
 	saveVector("Plotting/trueKalnajs.csv", kalnajs.oneDpotential(radiiVector, coefG));
	
 	std::cout << fitcoef << '\n';
 	


}
