#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <complex>
#include <chrono>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/Mestel.h"
#include "../Volterra_Solver/ExpansionCoeff.h"
#include "../Volterra_Solver/VolterraSolver.h"

#include "../Response_Matrix/ResponseMatrix.h"
#include "../Response_Matrix/ResponseMatrixThread.h"

void m1Stability() {
	Mestel df(1, 1, 0.377, 1, 1, 11.5, 4, 5);

	std::vector<int> vec {1,2,3,4,5,6,7};
	std::complex<double> omega(0.88,0.13); 
	std::vector<double> error;

	for (auto & m1 : vec) {
		ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, 2, m1, 251, 15);
		ResponseMatrix rm(test, df);
		error.emplace_back( abs(rm.det(omega)));} 

	for (auto & err : error) {std::cout << err <<'\n';}
}


/* Visualising Eigenmodes */ 

void generatingEigenMode(const std::string & modeFile, const std::string & outFile) {
	ExpansionCoeff coeff(modeFile, 1);
	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	coeff.write2dDensity2File(outFile,  PD, 1, 5, 301);
}

void unstableTimeEvolution(const std::string & densityFile, const std::string & coeffFile) {

	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 200, 0.5);
	solver.activeFraction(1);
	solver.kernelTesting(coeffFile, 0,true);

	solver.density2dEvolution(199, densityFile, PD, 5, 300);
}

/* Eigemode Search */ 
/* --------------- */ 

void uniformSeach() {
	Mestel df(1, 1, 0.377, 1, 1, 11.5, 4, 5);
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, 2, 4, 251, 15);

	ResponseMatrixThread rmt(test, df, 2);

	rmt.modeGridSearch("test.csv", 0.08, 0.16, 5, 0.6, 1.2, 11);
}