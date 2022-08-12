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

void testingKernelMethod() {
	EvolutionKernels kernels("Kernels/kalnajsN.out", 200);
	ResponseMatrixThread rmt(2);

	std::vector<std::complex<double>> omegas; omegas.emplace_back(1,-0.1), omegas.emplace_back(1,-0.2);
	auto t1 = std::chrono::high_resolution_clock::now();
	rmt.det(omegas, kernels); 
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "f() took "
              << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()
              << " milliseconds\n";

}

/* Visualising Eigenmodes */ 
/* ---------------------- */ 

template <class Tbf>
void removingResiduals(const std::string & residualFile, const std::string & timefile, ExpansionCoeff & mode, const Tbf & bf) {
	
	VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 200, 0.5);
	solver.activeFraction(1);
	solver.kernelTesting(0,true);

	ExpansionCoeff residuals(solver.maxRadialIndex()-1, 1); 
	residuals(0) = -(solver.responseCoef(199) - mode(0).dot(solver.responseCoef(199))*mode(0)); 
	std::cout << "Here\n";
	residuals.write2dDensity2File(0, residualFile,  bf, 5, 301);

	//for (int i = 0; i < 40; ++i) {std::cout << residuals(0)(i) << " "<< mode(0)(i) <<'\n';}
	std::cout << "We got here\n";

	mode(0) = solver.responseCoef(199);
	mode.write2dDensity2File(0, timefile,  bf, 5, 301);
}

void generatingEigenMode(const std::string & modeFile, const std::string & densityFile, const std::string & residualFile, const std::string & timefile) {
	ExpansionCoeff mode(modeFile, 1);
	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	mode.write2dDensity2File(densityFile, PD, 1, 5, 301);

	removingResiduals(residualFile, timefile, mode, PD);
}




/* Eigemode Search */ 
/* --------------- */ 

void uniformSeach(const std::string & filename, int nu) {
	/*Mestel df(1, 1, 0.377, 1, 1, 11.5, nu, 5);
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 40, 2, 7, 251, 15);


	//ResponseMatrixThread rmt(test, df, 2);

	ResponseMatrix rm(test, df);

	//std::vector<std::complex<double>> omegas; omegas.emplace_back(0.88, 0.13); omegas.emplace_back(0.88, 0.14); omegas.emplace_back(0.88, 0.12);
	std::vector<std::complex<double>> omegas; omegas.emplace_back(0.88, 0.13);omegas.emplace_back(0.87, 0.13); omegas.emplace_back(0.89, 0.13);
	for (auto o : omegas) {std::cout << abs(rm.det(o)) <<'\n';}

	

	//rmt.modeGridSearch(filename, 0.13, 0.17, 3, 0.75, 0.95, 21);*/

	Mestel df(1, 1, 0.377, 1, 1, 11.5, 4, 5);
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, 2, 4, 251, 15);

	ResponseMatrixThread rmt(test, df, 2);

	rmt.modeGridSearch("test.csv", 0.09, 0.17, 5, 0.6, 0.96, 7);
}

void unstableTimeEvolution(const std::string & stem) {

	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	std::vector<double> chi{0.5, 0.6, 0.7, 0.8, 0.9, 1};
	auto filename = [&stem] (double X) {return "Plotting/Modes_Data/" + stem + std::to_string(int(X*10))+".csv";};

		for (auto& X : chi)
		{
			VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 200, 0.5);
			solver.activeFraction(X);
			solver.kernelTesting(filename(X), 0,true);
			//solver.density2dEvolution(199, filename(X), PD, 5, 300);
		}
}

void unstableSearch(const std::string & filename) {
	EvolutionKernels kernels("Kernels/kalnajsShortN.out", 1200);
	
	ResponseMatrixThread rmt(2);

	rmt.modeGridSearch(filename, kernels, 0.09, 0.17, 5, 0.6, 0.96, 7); 
}
