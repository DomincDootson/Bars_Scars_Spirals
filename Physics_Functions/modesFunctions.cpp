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
#include "../DF_Class/ScarredMestel.h"
#include "../DF_Class/Scar.h"
#include "../DF_Class/AngularMomentumScar.h"

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
	EvolutionKernels kernels("Kernels/kalnajs_6_.out", 200);
	ResponseMatrixThread rmt(2);

	std::vector<std::complex<double>> omegas; omegas.emplace_back(0.90,0.22), omegas.emplace_back(0.91,0.22);
	auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<std::complex<double>> vec =  rmt.det(omegas, kernels); 
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "f() took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
    for (auto d : vec) {std::cout << abs(d) <<'\n';}
}

/* Visualising Eigenmodes */ 
/* ---------------------- */ 


void saveResponseMatrix(const std::string & filename, const std::string & kernelFile, const double omega0, const double eta, const double xi) {
	ResponseMatrix rm; std::complex<double> omega(omega0, eta); 
	EvolutionKernels kernels(kernelFile, 800);
	kernels.activeFraction(xi); 

	rm.responseMatrix(omega, kernels);
	rm.saveResponseMatrix(filename);
	std::cout << rm.det() << '\n';
}

void eigenvector2Density(const std::string & modeFile, const std::string & densityFile) {
	ExpansionCoeff mode(modeFile, 1);
	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat", 40);
	mode.write2dDensity2File(densityFile, PD, 1, 6, 301);
}


// Evolution of eigenen Modes

void saveComplexArray(std::ofstream & out, const Eigen::ArrayXXd & array){ 
	for (int i = 0; i < array.rows(); ++i) { 
		for (int j = 0; j < array.rows(); ++j) { 
			if ((i == (array.rows()-1)) && (j == (array.rows()-1))) {out << array(i,j) << '\n';}
			else {out << array(i,j) << ','; }
		}
	} 
} 

void eigenModeEvolution(const std::string & modeFile, const double omega0, const double eta, const std::string & stem) {
	double timeStep{0.1}, time;
	int nTimeStep{100}; 
	ExpansionCoeff posLmode(modeFile, 1); 
	ExpansionCoeff posL(modeFile, nTimeStep), negL(modeFile, nTimeStep);

	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat", 40);

	std::complex<double> omegaTimesi(eta,-omega0);  // check this 

	std::ofstream outPR(stem + "/pos_real.csv"), outPI(stem + "/pos_imag.csv"), outNR(stem + "/neg_real.csv"), outNI(stem + "/neg_imag.csv"); 

	for (int t = 0; t < nTimeStep; ++t){
		std::cout << "At time step: " << t << '\n'; 	 
		time = timeStep * t; 
		posL(t) = exp(omegaTimesi * time) * posLmode(0);
		negL(t) = std::conj(exp(omegaTimesi * time)) * posLmode(0); 

		Eigen::ArrayXXcd array = PD.densityArray(posL(t), 100, 5);
		saveComplexArray(outPR, array.real());
		saveComplexArray(outPI, array.imag());

		array = PD.densityArray(negL(t), 100, 5).conjugate();
		saveComplexArray(outNR, array.real());
		saveComplexArray(outNI, array.imag());
	}
	outPR.close(); outPI.close(); outNR.close(); outNI.close(); 
}

/* Eigemode Search */ 
/* --------------- */ 

void uniformSearchMatrix(const std::string & filename, int nu, double eLower, double oLower) {

	Mestel df(1, 1, 0.377, 1, 1, 11.5, nu, 5);/*
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, 2, 7, 251, 15);

	ResponseMatrixThread rmt(test, df, 2);

	rmt.modeGridSearch(filename, eLower, eLower+0.06, 4, oLower, oLower+0.4, 21);*/
	std::cout << df.diskMass() <<'\n';
}

void unstableTimeEvolution(const std::string & stem) {

	//PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	std::vector<int> chi{4};//,6,8
	auto filename = [&stem] (double X) {return "Plotting/Modes_Data/" + stem + std::to_string(int(ceil(X)))+".csv";};

		for (auto& X : chi)
		{
			VolterraSolver solver("Kernels/Kalnajs_"+std::to_string(X)+"_.out", 40, 2, 400, 0.1);
			//VolterraSolver solver("Kernels/KalnajsN.out", 48, 2, 400, 0.25);
			solver.kernelTesting(filename(X), 0,true);

			//solver.density2dEvolution(199, filename(X), PD, 5, 300);
		}
}


void uniformSearchKernel(const std::string & filename, int nu, double eLower, double oLower) {
	EvolutionKernels kernels("Kernels/kalnajs_4_.out", 800);
	//EvolutionKernels kernels("Kernels/Scarred_Kernels/AM_Scarred_Kernel_R_20_W_25_D_-95_G.out", 800);
	ResponseMatrixThread rmt(2);
	rmt.modeGridSearch(filename, kernels, eLower, eLower+0.15, 50, oLower, oLower+0.60, 60);
}

void uniformSearch(const std::string & filestem, const int taper) {
std::vector<double> xi {};
for (double x = 0.5; x <1.01; x += 0.05) {xi.push_back(x);}

	//EvolutionKernels kernels("Kernels/kalnajs_" + std::to_string(taper) + "_Long.out", 1600);
	EvolutionKernels kernels("Kernels/kalnajs_" + std::to_string(taper) + "_.out", 800);
	ResponseMatrixThread rmt(2);

	for (const auto c : xi) { 
		std::cout << "Finding Modes for xi: " << c << '\n'<<'\n';
		kernels.activeFraction(c);

		auto file = [filestem, taper] (const double x)  {return filestem + "Kernel_Mode_Searching_" + std::to_string(taper) +'_' +std::to_string(int (floor(100 * x))) + ".csv";};
		if (taper == 4)  {rmt.modeGridSearch(file(c), kernels, 0.01, 0.17, 80, 0.65, 0.85, 80);}
		else {rmt.modeGridSearch(file(c), kernels, 0.01, 0.30, 120, 0.70, 0.90, 80);}
		
		
		kernels.resetActiveFraction();
	}
}

/* Scars */
/* ----- */

std::string kernelFilename(const std::string & radius, const std::string & width, const std::string & depth, const bool isLong) {
	std::string kernelFile = "Kernels/Scarred_Kernels/AM_Scarred_Kernel_R_" + radius + "_W_" + width + "_D_" + depth;

	if (isLong) {kernelFile += "_Long";}

	return kernelFile + ".out";
}

std::string outFilename(const std::string & radius, const std::string & width, const std::string & depth) {
	return  "Plotting/Modes_Data/Kernel_Search/Single_Scar/AM_Scarred_Kernel_R_" +
	radius + "_W_" + width + "_D_" + depth + ".csv";
}

void uniformSearchScarredDensity(const std::string & radius, const std::string & width, const std::string & depth, bool isLong) {
	ResponseMatrixThread rmt(2); int nStep{800};

	if (isLong) {nStep = 800;}
	
	//EvolutionKernels kernels(kernelFilename(radius, width, depth, isLong), nStep);
	EvolutionKernels kernels("Kernels/AM_Scarred_Kernel_R_20_W_25_D_-95_Inner.out", nStep);
	kernels.activeFraction(0.5 * (11.6873/12));

	rmt.modeGridSearch(outFilename(radius, width, depth), kernels, 0.001, 0.2, 60, 0.2, 1.5, 30);
	//rmt.modeGridSearch("Plotting/Modes_Data/Kernel_Search/Single_Scar/InnerScarTest.csv", kernels, 0.02, 0.1, 60, 0.3, 1, 60);

}

/* Testing fitting */
/* --------------- */

void fittingMethod(const std::string & outFilename, const std::string & kernelFile, double xi) {
	VolterraSolver solver(kernelFile, 40, 2, 800, 0.5);
	solver.activeFraction(xi);
	solver.kernelTesting(0, true);
	solver.saveResponseCoeff(outFilename);
}

void fittingChi() {
	std::string dir{"Plotting/Modes_Data/Chi_Coeff/"};
	for (double chi = 0.2; chi <= 0.5; chi += 0.01) {
		std::cout << "Chi: " << chi << '\n'; 
		std::string filename = dir + "coeff_4_" + std::to_string(int(chi*100)) + ".csv";
		fittingMethod(filename, "Kernels/Kalnajs_4_.out", chi);

		filename = dir + "coeff_6_" + std::to_string(int(chi*100)) + ".csv";
		fittingMethod(filename, "Kernels/Kalnajs_6_.out", chi);

		filename = dir + "coeff_8_" + std::to_string(int(chi*100)) + ".csv";
		fittingMethod(filename, "Kernels/Kalnajs_8_.out", chi);
	}
}



/* Mode Evolution */
/* -------------- */

void saveModeEvolution(const std::string & filename, const std::string & eigenvalues, double omega, double timeStep) {
	ExpansionCoeff mode(eigenvalues, 100);
	std::complex<double> unitComplex(0,1);
	mode.rotateCoefficents(omega, timeStep);
	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat",40);
	
	mode.write2dDensity2File(filename, PD, 1,10, 201);
 	

}
