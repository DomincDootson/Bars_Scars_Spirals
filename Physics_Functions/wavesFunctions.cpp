#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../Wave/Wave.h"
#include "../Bar2D/Bar2D.h"

#include "../Volterra_Solver/VolterraSolver.h"
#include "../Volterra_Solver/ExpansionCoeff.h"



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


Eigen::VectorXcd deltaFunctionDensity(PotentialDensityPairContainer<KalnajsNBasis> & bf, double centre, double width = 0.5) {

	
	auto ic = [centre, width] (double radius) {return (0.5) * 1/sqrt(2 * M_PI*width*width)*exp(-0.5 *pow((radius-centre)/width , 2));}; 

	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(bf.maxRadialIndex()+1);
	double stepSize{0.01};
	for (int n = 0; n < coef.size(); ++n) {
		for (double r= stepSize; r < 15; r += stepSize) {coef(n) += - 2 * M_PI * bf.potential(r, n) * ic(r) * r;}
	} 

	return stepSize*coef; 
	//coef(0) = 1;
	// return coef;
}

Eigen::VectorXcd deltaFunctionPotential(PotentialDensityPairContainer<KalnajsNBasis> & bf,  double radius) {
	int nGrid{200}; double rMax{10}, spacing{(2*rMax) / ((double) nGrid)};

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(bf.maxRadialIndex()+1);
	Eigen::ArrayXXcd den = bf.potentialArray(coeff, nGrid, rMax);
	den((int) nGrid/2, (int) (nGrid/2 + radius/spacing)) = 1; 
	
	coeff = bf.potentialFitting(den, rMax);

	return coeff;
}

void checkingDeltaFunctionFitting() {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_10_2.dat");	

	auto densityCoef{deltaFunctionDensity(pd, 2.1)}, potentialCoef{deltaFunctionPotential(pd, 2.1)};
	auto arrayD = pd.densityArrayReal(densityCoef, 200, 10);
	auto arrayP = pd.potentialArrayReal(potentialCoef, 200, 10);

	std::ofstream outD("Plotting/Waves_Data/No_Scars/deltaFunctionDensity.csv");
	std::ofstream outP("Plotting/Waves_Data/No_Scars/deltaFunctionPotential.csv");

	for (int i = 0; i < 200; ++i) {
		for (int j = 0; j < 200; ++j) {
			if ((i!= 199) && (j!=199)) {
				outD << arrayD(i,j) <<','; 
				outP << arrayP(i,j) <<','; 
			}
			else {
				outD << arrayD(i,j);
				outP << arrayP(i,j);
			}
		}
	}
	outD.close();
	outP.close();

	for (int i =0; i < densityCoef.size(); ++i) {
		std::cout << i << " " << densityCoef(i) << " " << potentialCoef(i) <<'\n'; 
 	}
}





/* Wave Evolution Functions */
/* ------------------------ */ 

void waveTest(const std::string & filename) {
	ExpansionCoeff coeff(100, 48);
	PotentialDensityPairContainer<KalnajsNBasis> bf("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_10_2.dat");
	double speed{0.07}, R0{8};

	for (int i = 0; i < 100; ++i) {std::cout << i << '\n'; coeff(i) = deltaFunctionDensity(bf, -i*speed +8) + deltaFunctionDensity(bf, i*speed +1);}
	coeff.write2dDensity2File(filename, bf, 1, 10, 201);
}

void waveTestSpinning(const std::string & filename) {
	ExpansionCoeff coeff(100, 48);
	Eigen::VectorXcd bar0_c{Eigen::VectorXcd::Zero(48+1)}, bar1_c{Eigen::VectorXcd::Zero(48+1)}; 
	bar0_c(0) = 1; bar1_c(0) = 1;
	PotentialDensityPairContainer<KalnajsNBasis> bf("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_10_2.dat");

	Bar2D bar0(bar0_c, 2*M_PI), bar1(bar1_c, -4*M_PI); 

	for (int i =0; i < 100; ++i) {
		coeff(i) = bar0.barCoeff() + bar1.barCoeff();
		bar0.drift(0.02,0); bar1.drift(0.02,0);
	}
	coeff.write2dDensity2File(filename, bf, 1, 10, 201);
}


double patternSpeedFromCR(const double CR, const double l = 2) {
	double vc{1};
	double omega0{(l*vc)/CR};
	return omega0/l; 

}

void savingInitialWave() {
	double k{5};
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_10_2.dat", 35);

	Wave spiralLeading(pd, k,   5.5, 0.25);
	spiralLeading.density2dEvolution(0, pd, "Plotting/Spiral_Data/spiralLeading.csv", 10); 

	Wave spiralTrailing(pd, -k, 5.5, 0.25);
	spiralTrailing.density2dEvolution(0, pd, "Plotting/Spiral_Data/spiralTrailing.csv", 10); 

	std::ofstream out("Plotting/Spiralden1D.csv");

	for (double radius = 0.01; radius < 10; radius += 0.05) {out << radius <<',' << spiralTrailing.density(radius, 0) << ',' << spiralLeading.density(radius, 0) <<'\n';}
	out.close(); 
}

void turnOffBarFile(const std::string & filename, double turnOffTime) {
	int nStep{100}; double stepSize{0.5};

	std::ofstream out(filename);
	out << nStep <<'\n';

	for (double time =0; time < nStep*stepSize; time += stepSize) {
		if ( time <= turnOffTime) {out << time << " " << 1 <<'\n';}
		else {out << time << " " << 0 <<'\n';}
	}
	out.close();
}

void turnOnBarFile(const std::string & filename, double turnOnTime) {
	int nStep{100}; double stepSize{0.5};

	std::ofstream out(filename);
	out << nStep <<'\n';

	for (double time =0; time < nStep*stepSize; time += stepSize) {
		if ( time < turnOnTime) {out << time << " " << pow(sin(time *M_PI/(2 * turnOnTime)),2) <<'\n';}
		else {out << time << " " << 1 <<'\n';}
	}
	out.close();
}

void waveEvolutionTest(const std::string & densityfile, double CRposition, double perturberRadius, const std::string & growthFile) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat", 48);
	
	//turnOnBarFile(growthFile, 40);
	turnOffBarFile(growthFile, 150);
	Bar2D bar(deltaFunctionDensity(pd, 2), 0.2);  
	//bar.sormaniBar(pd);

	if (growthFile != "None") {bar.readInSize(growthFile);}

	VolterraSolver solver("Kernels/Waves/Kalnajs_2_WM.out", 48, 2, 100, 1.0);
	solver.activeFraction(0.5);
	//solver.barRotation(bar, true); 

	solver.barRotation(bar, "Plotting/Coeff_test.csv", "Plotting/Evolution_test.csv", true, false);
	solver.density2dEvolution(densityfile, pd, 2, 10);

}




