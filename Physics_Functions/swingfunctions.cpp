#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <numeric>

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


/* Axisymmetric Evolution */
/* ---------------------- */

Eigen::VectorXcd guassianRing(double centre, double width) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");

	auto ic = [centre, width] (double radius) {return   radius/sqrt(2 * M_PI*width*width)*exp(-0.5 *pow((radius-centre)/width , 2));}; 

	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(49);
	double stepSize{0.01};
	for (int n = 0; n < coef.size(); ++n) {
		for (double r= stepSize; r < 15; r += stepSize) {coef(n) += - 2 * M_PI * pd.potential(r, n) * ic(r) * r;}
	}

	return stepSize*coef; 
}

void ringEvolution(const std::string & filename, bool isSelfConsistent) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");
	VolterraSolver solver("Kernels/Kalnajs_0.out", 48, 0, 100, 0.35543); 
	solver.activeFraction(0.5);

	solver.setInitalPerturbation(guassianRing(8, 1)); 	
	if (isSelfConsistent) {solver.deltaPerturbationConsistent();}
	else {solver.deltaPerturbationTest();}

	std::vector<double> radii(100);
	std::iota(radii.begin(), radii.end(), 0); 
	std::for_each(radii.begin(), radii.end(), [] (double & n){n =(0.02*n)+4;});

	solver.density1dEvolution(filename, pd); 
}

/* Non-Axisymmetric Evolution */
/* -------------------------- */ 

void generateKernel() {
	int l{8};
	PotentialDensityPairContainer<KalnajsNBasis> PD("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_" +std::to_string(l) +".dat");

	Mestel DF(1, 1, 0.238, 1, 1, 11.5, 4, 5);
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, l, 4, 251, 15); 
	//test.scriptW(PD, DF, "KalnajsN"); // Use file function name here

	VolterraSolver solver2(48, l, 100, 0.35543);
	solver2.generateKernel("Kernels/Kalnajs_kernel_" +std::to_string(l) +".out", DF, test); 
}

Eigen::VectorXcd perturbationCoef(double centre, double width) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_8.dat", 40);

	auto ic = [centre, width] (double radius) {return (0.5) * 1/sqrt(2 * M_PI*width*width)*exp(-0.5 *pow((radius-centre)/width , 2));}; 

	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(41);
	double stepSize{0.01};
	for (int n = 0; n < coef.size(); ++n) {
		for (double r= stepSize; r < 15; r += stepSize) {coef(n) += - 2 * M_PI * pd.potential(r, n) * ic(r) * r;}
	}

	/*std::ofstream out("Plotting/test.csv");
	for (double r = 0.01; r < 15; r+=0.05) {
		out << r << ',' << ic(r) <<',';
		std::complex<double> sum{0};
		for (int i = 0; i < pd.maxRadialIndex()+1; ++i) {sum += stepSize*coef(i)*pd.density(r, i);}
		out << sum.real() <<'\n';
	}
	out.close();*/

	return stepSize*coef; 
}


void discComparison(const std::string & filename, bool isSelfConsistent) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_8.dat");

	VolterraSolver solver("Kernels/Kalnajs_kernel_8.out", 48, 8, 100, 0.35543);
	solver.activeFraction(0.5); 
	solver.setInitalPerturbation(perturbationCoef(8, 0.5)); 
	
	// solver.deltaPerturbationTest(); 
	// solver.density1dCorotating("Plotting/test.csv", pd, 0.125, 4, 12, 1); 
	
	solver.deltaPerturbationConsistent(); 
	solver.density1dCorotating("Plotting/cons.csv", pd, 0.125, 4, 12, 1); 
}

/* Generating Multiple Kernels */ 
/* --------------------------- */ 


void generateIndividualSwingKernel(int l, const double timeStep = 0.3, const int nStep = 500, const int nMax = 40, const double temp = 0.24) {
	const std::string filename{"Kernels/Swing_Kernels/SW_Ker_Cold_" + std::to_string(l) + ".csv"};
	Mestel DF(1, 1, temp, 1, 1, 11.5, 4, 5);

	
	VolterraSolver solver(nMax, l, nStep, timeStep);	
	ActionAngleBasisContainer pd("KalnajsN/Swing_Amplification", "KalnajsN", nMax, l, 4, 251, 15);

	solver.generateKernel(filename, DF, pd); 

}

void generateIndividualSwingKernel_name(const std::string filename, int l, const double timeStep = 0.3, const int nStep = 500, const int nMax = 40, const double temp = 0.24) {
	Mestel DF(1, 1, temp, 1, 1, 11.5, 4, 5);
	VolterraSolver solver(nMax, l, nStep, timeStep);	
	ActionAngleBasisContainer pd("KalnajsN/Swing_Amplification", "KalnajsN", nMax, l, 4, 251, 15);
	solver.generateKernel(filename, DF, pd); 
}


void generateSwingKernels(int lMax) {
	for (int i = 0; i < 15; ++i) {
	std::cout << "Calculating kernels for l: " << i <<'\n';
	generateIndividualSwingKernel(i);}
}



/* Fixed Radius Amplification */
/* -------------------------- */


void writeVector2file(const std::vector<double> & values, const std::string & filename, int startHarmonic = 0, int endHarmonic = 15) {
	std::ofstream out(filename);
	for (int l = 0; l < endHarmonic; ++l) {out << l <<',' << values[l] << '\n';}
	out.close();
} 

void writeVector2file(const std::vector<double> & x_vec, const std::vector<double> & y_vec, const std::string & filename) {
	std::ofstream out(filename);
	for (int i =0; i < x_vec.size(); ++i) {out << x_vec[i] <<',' << y_vec[i] << '\n';}
	out.close();
} 

void saveRotationFile(const std::string filename, double timeAlive) {
	std::ofstream out(filename); int nStep{500}; double timeStep(400/((double) nStep));
	out << nStep <<'\n';
	for (double time = 0; time<timeAlive; time += timeStep) {out << time << " " << 1 <<'\n';}
	for (double time = timeAlive; time < 400; time += timeStep) {out << time << " " << 0 <<'\n';}
	//for (double time =0; time<1; time += timeStep) {out << time << " " << 1 <<'\n';}
	//for (double time =1; time<400; time += timeStep) {out << time << " " << 0 <<'\n';}
	out.close();
}

double amplificationFixedRadiusGivenl(const double rad, const int l, const bool isSelfConsistent, const double isClockwise) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_" +std::to_string(l) +".dat", 40);
	Eigen::VectorXcd coef(40+1);

	for (int n = 0; n <= pd.maxRadialIndex(); ++n) {coef(n) = - 2 * M_PI * pd.potential(rad, n);}
	Bar2D bar(coef, isClockwise/rad, "Bar2D/Swing_Files/Radius_"+std::to_string((int) rad)+".out", l); 

	const std::string kernelFile{"Kernels/Swing_Kernels/SW_Ker_Cold_" + std::to_string(l) + ".csv"};
	VolterraSolver solver(kernelFile, 40, l, 500, 0.3);
	
	solver.activeFraction(0.5); 
	solver.barRotationUnsaving(bar, isSelfConsistent, false, true);

	return solver.maxDensity(pd); 
}

void amplificationFixedRadius(double rad, int startHarmonic, int endHarmonic) {
	saveRotationFile("Bar2D/Swing_Files/Radius_"+std::to_string((int) rad)+".out", 0.5*M_PI * rad);  // Keep alive for quarter a turn
	
	std::vector<double> cc, ca, tc, ta;
	for (int l = 0; l < 15; ++l) {
		cc.emplace_back(amplificationFixedRadiusGivenl(rad, l, true, 1));
		ca.emplace_back(amplificationFixedRadiusGivenl(rad, l, true, -1));
		tc.emplace_back(amplificationFixedRadiusGivenl(rad, l, false, 1));
		ta.emplace_back(amplificationFixedRadiusGivenl(rad, l, false, -1));
	}

	writeVector2file(cc, "Plotting/Swing_Data/Amplification_Fixed_5/ConsistentClockwise_Cold.csv", startHarmonic, endHarmonic); 
	writeVector2file(ca, "Plotting/Swing_Data/Amplification_Fixed_5/ConsistentAntiClockwise_Cold.csv", startHarmonic, endHarmonic); 
	writeVector2file(tc, "Plotting/Swing_Data/Amplification_Fixed_5/TestClockwise_Cold.csv", startHarmonic, endHarmonic); 
	writeVector2file(ta, "Plotting/Swing_Data/Amplification_Fixed_5/TestAntiClockwise_Cold.csv", startHarmonic, endHarmonic); 
}

/* Evolution */
/* --------- */ 

void evolutionFixedRadiusGivenl(const std::string & filestem, const double rad, const int l, const bool isSelfConsistent, const double isClockwise, const int skip = 50) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_" +std::to_string(l) +".dat", 40);
	Eigen::VectorXcd coef(40+1);

	for (int n = 0; n <= pd.maxRadialIndex(); ++n) {coef(n) = - 2 * M_PI * pd.potential(rad, n);}
	Bar2D bar(coef, isClockwise/rad, "Bar2D/Swing_Files/Radius_"+std::to_string((int) rad)+".out", l); 

	const std::string kernelFile{"Kernels/Swing_Kernels/SW_Ker_" + std::to_string(l) + ".csv"};
	VolterraSolver solver(kernelFile, 40, l, 500, 0.3);
	
	solver.activeFraction(0.5); 
	solver.barRotationUnsaving(bar, isSelfConsistent, false, true);

	std::string filename = "Plotting/Swing_Data/Amplification_Fixed_5/" + filestem + std::to_string(l) + ".csv"; 
	solver.density2dEvolution(filename , pd, skip, 8);  

}

void densityEvolutionFixedRadius(double rad, int startHarmonic, int endHarmonic) {
	saveRotationFile("Bar2D/Swing_Files/Radius_"+std::to_string((int) rad)+".out", 0.5*M_PI * rad);  // Keep alive for half a turn
	
	
	for (int l = 0; l < 15; ++l) {
		evolutionFixedRadiusGivenl("CC_Evolution_", rad, l, true, 1);
		evolutionFixedRadiusGivenl("CA_Evolution_", rad, l, true, -1);
		evolutionFixedRadiusGivenl("TC_Evolution_", rad, l, false, 1);
		evolutionFixedRadiusGivenl("TA_Evolution_", rad, l, false, -1);
	}
}

/* Cold Disc limit */ 
/* --------------- */ 

void generateSwingKernelsTemp(int l, double Q) {
	std::vector<double> chi_d {0.05, 0.10, 0.15,0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
	std::vector<std::string> chi_s {"05", "10", "15", "20", "25", "30", "35", "40", "45", "50"};

	std::string Q_s{std::to_string(int(Q*10))};

	auto sigma = [Q] (double xi) {return Q * xi * (3.36/(2*M_PI*sqrt(2)));};

	for (int i =0; i < chi_d.size(); ++i) {
		std::string filename = "Kernels/Swing_Kernels/Diff_Temp/SW_Ker_" + std::to_string(l) +"_" + Q_s +"_"+chi_s[i]+ ".out";
		std::cout << filename <<'\n';
		generateIndividualSwingKernel_name(filename, l, 1, 150, 40, sigma(chi_d[i]));
	}
}

double amplificationGivenChi(const double rad, const int l, const bool isSelfConsistent, const double activeFraction, const std::string & chi_s) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_" +std::to_string(l) +".dat", 40);
	Eigen::VectorXcd coef(40+1);

	for (int n = 0; n <= pd.maxRadialIndex(); ++n) {coef(n) = - 2 * M_PI * pd.potential(rad, n);}
	Bar2D bar(coef, 1/rad, "Bar2D/Swing_Files/Radius_"+std::to_string((int) rad)+".out", l); 

	  

	std::string kernelFile = "Kernels/Swing_Kernels/Diff_Temp/SW_Ker_" + std::to_string(l) +"_13_"+chi_s+ ".out";
	VolterraSolver solver(kernelFile, 40, l, 150, 1);
	
	solver.activeFraction(activeFraction); 
	solver.barRotationUnsaving(bar, isSelfConsistent, false, true);

	return solver.maxDensity(pd); 
}

void densityEvolutionChi(double rad, int harmonic) {
	saveRotationFile("Bar2D/Swing_Files/Radius_"+std::to_string((int) rad)+".out", 0.5*M_PI * rad);  // Keep alive for half a turn
	
	std::vector<double> chi_d {0.05, 0.10, 0.15,0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
	std::vector<std::string> chi_s {"05", "10", "15", "20", "25", "30", "35", "40", "45", "50"};
	
	std::vector<double> c, t;
	for (int i = 0; i < chi_s.size(); ++i) {
		std::cout << "Active Fraction: " << chi_d[i] << '\n';
		c.emplace_back(amplificationGivenChi(5, harmonic, true, chi_d[i], chi_s[i]));
		t.emplace_back(amplificationGivenChi(5, harmonic, false, chi_d[i], chi_s[i]));
		
		
	}

	
	
	
	
	writeVector2file(chi_d, c, "Disc_Consistent_13_" +std::to_string(harmonic)+".csv");
	writeVector2file(chi_d, t, "Disc_Test_13_" +std::to_string(harmonic)+".csv");


}



