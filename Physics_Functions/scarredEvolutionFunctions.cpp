#include <cmath>
#include <iostream>
#include <fstream>

#include <chrono>

#include <boost/math/special_functions/ellint_1.hpp>

#include "../DF_Class/Mestel.h"
#include "../DF_Class/ScarredMestel.h"
#include "../DF_Class/Scar.h"
#include "../DF_Class/AngularMomentumScar.h"

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../Volterra_Solver/VolterraSolver.h"

void addScars(ScarredMestel<Scar> & df, bool scar1 = true, bool scar2 = true, bool scar3 = true) { 
	if (scar1) {Scar scar1(1.20, 11, 0.29, 0.474, -0.044, -5.219, -44.796); df.addScar(scar1);}
	//if (scar2) {Scar scar2(1.55,  6, 0.10, 0.562,  0.288, -3.371, -15.619); df.addScar(scar2);}
	//if (scar3) {Scar scar3(1.90,  6, 0.12, 0.081,  0.869, -0.486,  -3.242); df.addScar(scar3);}
}

void addScars(ScarredMestel<AngularMomentumScar> & df, bool scar1 = true, bool scar2 = true, bool scar3 = true) { 
	if (scar1) {AngularMomentumScar scar1(2, 0.1, 5); df.addScar(scar1);}
	//if (scar2) {Scar scar2(1.55,  6, 0.10, 0.562,  0.288, -3.371, -15.619); df.addScar(scar2);}
	//if (scar3) {Scar scar3(1.90,  6, 0.12, 0.081,  0.869, -0.486,  -3.242); df.addScar(scar3);}
}


void savingEvolutionKernel(bool toAddScars, double sigma = 0.35, int nTimeStep = 200, double timeStep = 0.1, const std::string & filename = "Kernels/kalnajsScarred.out", const bool generateBF = false) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat"); int l{0};
	//ScarredMestel<Scar> DF(1, 1, sigma); 
	ScarredMestel<AngularMomentumScar> DF(1, 1, sigma); 
	double unscarredDiskMass{DF.diskMass()};
	if (toAddScars) {addScars(DF);} 
	double scarredDiskMass{DF.diskMass()}; std::cout << "Ratio of disk masses: " << unscarredDiskMass/scarredDiskMass <<'\n';
	if (generateBF) {
		ActionAngleBasisContainer test("KalnajsN", 48, l, 7, 251, 20); test.scriptW(pd, DF, "KalnajsN");}
	
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, l, 7, 251, 20);

	VolterraSolver solver1(48, l, nTimeStep, timeStep);
	solver1.generateKernel(filename, DF, test, unscarredDiskMass/scarredDiskMass);
}



void circularCutThrought() {
	ScarredMestel<Scar> df; 

	std::ofstream out("Plotting/Scar_Data/withScar.csv");
	for (double r = 0.1; r<3; r += 0.01) {
		
		out << df.distFunc(0.5 + log(r), r) << '\n';
	}
	out.close(); 

	addScars(df); 

	std::ofstream out1("Plotting/Scar_Data/withoutScar.csv");
	for (double r = 0.1; r<3; r += 0.01) {
		
		out1 << df.distFunc(0.5 + log(r), r) << '\n';
	}
	out1.close(); 

}

void scarredDensity(const double temp, const std::string & filename) {
	ScarredMestel<AngularMomentumScar> df(1,1, temp); 

	/*std::ofstream out("Plotting/Scar_Data/unscarredDensity.csv");
	for (double r = 0.1; r<10; r += 0.01) {std::cout << "Density for point: " << r << '\n'; out << r << ',' << df.density(r) << '\n';}
	out.close(); */ 

	addScars(df); 
	std::ofstream out1(filename);
	for (double r = 0.1; r<10; r += 0.01) {std::cout << "Density for point: " << r << '\n'; out1 << r << ',' << df.density(r) << '\n';}
	out1.close();  
}

void scarredPotential() {
	
	ScarredMestel<Scar> df(1,1,0.35);
 
	std::ofstream out("Plotting/Scar_Data/unscarredPotential.csv");
	for (double r = 0.1; r<15; r += 0.01) {std::cout << "Potential for point: " << r << '\n'; out << r << ',' << df.potentialFromVector(r) << '\n';}
	out.close(); 

	Scar scar1(1.20, 11, 0.29, 0.474, -0.044, -5.219, -44.796); 
	ScarredMestel<Scar> dfS({scar1}, 1,1,0.35);
 
	std::ofstream out1("Plotting/Scar_Data/scarredPotential.csv");
	for (double r = 0.1; r<15; r += 0.01) {std::cout << "Potential for point: " << r << '\n'; out1 << r << ',' << dfS.potentialFromVector(r) << '\n';}
	out1.close(); 

}


void backgroundDF() {
	ScarredMestel<Scar> DF(1, 1, 0.35); 
	addScars(DF);
	int nL{2000}, nE{3000};
	double upperE{0.5+log(8)}, upperL{5}, stepE{upperE/((double) nE)}, stepL{upperL/((double) nL)}, E, L; 
	std::string nan = "nan";
	auto lessThanCircular  = [nan, DF] (double E, double L) {if (E < (0.5 +log(L))) {return nan;} else {return std::to_string(DF.distFunc(E,L));}};

	std::ofstream out("Plotting/Scar_Data/scarredDF.csv"); 
	for (int i = 0; i < nE; ++i) {
		for (int j = 0; j < nL; ++j) {
			E = i * stepE; L = j * stepL; 
			if (j == nL-1) {out << lessThanCircular(E,L) << '\n';}
			else {out <<  lessThanCircular(E,L)<< ',';}
		}
	}

	out.close();
}

void scarredModes() {
	VolterraSolver solver1("Kernels/kalnajsScarred0.out", 48, 0, 200, 0.1);
	solver1.activeFraction(0.5); 
	solver1.kernelTesting("Plotting/Scar_Data/scarredEvolutionModesT.csv", 48, false);
}

// Bouncing off the Scars //
// ---------------------- // 

void inFallingCoefficents(double timeStep, int nTimeStep) {
	std::vector<double> positions, velocity; 
	positions.push_back(10); velocity.push_back(0);
	
	auto accel = [] (double x) {return -1/x;};
	do
	 {
	 	double newPosition = positions.back() + timeStep * velocity.back() + 0.5 * timeStep * timeStep * accel(positions.back());
	 	double newVelocity = velocity.back() + 0.5 * timeStep * (accel(newPosition) + accel(positions.back()));

	 	positions.push_back(newPosition); velocity.push_back(newVelocity);
	 } while (positions.back() > 3); 

	 PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");
	 ExpansionCoeff coef(200, pd.maxRadialIndex()); 

	 double mass{0.01}; 
	 for (int i = 0; i < positions.size(); ++i) {
	 	for (int n = 0; n <= pd.maxRadialIndex(); ++n) {coef(i)(n) = (mass*positions[i]) * pd.potential(positions[i], n);} // include the position factor as the delta is in polar coords
	 }
	 for (int i = positions.size(); i < coef.nTimeStep(); ++i) {coef(i) = 0 * coef(0);}
	 coef.writePerturbation2File("Perturbations/infallingSatellite.out"); 
}



void circularInfall() {
	//savingEvolutionKernel(false, 0.15, 200, 0.1, "Kernels/kalnajs0.out", false); 
	inFallingCoefficents(0.1, 200);
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");

	std::vector<std::string> filenames = {"Plotting/Scar_Data/infallingRingSS.csv", "Plotting/Scar_Data/infallingRingST.csv", "Plotting/Scar_Data/infallingRingUS.csv", "Plotting/Scar_Data/infallingRingUT.csv"};
	std::vector<std::string> kernels = {"Kernels/kalnajsScarred0.out", "Kernels/kalnajsScarred0.out", "Kernels/kalnajs0.out", "Kernels/kalnajs0.out"};
	std::vector<bool> selfConsistent = {true, false, true, false};


	for (int i =0; i < filenames.size(); ++i) { 
		VolterraSolver solver1(kernels[i], 48, 2, 200, 0.1); // kalnajsScarred0.out
		solver1.activeFraction(0.35); 

		solver1.solveVolterraEquationPerturbation("Perturbations/infallingSatellite.out", selfConsistent[i]); 

		//solver1.density2dEvolution("Plotting/Scar_Data/infallingRingUS.csv", pd); 
		solver1.density1dEvolution(filenames[i], pd);
	}
}

// Anuglar Momentum Scars //
// ---------------------- //

void angularMomentumScar() {
	AngularMomentumScar scar(2, 0.1, 5);
	ScarredMestel<AngularMomentumScar> dfs;

	dfs.addScar(scar);

	std::ofstream out("Plotting/density.csv");

	for (double r = 0.01; r < 5; r += 0.05) {out << dfs.density(r) << '\n';}
	out.close();
}