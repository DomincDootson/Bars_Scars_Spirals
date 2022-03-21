#include <cmath>
#include <iostream>
#include <fstream>

#include <chrono>

#include <boost/math/special_functions/ellint_1.hpp>

#include "../DF_Class/Mestel.h"
#include "../DF_Class/ScarredMestel.h"
#include "../DF_Class/Scar.h"

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../Volterra_Solver/VolterraSolver.h"

void addScars(ScarredMestel & df, bool scar1 = true, bool scar2 = true, bool scar3 = true) { 
	if (scar1) {Scar scar1(1.20, 11, 0.29, 0.474, -0.044, -5.219, -44.796); df.addScar(scar1);}
	//if (scar2) {Scar scar2(1.55,  6, 0.10, 0.562,  0.288, -3.371, -15.619); df.addScar(scar2);}
	//if (scar3) {Scar scar3(1.90,  6, 0.12, 0.081,  0.869, -0.486,  -3.242); df.addScar(scar3);}

}

void savingEvolutionKernel(bool toAddScars, double sigma = 0.35, int nTimeStep = 200, double timeStep = 0.1, const std::string & filename = "Kernels/kalnajsScarred.out", const bool generateBF = false) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat"); int l{0};
	ScarredMestel DF(1, 1, sigma); 
	if (toAddScars) {addScars(DF);}

	if (generateBF) {
		ActionAngleBasisContainer test("KalnajsN", 48, l, 7, 251, 20); test.scriptW(pd, DF, "KalnajsN");}
	
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, l, 7, 251, 20);

	VolterraSolver solver1(48, l, nTimeStep, timeStep);
	solver1.generateKernel(filename, DF, test);
}



void circularCutThrought() {
	ScarredMestel df; 

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

void scarredDensity() {
	ScarredMestel df(1,1,0.284); 

	std::ofstream out("Plotting/Scar_Data/unscarredDensity.csv");
	for (double r = 0.1; r<10; r += 0.01) {std::cout << "Density for point: " << r << '\n'; out << r << ',' << df.density(r) << '\n';}
	out.close(); 

	addScars(df); 
	std::ofstream out1("Plotting/Scar_Data/scarredDensity.csv");
	for (double r = 0.1; r<10; r += 0.01) {std::cout << "Density for point: " << r << '\n'; out1 << r << ',' << df.density(r) << '\n';}
	out1.close();  
}

void scarredPotential() {
	
	ScarredMestel df(1,1,0.35);
 
	std::ofstream out("Plotting/Scar_Data/unscarredPotential.csv");
	for (double r = 0.1; r<15; r += 0.01) {std::cout << "Potential for point: " << r << '\n'; out << r << ',' << df.potentialFromVector(r) << '\n';}
	out.close(); 

	Scar scar1(1.20, 11, 0.29, 0.474, -0.044, -5.219, -44.796); 
	ScarredMestel dfS({scar1}, 1,1,0.35);
 
	std::ofstream out1("Plotting/Scar_Data/scarredPotential.csv");
	for (double r = 0.1; r<15; r += 0.01) {std::cout << "Potential for point: " << r << '\n'; out1 << r << ',' << dfS.potentialFromVector(r) << '\n';}
	out1.close(); 

}


void backgroundDF() {
	ScarredMestel DF(1, 1, 0.35); 
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
	VolterraSolver solver1("Kernels/scarredMestel.out", 10, 2, 500, 0.25);
	solver1.activeFraction(0.06); 
	solver1.kernelTesting("Plotting/scarredEvolution.csv", 10, true);
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
	//savingEvolutionKernel(true, 0.35, 200, 0.1, "Kernels/kalnajsScarred0.out", false); 
	inFallingCoefficents(0.1, 200);
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");

	VolterraSolver solver1("Kernels/kalnajs0.out", 48, 2, 200, 0.1); // kalnajsScarred0.out
	solver1.activeFraction(0.06); 

	solver1.solveVolterraEquationPerturbation("Perturbations/infallingSatellite.out", true); 

	solver1.density2dEvolution("Plotting/Scar_Data/infallingRingUS.csv", pd); 

}


