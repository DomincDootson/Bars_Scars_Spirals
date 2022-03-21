#include <iostream>

#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../DF_Class/Mestel.h"

#include "../Volterra_Solver/VolterraSolver.h"


#include "../Spiral2D/Spiral2D.h"

/* General Function */

void savingVector(std::ofstream & out, std::vector<double> & vec, double value) {
	out << value << ',';
	for (int i = 0; i < vec.size()-1; ++i){
		out << vec[i] << ',';
	}
	out << vec.back() << '\n';
}


void savingInitialSpiral(const std::string & filename) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	Spiral2D spiral(pd);
	spiral.density2dEvolution(pd, filename, 1); 
}


void savingEvolutionKernel(double sigma = 0.35, int nTimeStep = 200, double timeStep = 0.1, const std::string & filename = "Kernels/kalnajsN.out", const bool generateBF = false) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");	
	Mestel DF(1, 1, sigma); 

	if (generateBF) {
		ActionAngleBasisContainer test("KalnajsN", 48, 2, 7, 251, 20); test.scriptW(pd, DF, "KalnajsN");}
	
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", 48, 2, 7, 251, 20);

	VolterraSolver solver1(48, 2, nTimeStep, timeStep);
	solver1.generateKernel(filename, DF, test);
}


/* Density Evolution */
/* ----------------- */ 

void densityEvolution(double k) {
	//savingEvolutionKernel(.35, 200, 0.25); 
	VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 800, 0.25);
	solver.activeFraction(.4);

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	Spiral2D spiral(pd, k);

	solver.spiralEvolution(spiral); 

	spiral.density2dEvolution(pd, "Plotting/spiral.csv", 1, 10) ;
}




/* Code Varying */
/* ------------ */


void varyingKEvolution() {
	//savingEvolutionKernel(.35, 250, 0.25, "Kernels/kalnajsN.out", false); 
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	

std::vector<double> kVector = {.5, .75, 1, -.5, -.75, -1}; 

	for (auto k : kVector) { 
		auto filename = [k] () {return "Plotting/Spiral_Data/Spiral_" + std::to_string((int) round(100*k)) + ".csv";};
		std::cout << "Saving to: " << filename() << '\n';
		std::cout << '\n' << "Evolution for k = " << k << '\n' << '\n';
		VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 250, 0.25);
		solver.activeFraction(.4);
		
		Spiral2D spiral(pd, k, 5);
		solver.spiralEvolution(spiral); 

		//std::vector<double> vec = spiral.densityPower(pd);
		//savingVector(out, vec, k);

		
		spiral.density2dEvolution(pd, filename(), 1, 10);
	}
	//out.close(); 
}


void varyingSigmaEvolution() {

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	
	std::ofstream out("Plotting/VaryingSpiralTemperature.csv");

	std::vector<double> sigma = {.2, .30, .40}; 

	for (auto s : sigma) { 
		std::cout << '\n' << "Evolution for sigma = " << s << '\n' << '\n';
		savingEvolutionKernel(s, 200, 0.25); 
		VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 200, 0.25);
		solver.activeFraction(.4);
		
		Spiral2D spiral(pd, -1);
		solver.spiralEvolution(spiral); 

		std::vector<double> vec = spiral.densityPower(pd);
		

		savingVector(out, vec, s);
	}
	out.close(); 
}

void smallSpiral() {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	double taper{30};
	Spiral2D spiral(pd, -1, 5, 0.01, taper);

	VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 250, 0.25);
	solver.activeFraction(.5);
	solver.spiralEvolution(spiral); 

	spiral.removeIC(); 
	spiral.density2dEvolution(pd, "Plotting/Spiral_Data/TightlyWoundNoIC.csv", 1, 10);
}
