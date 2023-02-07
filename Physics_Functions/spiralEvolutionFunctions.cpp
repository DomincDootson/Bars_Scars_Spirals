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
	double k{1};
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	Spiral2D spiralLeading(pd, k);
	spiralLeading.density2dEvolution(0, pd, "Plotting/Spiral_Data/spiralLeading.csv", 10); 

	Spiral2D spiralTrailing(pd, -k);
	spiralTrailing.density2dEvolution(0, pd, "Plotting/Spiral_Data/spiralTrailing.csv", 10); 

	std::ofstream out("Plotting/Spiralden1D.csv");

	for (double radius = 0.01; radius < 10; radius += 0.05) {out << radius <<',' << spiralTrailing.density(radius, 0) << ',' << spiralLeading.density(radius, 0) <<'\n';}
	out.close(); 
}


void savingEvolutionKernel(double sigma = 0.35, int nTimeStep = 200, double timeStep = 0.1, const std::string & filename = "Kernels/kalnajsN.out", const bool generateBF = false, const int nMax = 48) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat", nMax);	
	Mestel DF(1, 1, sigma); 

	if (generateBF) {
		ActionAngleBasisContainer test("KalnajsN", nMax, 2, 7, 251, 20); test.scriptW(pd, DF, "KalnajsN");}
	
	ActionAngleBasisContainer test("KalnajsN", "KalnajsN", nMax, 2, 7, 251, 20);

	VolterraSolver solver1(nMax, 2, nTimeStep, timeStep);
	solver1.generateKernel(filename, DF, test);
}


/* Density Evolution */
/* ----------------- */ 

void densityEvolution(double k) {
	//savingEvolutionKernel(.35, 200, 0.25); 
	//VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 800, 0.25);
	//solver.activeFraction(.4);

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	Spiral2D spiral(pd, k);

	//solver.spiralEvolution(spiral); 

	//spiral.density2dEvolution(pd, "Plotting/spiral.csv", 1, 10) ;
}


/* To Compare with N body */ 
/* ---------------------- */

void spiralEvolution() {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", 40);
	Spiral2D spiral(pd, 1, 5, 1);

	VolterraSolver solver("Kernels/Swing_Kernels/SW_Ker_2.csv", 40, 2, 500, 0.3);
	solver.activeFraction(.5);

	solver.spiralEvolution(spiral); 
	solver.density2dEvolution("Plotting/test_spiral.csv", pd, 1, 15);
	
}



/* Code Varying */
/* ------------ */


void varyingKEvolution() {
	//savingEvolutionKernel(.35, 500, 0.50, "Kernels/kalnajsN.out", false); 
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	
	//std::ofstream out("Plotting/Spiral_Data/VaryingK.csv"); 
	std::vector<double> kVector = {-.5, -.75, -1, .5, .75, 1};  //  

	for (auto k : kVector) { 
		auto filename = [k] () {return "Plotting/Spiral_Data/Spiral_" + std::to_string((int) round(100*k)) + ".csv";};
		std::cout << "Saving to: " << filename() << '\n';
		std::cout << '\n' << "Evolution for k = " << k << '\n' << '\n';
		VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 500, 0.50);
		solver.activeFraction(.4);
		
		Spiral2D spiral(pd, k, 5);
		solver.spiralEvolution(spiral); 

		//std::vector<double> vec = spiral.densityPower(pd);
		//savingVector(out, vec, k);
		
		//spiral.density2dFinal(pd, filename(), 250);
		spiral.density2dEvolution(pd, filename(), 1, 10);
	}
	//out.close(); 
}


void varyingSigmaEvolution() {

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");
	
	std::ofstream out("Plotting/Spiral_Data/VaryingSpiralTemperature.csv");

	std::vector<double> sigma = {.20, .25, .30, .35}; 

	for (auto s : sigma) { 
		std::cout << '\n' << "Evolution for sigma = " << s << '\n' << '\n';
		savingEvolutionKernel(s, 400, 0.25, "Kernels/kalnajsTempN.out"); 
		auto filename = [s] () {return "Plotting/Spiral_Data/Spiral_Sigma_" + std::to_string((int) round(100*s)) + ".csv";};
		VolterraSolver solver("Kernels/kalnajsTempN.out", 48, 2, 400, 0.25);
		solver.activeFraction(.4);
		
		Spiral2D spiral(pd, 1);
		solver.spiralEvolution(spiral); 

		std::vector<double> vec = spiral.densityPower(pd);

		//spiral.density2dFinal(pd, filename());

		savingVector(out, vec, s);
	}
	out.close(); 
}

std::string smallSpiralFilename(double taper, double k) {
	if (k > 0) {return "Plotting/Spiral_Data/TightlyWound_" + std::to_string((int) taper) + "_trailing.csv";}
	else       {return "Plotting/Spiral_Data/TightlyWound_" + std::to_string((int) taper) + "_leading.csv";}
}

void smallSpiral() {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");


	std::vector<double> tapers{10, 20, 30, 10, 20, 30}, ks{1, 1, 1, -1, -1, -1};

    for (int i = 0; i < tapers.size(); ++i) {		
			Spiral2D spiral(pd, ks[i], 5, 0.01, tapers[i]);
		
			VolterraSolver solver("Kernels/kalnajsN.out", 48, 2, 250, 0.25);
			solver.activeFraction(.5);
			solver.spiralEvolution(spiral); 
		
			//spiral.removeIC(); 
			spiral.density2dEvolution(pd, smallSpiralFilename(tapers[i], ks[i]), 1, 10);}
}


void varyingNumberBasisFunctions() {
	
	std::vector<double> nMax{24, 34, 48};

	std::ofstream out("Plotting/Spiral_Data/VaryingK.csv"); 
	  

	for (auto n : nMax) { 
		PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat", n);
		std::cout << pd.maxRadialIndex() << '\n';
		savingEvolutionKernel(0.35, 500, 0.5, "Kernels/kalnajsN.out", false, n);
		auto filename = [n] () {return "Plotting/Spiral_Data/Spiral_" + std::to_string((int) round(n)) + ".csv";};
		std::cout << "Saving to: " << filename() << '\n';
		std::cout << '\n' << "Evolution for nMax = " << n << '\n' << '\n';
		VolterraSolver solver("Kernels/kalnajsN.out", n, 2, 500, 0.50);
		solver.activeFraction(.4);
		
		Spiral2D spiral(pd, 1, 5);
		solver.spiralEvolution(spiral); 

		std::vector<double> vec = spiral.densityPower(pd);
		savingVector(out, vec, n);
		
		spiral.density2dFinal(pd, filename(), 250);
		//spiral.density2dEvolution(pd, filename(), 1, 10);
		
	}
	out.close(); 
} 
