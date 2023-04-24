#include <iostream>

#include "swingfunctions.h"


int main() {
	std::cout << "Let's get swinging\n";

	//ringEvolution("Plotting/test.csv", false);
	//ringEvolution("Plotting/selfconsistent.csv", true);

	//generateKernel(); 
	//discComparison("test.csv", false); 

	//generateSwingKernels(15); 
	//amplificationFixedRadius(5);

	//densityEvolutionFixedRadius(5);

	// generateSwingKernelsTemp(2, 1.3);
	// generateSwingKernelsTemp(4, 1.3);
	// generateSwingKernelsTemp(8, 1.3);
	densityEvolutionChi(5, 2);
	densityEvolutionChi(5, 4);
	densityEvolutionChi(5, 8);

	return 0; 
}