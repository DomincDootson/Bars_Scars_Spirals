#include <iostream>

#include "swingfunctions.h"


int main() {
	std::cout << "Let's get swinging\n";

	//ringEvolution("Plotting/test.csv", false);
	//ringEvolution("Plotting/selfconsistent.csv", true);

	//generateKernel(); 
	//discComparison("test.csv", false); 

	generateSwingKernels(15); 
	amplificationFixedRadius(5);

	//densityEvolutionFixedRadius(5);

	return 0; 
}