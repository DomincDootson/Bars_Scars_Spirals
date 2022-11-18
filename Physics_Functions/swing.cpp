#include <iostream>

#include "swingfunctions.h"


int main() {
	std::cout << "Let's get swinging\n";

	//ringEvolution("Plotting/test.csv", false);
	ringEvolution("Plotting/selfconsistent.csv", true);

	return 0; 
}