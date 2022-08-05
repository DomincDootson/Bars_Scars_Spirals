#include "wavesFunctions.h"

int main() {
	/*selfConsistentWaves(10);
	perturbationWaves(10);

	selfConsistentWaves(20);
	perturbationWaves(20);

	selfConsistentWaves(30);
	perturbationWaves(30);*/ 
	
	selfConsistentDensity(2);
	selfConsistentDensity(5);
	selfConsistentDensity(10);
	selfConsistentDensity(13);

	pullingDensity(2);
	pullingDensity(5);
	pullingDensity(10);
	pullingDensity(13);

	// selfConsistentPotential(2);
	// selfConsistentPotential(5);
	// selfConsistentPotential(10);
	// selfConsistentPotential(13);

	// pullingPotential(2);
	// pullingPotential(5);
	// pullingPotential(10);
	// pullingPotential(13);

	selfConsistentQuadDensity(5); 

	//generatingKalnajsKernels(0, 38);
	return 0; 
}