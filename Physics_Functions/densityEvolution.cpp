#include "densityEvolutionFunctions.h"
#include <iostream>

int main()
{
	
	//maxDensityRadii();

	// First Generate the new Kernels!
	
	plottingPerturbations();
	/*kalnajBF();
	kalnajsKernelsVaryingSigma(0);
	kalnajsKernelsVaryingSigma(1);
	kalnajsKernelsVaryingSigma(2);
	
	makeSelfConsistent();
	energyEvolution("Disk_Kicking/Energy_Evolution/KalnajsEnergyEvolution_10_15.csv");
	makeTestParticle();
	energyEvolution("Disk_Kicking/Energy_Evolution/KalnajsEnergyEvolution_10_15_Test.csv");*/
	//coefficentEvolution();
	//diskKicking();
	//energyEvolution("Disk_Kicking/Energy_Evolution/KalnajsEnergyEvolution_01_10.csv");
	//energyEvolution("Disk_Kicking/Energy_Evolution/KalnajsEnergyEvolution_05_10.csv");
	//energyEvolution("Disk_Kicking/Energy_Evolution/KalnajsEnergyEvolution_1_10.csv");
	
	//energyEvolution("Disk_Kicking/Energy_Evolution/KalnajsEnergyEvolution_2_10.csv");
	return 0;
}
// Could we write a function that would clean up the perturbation files once they have been used? 