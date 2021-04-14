#include "densityEvolutionFunctions.h"
#include <iostream>

int main()
{
	
	//maxDensityRadii();

	// First Generate the new Kernels!
	
	// first we need the new BF

	/*generatingSpiralBF("GaussianLog_20_150", 2, 15);
	generatingSpiralBF("GaussianLog_5_150", 2, 15);
	generatingSpiralBF("GaussianLog_10_100", 1, 10);
	generatingSpiralBF("GaussianLog_10_175", 1, 17.5);*/



	GaussianLogKernelsVaryingSigma(0, 2, 15);
	GaussianLogKernelsVaryingSigma(1, 2, 15);
	GaussianLogKernelsVaryingSigma(2, 2, 15);
	diskKickingLGEnergy("GaussianLogEnergy_20_150.csv");

	
	GaussianLogKernelsVaryingSigma(0, 0.5, 15);
	GaussianLogKernelsVaryingSigma(1, 0.5, 15);
	GaussianLogKernelsVaryingSigma(2, 0.5, 15);
	diskKickingLGEnergy("GaussianLogEnergy_5_150.csv");

	
	GaussianLogKernelsVaryingSigma(0, 1, 10);
	GaussianLogKernelsVaryingSigma(1, 1, 10);
	GaussianLogKernelsVaryingSigma(2, 1, 10);
	diskKickingLGEnergy("GaussianLogEnergy_10_100.csv");

	
	GaussianLogKernelsVaryingSigma(0, 1, 17.5);
	GaussianLogKernelsVaryingSigma(1, 1, 17.5);
	GaussianLogKernelsVaryingSigma(2, 1, 17.5);
	diskKickingLGEnergy("GaussianLogEnergy_10_175.csv");






	//plottingPerturbations();
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