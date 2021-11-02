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
	

	//GaussianLogKernelsVaryingSigma(2, 1, 10);
	//density2DGaussian(.25, 1.95, 2);

	comparisonDensityEvolution();

	//density2DGaussian(0.35, 2, 2);





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