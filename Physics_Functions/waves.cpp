#include "wavesFunctions.h"
#include <cmath>

int main() {
	//scarredTappingCoefficents("test");
	//checkingDeltaFunctionFitting(); 

	//savingInitialWave();
	//waveEvolutionTest();
	//waveTest("Plotting/IngoingOutgoingWave.csv");
	//waveTestSpinning("Plotting/Spinning_wave.csv"); 

	//turnOffBarFile("Bar2D/barSizeTurnOff.out", 1);

	
	//waveEvolutionTest("Plotting/Waves_Data/Stirring/CR_5_ILR_R_P.csv", 5, 5 * (1-0.5*sqrt(2)));
	//waveEvolutionTest("Plotting/Waves_Data/Stirring/CR_5_CR_R_P.csv", 5, 5);
	waveEvolutionTest("Plotting/Test_density.csv", 5, 5 * (1-0.5*sqrt(2)), "Bar2D/barSizeTurnOnSlow.out");

	//waveEvolutionTest("Plotting/Waves_Data/Stirring/CR_5_ILR_W_P.csv", 5, 5 * (1-0.5*sqrt(2)), "None");
	//waveEvolutionTest("Plotting/Waves_Data/Stirring/CR_5_CR_W_P.csv", 5, 5, "None");
	
	//waveEvolutionTest("Plotting/Waves_Data/Stirring/CR_5_OLR_W_N.csv", 5, 5 * (1+0.5*sqrt(2)), "None");
	return 0; 
	
}