#ifndef WAVESFUNCTIONS
#define WAVESFUNCTIONS

#include <string>

void checkingDeltaFunctionFitting();

void waveTest(const std::string & filename); 
void waveTestSpinning(const std::string & filename);

void savingInitialWave();

void turnOffBarFile(const std::string & filename, double turnOffTime); 
void waveEvolutionTest(const std::string & densityfile, double CRposition, double perturberRadius, const std::string & growthFile = "Bar2D/barSizeTurnOff.out");

#endif