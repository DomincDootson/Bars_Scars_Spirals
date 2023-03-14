#ifndef SPIRALEVOLUTIONFUNCTIONS
#define SPIRALEVOLUTIONFUNCTIONS

#include <string>

void savingInitialSpiral(const std::string & filename = "Plotting/spiral.csv");
void savingEvolutionKernel(double sigma, int nTimeStep, double timeStep, const std::string & filename , const bool generateBF); 

void densityEvolution(double k); 

void spiralEvolution(); 

void varyingKEvolution();
void varyingSigmaEvolution(); 
void smallSpiral(); 

void varyingNumberBasisFunctions();
void spiralGaussianTest();

#endif