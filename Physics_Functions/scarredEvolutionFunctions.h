#ifndef SCARREDEVOLUTIONFUNCTIONS
#define SCARREDEVOLUTIONFUNCTIONS 

#include <string>

void savingEvolutionKernel(double scarRadius, int nTimeStep = 200, double timeStep = 0.1, const std::string & filename = "Kernels/kalnajsScarred.out", const bool generateBF = false);

void circularCutThrought(); 
void scarredDensity(double temp = 0.35, const std::string & filename = "Plotting/Scar_Data/scarredDensity.csv");
void scarredPotential();
void backgroundDF();
 
void testingScaredMestelEvolution();
void scarredModes();

void circularInfall(); 
void angularMomentumScar();

void savingDensityFile(double scarRadius, const std::string & filename);
#endif