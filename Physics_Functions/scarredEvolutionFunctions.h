#ifndef SCARREDEVOLUTIONFUNCTIONS
#define SCARREDEVOLUTIONFUNCTIONS 

#include <string>

void circularCutThrought(); 
void scarredDensity(double temp = 0.35, const std::string & filename = "Plotting/Scar_Data/scarredDensity.csv");
void scarredPotential();
void backgroundDF();
 
void testingScaredMestelEvolution();
void scarredModes();

void circularInfall(); 
void angularMomentumScar();
#endif