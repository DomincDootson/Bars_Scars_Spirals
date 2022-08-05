#ifndef NBODYFUNCTIONS
#define NBODYFUNCTIONS
#include <string>

void barEvolutionGaussian(); 
void barEvolutionKalnajs(const std::string & stem, const bool isSelfConsistent, const double littleSigma);

void checkingConservedJacobi();
void orbitSection();  

void kalanajTest();

void spiralTesting();
void calculateDiscAM();
#endif