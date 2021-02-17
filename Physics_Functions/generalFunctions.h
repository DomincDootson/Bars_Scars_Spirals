#ifndef GENERALFUNCTIONS
#define GENERALFUNCTIONS

#include <string>

void gaussianScriptE(int m2); 

void generatingKalnajsBF(int m2);
void generatingSpiralBF(int m2);

void generatingKalnajsKernels(int m2);
void generatingGaussianKernels(int m2);
void kernelFlipped(); // Saves a flipped kernel so that it can be used in the old code

void testEvolutionKalanajs(int m2);
void testEvolutionGaussian(int m2); // Does evolution for singular l harmonic
void spiralTestEvolution(); // Does evolution for the test plot in paper 1

void savingPotentialArrays(std::string dir);

#endif