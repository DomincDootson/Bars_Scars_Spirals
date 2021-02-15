#ifndef PHYSICS
#define PHYSICS 

#include <string>
// Could we possibly put comments to tell us what each of these functions do? 

void gaussianScriptE(int m2); 

// Basis Function Generation

void generatingKalnajsBF(int m2);
void generatingSpiralBF(int m2);

void kalnajBFVaryingK();
void kalnajBFVaryingR();


// Kernel Generation
void kernelFlipped(); // Saves a flipped kernel so that it can be used in the old code

void generatingKalnajsKernels(int m2);
void generatingGaussianKernels(int m2);

void kalnajsKernelsVaryingK();
void kalnajsKernelsVaryingR();
void kalnajsKernelsVaryingSigma(int l);

// Provides some test evolutions

void testEvolutionKalanajs(int m2);

void testEvolutionGaussian(int m2); // Does evolution for singular l harmonic
void spiralTestEvolution(); // Does evolution for the test plot in paper 1


// Paper 1 Functions //
void barKickingPerturbations();


// Paper 2 functions

void barVaryingAngularSpeed(); // Genereates data for different rotational speeds of bars

void barVaryingKka();
void barVaryingRka();




void kalnajBasisFunctionsVaryingK(); // 


// Misc 
void savingPotentialArrays(std::string dir);




#endif