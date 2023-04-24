#ifndef SWINGFUNCTIONS
#define SWINGFUNCTIONS
#include <string>


void ringEvolution(const std::string & filename, bool isSelfConsistent = true);
void generateKernel();

void discComparison(const std::string & filename, bool isSelfConsistent); 

void generateSwingKernels(int lMax = 15);
void amplificationFixedRadius(double rad, int startHarmonic = 0, int endHarmonic = 15);

void densityEvolutionFixedRadius(double rad, int startHarmonic = 0, int endHarmonic = 15);
void generateSwingKernelsTemp(int l, double Q);
void densityEvolutionChi(double rad, int harmonic);
#endif