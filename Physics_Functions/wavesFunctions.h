#ifndef WAVESFUNCTIONS
#define WAVESFUNCTIONS

void generatingKalnajsKernels(int m2, int nMax = 48);

void selfConsistentWaves(int nMode);
void perturbationWaves(int nMode);

void selfConsistentDensity(double radius);
void pullingDensity(double radius);

void selfConsistentPotential(double radius);
void pullingPotential(double radius);

void selfConsistentQuadDensity(double radius); 

#endif