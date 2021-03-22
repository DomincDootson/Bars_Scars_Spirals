#ifndef DENSITYEVOLUTIONFUNCTIONS
#define DENSITYEVOLUTIONFUNCTIONS

#include <string>

void makeSelfConsistent();
void makeTestParticle();

void kalnajBF();

void kalnajsKernelsVaryingSigma(int l);

void maxDensityRadii();

void plottingPerturbations();

void coefficentEvolution();

void diskKicking();
void density2D(double littleSigma, double radius, int angHarmonic);

void energyEvolution(const std::string & energyFilename);

#endif