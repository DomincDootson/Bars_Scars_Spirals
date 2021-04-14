#ifndef DENSITYEVOLUTIONFUNCTIONS
#define DENSITYEVOLUTIONFUNCTIONS

#include <string>

void makeSelfConsistent();
void makeTestParticle();

void kalnajBF();

void kalnajsKernelsVaryingSigma(int l);
void GaussianLogKernelsVaryingSigma(int l, double rInner, double rOuter, const std::string & dir = "GaussianLog");

void maxDensityRadii();

void plottingPerturbations();

void coefficentEvolution();

void diskKicking();
void density2D(double littleSigma, double radius, int angHarmonic);

void energyEvolution(const std::string & energyFilename);


// Gaussian Functions //

void generatingSpiralBF(const std::string & dir, const double innerTaper, const double outerTaper);
void diskKickingLGEnergy(const std::string & energyFilename);

#endif