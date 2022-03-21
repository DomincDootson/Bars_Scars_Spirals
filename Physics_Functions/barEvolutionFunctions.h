#ifndef BAREVOLUTIONFUNCTIONS
#define BAREVOLUTIONFUNCTIONS

#include <Eigen/Dense>


void kalnajsTorque(int nMax);

void kalnajsBarTest();

void kalnajBFVaryingK();
void kalnajBFVaryingR();

void kalnajsKernelsVaryingK();
void kalnajsKernelsVaryingR();


void barVaryingAngularSpeed(); 
void barVaryingKka();
void barVaryingRka();

void barVaryingTurnOn();
void barVaryingActiveFraction();

void kalnajsKernelsDiffTemp(); 

void gaussianBarEvolution();
void guassianDifferentTemps(); 
void gaussianBarDiffGrowthRate();

void gaussianBarModel();
void kalnajsBarModel();
#endif