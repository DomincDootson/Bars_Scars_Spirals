#ifndef GENERALFUNCTIONS
#define GENERALFUNCTIONS

#include <string>
void freqCheck();
void savingKalnajsFunctions(const std::string & filename); 
void generatingKalnajsBF(int m2);
void generatingSpiralBF(int m2);
void getSpiralParam();

void generatingKalnajsKernels(const std::string & filename, int m2, int nMax = 10, double rInner = 1);
void generatingGaussianKernels(int m2);

void generatingSpiralBFDiffTemp(int m2);


void testingBarTorque();
void testingFitting();

void saveKalnajs(); 

void kalnajTest();
void spiralAA();

void spiralBasisComparions();
void spiralTestAnalytic(); 
void spiralTestQuasi();
void spiralTestTrue();
void spiralEvolution();


void differentInnerTapers(); 
void energyTapping(int nMax, int rInner);

void waveTesting();
#endif
