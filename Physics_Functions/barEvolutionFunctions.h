#ifndef BAREVOLUTIONFUNCTIONS
#define BAREVOLUTIONFUNCTIONS

#include <Eigen/Dense>


void kalnajsTorque(int nMax);

void kalnajsBarTest();

void barVaryingTurnOn();

void kalnajsKernelsDiffTemp(); 

void saveFittedSormani(const std::string & filename, const std::string & barFile = "Bar2D/Bar_Potentials/Sormani_Medium.out");

void sormaniBarEvolution(const std::string & barFile = "Bar2D/Bar_Potentials/Sormani_Large.out");
void differentTempKernels(); 


void sormaniPatternSpeed(const std::string & dir = "Plotting/Bar_Data/Pattern_Speed/"); 
void sormaniBarShape(const std::string & dir = "Plotting/Bar_Data/Kalnajs_Shape/");
void sormaniConstantRatio(const std::string & dir = "Plotting/Bar_Data/Sormani_Ratio/");



#endif