#ifndef MODESFUNCTIONS
#define MODESFUNCTIONS 

#include <string>

void m1Stability();
void testingKernelMethod();

void generatingEigenMode(const std::string & modeFile, const std::string & densityFile, const std::string & residualFile, const std::string & timefile); 
void unstableTimeEvolution(const std::string & stem = "Coefficent_JB_Mode_");

void uniformSeach(const std::string & filename, int nu = 6);
void unstableSearch(const std::string & filename);

#endif