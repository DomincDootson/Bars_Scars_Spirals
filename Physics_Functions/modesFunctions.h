#ifndef MODESFUNCTIONS
#define MODESFUNCTIONS 

#include <string>

void m1Stability();

void generatingEigenMode(const std::string & modeFile, const std::string & outFile); 
void unstableTimeEvolution(const std::string & densityFile, const std::string & coeffFile);

void uniformSeach();

#endif