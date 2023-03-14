#ifndef MODESFUNCTIONS
#define MODESFUNCTIONS 

#include <string>

void m1Stability();
void testingKernelMethod();

void saveResponseMatrix(const std::string & filename, const std::string & kernelFile, const double omega0, const double eta, const double xi = 0.5);
void eigenvector2Density(const std::string & modeFile, const std::string & densityFile); 
void eigenModeEvolution(const std::string & modeFile, const double omega0, const double eta, const std::string & stem);

void unstableTimeEvolution(const std::string & stem = "Coefficent_JB_Mode_");

void uniformSearchMatrix(const std::string & filename, int nu = 4, double eLower = 0.09, double oLower = 0.6);
void uniformSearchKernel(const std::string & filename, int nu = 4, double eLower = 0.09, double oLower = 0.6);
void uniformSearch(const std::string & filestem, const int taper);

void uniformSearchScarredDensity(const std::string & radius, const std::string & width, const std::string & depth, bool isLong = false); 
void fittingMethod(const std::string & outFilename = "Plotting/Modes_Data/fitting_test.csv", const std::string & kernelFile = "Kernels/Kalnajs_6_.out", double xi = 1);
void fittingChi();

void saveModeEvolution(const std::string & filename, const std::string & eigenvalues, double omega, double timeStep = 0.25);
#endif