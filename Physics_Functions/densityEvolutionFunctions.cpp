#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <cstring>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/Mestel.h"

#include "../Volterra_Solver/VolterraSolver.h"
#include "../Volterra_Solver/ExpansionCoeff.h"

#include "../Bar2D/Bar2D.h"

const int NUMBTIMESTEPS{200}; // some global paramters for the Volterra Solver
const double DELTAT{0.25};
double XI{.10};

bool SELF{true}; 

void makeSelfConsistent() {SELF = true;}
void makeTestParticle()   {SELF = false;}


std::string isTestParticle()
{
	if (!SELF)
	{
		return "_test";
	}
	return "";
}

std::string perturbationFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Perturbation_" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() + ".out";
}

std::string kickingCoefFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Coeff" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() +".csv";
}

std::string kickingDensityFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Density" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() +".csv";
}

std::string density2dFilename(const double littleSigma, const double radius, const int fourierHarmonic){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Density2D" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + isTestParticle() + ".csv";
}

std::string kernelFilename(const double littleSigma, const int fourierHarmonic, const std::string & bfName = "Kalnajs"){
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/" + bfName + "_" +std::to_string(fourierHarmonic) + ".out";
}



std::vector<double> littleSigmas(){
	return {0.25, 0.35, 0.50};
}

std::vector<double> perturbationRadii() {
	std::vector<double> holding;
	for (double radius = .5; radius <= 7; radius +=.1){
		holding.push_back(radius);
	}
	return holding;
}

void kalnajBF()
{
	Mestel DF;
	double Rka{15};
	std::vector<double> params{4, Rka};
	for (int i =0; i<=2; ++i){
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, i);

		ActionAngleBasisContainer test("Kalnajs", 10, i, 5, 101, 10);
		std::string file = "Kalnajs/Kalnajs_4_" + std::to_string((int) Rka); // Use file function name here
		test.scriptW(PD, DF, file);
	} 	
}


void kalnajsKernelsVaryingSigma(int l)
{
	std::vector<double> littleSigma{littleSigmas()};
	
	std::string file = "Kalnajs/Kalnajs_4_20";
	ActionAngleBasisContainer test(file, "Kalnajs", 10, l, 5, 101, 20);

	VolterraSolver solver(10, l, NUMBTIMESTEPS, DELTAT);
	for (int i = 0; i < littleSigma.size(); ++i){
		
		std::cout << "Calculating kernel for littleSigma: " << littleSigma[i] << '\n';
		Mestel DF(1, 1, littleSigma[i]);
		
		std::string kernel{kernelFilename(littleSigma[i], l)};
		
		solver.generateKernel(kernel, DF, test);
	}
} 


void maxDensityRadii(){
	std::ofstream out("../maxDensity.csv");
	std::vector<double> littleSigma{littleSigmas()};
	for (int i = 0; i < littleSigma.size()-1; ++i) {out << littleSigma[i] << ',';}
	out << littleSigma.back() << '\n';

	for (double rInner = .1; rInner < 5; rInner += .1){
		out << rInner << ',';
		std::cout << rInner << '\n';
		for (auto i = 0; i < littleSigma.size(); ++i){
			Mestel DF(1, 1, littleSigma[i], 1, rInner);
			if (i < littleSigma.size() -1){ out << DF.maxDensityRadius() << ',';}
			else {out << DF.maxDensityRadius() << '\n';}
		}
	}
	out.close(); 
}

// Perturbation Functions // 
// ---------------------- //

template <typename T>
void outPutPerturbation(const T & pd, double littleSigma, double radius, int numbTimeSteps = NUMBTIMESTEPS)
{
	Eigen::VectorXcd t0Coeff(pd.maxRadialIndex()+1);
	for (int i = 0; i <= pd.maxRadialIndex(); ++i){
		t0Coeff[i] = pd.potential(radius, i);
	}
	t0Coeff = - (pd.scriptE()).inverse() * t0Coeff; 
	ExpansionCoeff holding(t0Coeff, numbTimeSteps, pd.maxRadialIndex());
	holding.writePerturbation2File(perturbationFilename(littleSigma, radius, pd.fourierHarmonic()));
}
 
void cleaningPerturbations(const std::vector<double> & radii)
{
	std::vector<double> littleSigma{littleSigmas()}, angHarmonic{0,1,2};
	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto r = radii.begin(); r != radii.end(); ++r){
			for (auto l = angHarmonic.begin(); l != angHarmonic.end(); ++l){
				if (!remove(perturbationFilename(*s, *r, *l).c_str()))
				{
					std::cout << "Sucessfully deleted: " << perturbationFilename(*s, *r, *l) << '\n';
				}
			}
		}
	}
}

void diskKickingPerturbations()
{
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0(params, 10,0), pd1(params, 10, 1), pd2(params, 10, 2);

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto r = radii.begin(); r != radii.end(); ++r){
			outPutPerturbation(pd0, *s, *r);	
			outPutPerturbation(pd1, *s, *r);	
			outPutPerturbation(pd2, *s, *r);	
		}
	}
}
// SPIRAL EQU //
// ---------- //



void greensFunctions() {

	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_15_2.dat", 48);

	VolterraSolver solver("KalnajsNKernel.out", 48, 2, 100, 0.5); 

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(48+1); 

	for (int n =0; n <= pd.maxRadialIndex(); ++n) {coeff(n) = -pd.potential(2, n);}

	solver.activeFraction(0.5); 
	solver.setInitalPerturbation(coeff);
	solver.deltaPerturbationTest(); 
	solver.density2dEvolution("Plotting/Greens_Data/Kalnajs_Test_Green_2.csv", pd, 13, 5);


	VolterraSolver solverC("KalnajsNKernel.out", 48, 2, 100, 0.5); 
	solverC.activeFraction(0.5); 
	solverC.setInitalPerturbation(coeff);
	solverC.deltaPerturbationConsistent(); 
	solverC.density2dEvolution("Plotting/Greens_Data/Kalnajs_Green_2.csv", pd,13, 5);

}


// Coefficent Evolution //
// -------------------- //

void indivdualCoeffEvolution(const int nRadialHarmonics, const int angHarmonic, const double littleSigma, const double radius){
	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	std::string coeff{kickingCoefFilename(littleSigma, radius, angHarmonic)};
	std::string kernel{kernelFilename(littleSigma, angHarmonic)};
	std::cout << "Solving: " << littleSigma << " " << radius << " " << angHarmonic <<'\n';
	

	VolterraSolver solver(kernel, nRadialHarmonics , angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(XI);
	solver.coefficentEvolution(coeff, perturbationFile, SELF);
}

void coefficentEvolution()
{
	diskKickingPerturbations();
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto radius = radii.begin(); radius != radii.end(); ++radius){
			indivdualCoeffEvolution(10, 0, *s, *radius);
			indivdualCoeffEvolution(10, 1, *s, *radius);
			indivdualCoeffEvolution(10, 2, *s, *radius);
		}
	}
	cleaningPerturbations(radii);
}


// 2D Density Evolution //
// -------------------- //

void density2DKalnajs(double littleSigma, double radius, int angHarmonic){
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd{params, 10, angHarmonic};	

	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	outPutPerturbation(pd, littleSigma, radius);

	std::string kernel{kernelFilename(littleSigma, angHarmonic)};
	VolterraSolver solver(kernel, pd.maxRadialIndex(), angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(XI);

	solver.solveVolterraEquationPerturbation(perturbationFile, SELF);

	std::string outFilename{density2dFilename(littleSigma, radius, angHarmonic)};
	solver.density2dEvolution(outFilename, pd);
}

void density2DGaussian(double littleSigma, double radius, int angHarmonic) {
	std::vector<double> params{48, .5, 15};
	PotentialDensityPairContainer<GaussianLogBasis> pd{params, 48, angHarmonic};	

	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	outPutPerturbation(pd, littleSigma, radius);

	std::string kernel{"Kernels/gaussianComparison.out"};
	VolterraSolver solver(kernel, pd.maxRadialIndex(), angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(.40);

	solver.solveVolterraEquationPerturbation(perturbationFile, SELF);

	std::string outFilename{density2dFilename(littleSigma, radius, angHarmonic)};
	solver.density2dEvolution(outFilename, pd);
}


// Energy Evolution //
// ---------------- // 

void outPutEnergy(std::ofstream & out, const std::vector<double> & energies, double littleSigma, int angHarmonic, double radius){
	int skip{1};
	out << littleSigma << ',' << angHarmonic << ',' << radius << ',';
	for (int i = 0; i<energies.size()-skip; i += skip){
		out << energies[i] <<',';
	}
	out << energies[energies.size() - skip] << '\n';
}

template <class Tbf>
std::vector<double> individualEnergyEvolution(const Tbf & pd, double littleSigma, double radius, const std::string & bfName = "Kalnajs")
{
	int angHarmonic{pd.fourierHarmonic()};
	std::string perturbationFile{perturbationFilename(littleSigma, radius, angHarmonic)};
	std::string kernel{kernelFilename(littleSigma, angHarmonic, bfName)};
	std::cout << "Solving: " << littleSigma << " " << radius << " " << angHarmonic <<'\n';

	VolterraSolver solver(kernel, pd.maxRadialIndex(), angHarmonic, NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(XI);
	return solver.energyEvolution(pd, perturbationFile, SELF);
}


void energyEvolution(const std::string & energyFilename)
{
	diskKickingPerturbations();
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0{params, 10, 0}, pd1{params, 10, 1}, pd2{params, 10, 2};
	std::vector<PotentialDensityPairContainer<KalnajsBasis>> potentialDensityPairs{pd0, pd1, pd2};
	
	std::vector<double> littleSigma{littleSigmas()}, radii{perturbationRadii()};
	std::ofstream out(energyFilename);

	for (auto s = littleSigma.begin(); s != littleSigma.end(); ++s){
		for (auto pd = potentialDensityPairs.begin(); pd != potentialDensityPairs.end(); ++pd){
			for (auto r = radii.begin(); r != radii.end(); ++r){
				std::vector<double> energies = individualEnergyEvolution(*pd, *s, *r);
				outPutEnergy(out, energies, *s, pd -> fourierHarmonic(), *r);
			}
		}
	} 	
	out.close();
	cleaningPerturbations(radii);
}


// Basis Function Comparison // 
// ------------------------- //

double timeProfile(double time, double timescale = 50) {return sin(M_PI * (time/timescale));} 

void savePerturbation(const std::string & filename, const Eigen::VectorXcd & coeff, int nStep = NUMBTIMESTEPS, double timeStep = DELTAT) {
	ExpansionCoeff holding(nStep, coeff.size()+1);

	for (int time = 0; time < holding.nTimeStep(); ++time) {holding(time) = (timeProfile(time*timeStep)) * coeff;}
	holding.writePerturbation2File(filename); 
}
	
void saveVector(const std::string & filename, const std::vector<double> & vec) {
	std::ofstream out(filename);
	for (auto i : vec) {out << i << '\n';}
	out.close();
}
std::vector<double> radiiVector() {
	std::vector<double> radii;
 	for (int i =1; i < 2000; ++i) {radii.push_back(i*0.01);}
 	return radii;
 }


// Some function that call the above function and then the two d density evolution functiuns for both BF. 

void gaussianRepresentation (int nMax = 24, int index = 15) {
	/*std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> kalnajs(params, 10, 2);*/
	PotentialDensityPairContainer<KalnajsNBasis> kalnajs("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	
	
	std::vector<double> params1{static_cast<double>(nMax), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params1, nMax, 2);

 	Eigen::VectorXcd coefG = Eigen::VectorXcd::Zero(nMax+1);
 	//coefG(index-2) = -0.2;
 	//coefG(index-1) = -0.15;
 	coefG(index) = -0.1;
 	//coefG(index+1) = -0.07;
 	//coefG(index+1) = -0.03;
 	coefG *=0.03;
 
 	savePerturbation("Physics_Functions/gaussianPertGaussianRep.csv", coefG);
 	saveVector("Plotting/BF_Comparison/gaussianPerturbationG.csv", PD.oneDpotential(radiiVector(), coefG));
 	

 	Eigen::ArrayXXcd potentialG = PD.potentialArray(coefG, 1500, 20);
 	Eigen::VectorXcd coefK = kalnajs.potentialFitting(potentialG, 20);
 	savePerturbation("Physics_Functions/gaussianPertKalnajsRep.csv", coefK);
 	saveVector("Plotting/BF_Comparison/gaussianPerturbationK.csv", kalnajs.oneDpotential(radiiVector(), coefK));


 	/*potentialG = kalnajs.potentialArray(coefK, 1500, 20);
 	coefG = PD.potentialFitting(potentialG, 20);
 	savePerturbation("Physics_Functions/gaussianPertGaussianRep.csv", coefG);
 	saveVector("Plotting/BF_Comparison/gaussianPerturbationG.csv", PD.oneDpotential(radiiVector(), coefG));*/ 
}

void kalnajsRepresentation (int nMax = 24) {
	/*std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> kalnajs(params, 10, 2);*/
	PotentialDensityPairContainer<KalnajsNBasis> kalnajs("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	
	
	std::vector<double> params1{static_cast<double>(nMax), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params1, nMax, 2);

 	Eigen::VectorXcd coefK = Eigen::VectorXcd::Zero(kalnajs.maxRadialIndex() +1);
 	coefK(0) = 0.01;
 	savePerturbation("Physics_Functions/kalnajsPertKalnajsRep.csv", coefK);
 	saveVector("Plotting/BF_Comparison/kalnajsPerturbationK.csv", kalnajs.oneDpotential(radiiVector(), coefK));
	
 	Eigen::ArrayXXcd potentialK = kalnajs.potentialArray(coefK, 1500, 30);
 	Eigen::VectorXcd coefG(nMax+1);
 	coefG = PD.potentialFitting(potentialK, 30);	
 	
 	//coefG(nMax) = 1;
 	savePerturbation("Physics_Functions/kalnajsPertGaussianRep.csv", coefG);
 	saveVector("Plotting/BF_Comparison/kalnajsPerturbationG.csv", PD.oneDpotential(radiiVector(), coefG));	
}

// Shall we have a function that does the 2D evolution for a given perturbation? 

void kalnajs2D(const std::string & perturbationFile, const std::string & outFile) {
	/*std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> kalnajs(params, 10, 2);*/
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical.dat");

	VolterraSolver solver("Kernels/kalnajsComparison.out", pd.maxRadialIndex(), pd.fourierHarmonic(), NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(.50);
	solver.solveVolterraEquationPerturbation(perturbationFile, false);

	solver.density2dEvolution(outFile, pd); 
}

void guassian2D(const std::string & perturbationFile, const std::string & outFile, int nMax = 24) {
	std::vector<double> params1{static_cast<double>(nMax), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> pd(params1, nMax, 2);

 	VolterraSolver solver("Kernels/gaussianComparison.out", pd.maxRadialIndex(), pd.fourierHarmonic(), NUMBTIMESTEPS, DELTAT);
	solver.activeFraction(.50);
	solver.solveVolterraEquationPerturbation(perturbationFile, false);

	solver.density2dEvolution(outFile, pd); 
}

void gaussKernel() {
	Mestel DF; int nMax{48};
	
	/*std::vector<double> params{static_cast<double>(nMax), .15, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params, 48,2);*/ 

	ActionAngleBasisContainer test("GaussianLog", "GaussianLog", nMax, 2, 7, 251, 20); 
	//test.scriptW(PD, DF, "GaussianLog");

	VolterraSolver solver1(48, 2, 200, 0.25);
	std::string kernel1 = "Kernels/gaussianComparison.out";
	solver1.generateKernel(kernel1, DF, test);

}

void comparisonDensityEvolution() {
	//gaussianRepresentation(48, 40); 
	//gaussKernel(); 
	//kalnajsRepresentation(48); 
	
	kalnajs2D("Physics_Functions/kalnajsPertKalnajsRep.csv", "Plotting/BF_Comparison/kalnajsPertKalnajsRepDensity.csv");
	guassian2D("Physics_Functions/kalnajsPertGaussianRep.csv", "Plotting/BF_Comparison/kalnajsPertGaussianRepDensity.csv", 48); 

	//kalnajs2D("Physics_Functions/gaussianPertKalnajsRep.csv", "Plotting/BF_Comparison/gaussianPertKalnajsRepDensity.csv");
	//guassian2D("Physics_Functions/gaussianPertGaussianRep.csv", "Plotting/BF_Comparison/gaussianPertGaussianRepDensity.csv", 48);
}


void densityTest() {
	Mestel DF;
	std::ofstream out("Plotting/Density.csv");

	for (double radius = 0.01; radius < 15; radius += 0.05) {out << radius<<','<< DF.density(radius) << ',' << 1/(2*M_PI*radius) << '\n';}


}