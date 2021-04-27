#include <iostream>
#include <typeinfo>
#include <Eigen/Dense>
#include <vector>
#include <complex>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/Mestel.h"

#include "../Volterra_Solver/VolterraSolver.h"
#include "../Volterra_Solver/ExpansionCoeff.h"

#include "../Bar2D/Bar2D.h"

#include "barEvolutionFunctions.h"


std::string kernelName(std::string dir, std::string stem, int Kka, int Rka, int m2)
{
	return dir + "/" + stem + "_" + std::to_string(Kka) + "_" + std::to_string(Rka) + "_" + std::to_string(m2) + ".out";
}

std::string kernelName(std::string dir, std::string stem, int Kka, double Rka, int m2) 
{
	return dir + "/" + stem + "_" + std::to_string(Kka) + "_" + std::to_string((int) round(Rka)) + "_" + std::to_string(m2) + ".out";
}

void kalnajsBarTest() {


	VolterraSolver solver("test200.out", 10, 2, 200, 0.1);

	solver.activeFraction(.25);	
	
	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	coef[0] = .01; 
	
	Bar2D bar(coef, 0.5);
	solver.barRotation(bar, "linearCoeff.csv", "linearEvolution.csv", false, false); 
	
}

void kalnajBFVaryingK()
{
	Mestel DF;
	std::vector<double> Kka{4, 5, 6, 7};
	for (int i = 0; i < 4; ++i){
		
		std::cout << "Calculationg BF for Kka: " << Kka[i] << '\n';

		std::vector<double> params{Kka[i], 10};
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2);

		ActionAngleBasisContainer test("Kalnajs", 10, 2, 5, 101, 10);
		std::string file = "Kalnajs/Kalnajs_" + std::to_string((int) Kka[i]) + "_10"; // Use file function name here
		test.scriptW(PD, DF, file);
	}
}

void kalnajBFVaryingR()
{
	Mestel DF;
	std::vector<double> Rka{5, 10, 15, 20};
	for (int i = 0; i < 4; ++i){
		
		
		std::cout << "Calculationg BF for Rka: " << Rka[i] << '\n';
		std::vector<double> params{4, Rka[i]};
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2);

		ActionAngleBasisContainer test("Kalnajs", 10, 2, 5, 101, 10);
		std::string file = "Kalnajs/Kalnajs_4_" + std::to_string((int) Rka[i]); // Use file function name here
		test.scriptW(PD, DF, file);
	}
}


void kalnajBasisFunctionsVaryingK() 
{
	std::vector<double> params{4, 10};
	std::ofstream out("quad_Density_Functions.csv");

	for (int i = 0; i<999; ++i)
	{
		out << i*0.01 << ',';
	}
	out << 999*0.01 << '\n';

	for (int k = 3; k <9;++k)
	{
		params[0] = k;
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2); // Remake the plots
		for (int i = 0; i<999; ++i)
		{
			out << PD.density(i*0.01,0) << ',';
		}
		out << PD.density(999*0.01,0)  << '\n';

	}
	out.close();
}

std::string evolutionFileName(std::string dir, double omega0)
{
	return dir + "/evolution" +std::to_string((int) round(100*omega0)) + ".csv";
}

std::string coeffFileName(std::string dir, double omega0)
{
	return dir + "/coeff" +std::to_string((int) round(100*omega0)) + ".csv";		
}

void barVaryingAngularSpeed()
{
	std::vector<double> angSpeed = {.04, .08, .12, .16, .20};

	std::string kernelFileName = kernelName("Kernels", "Kalnajs", 4, 10, 2);
	VolterraSolver solver(kernelFileName, 10, 2, 2000, 0.025);

	solver.activeFraction(.25);	
	



	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	coef[0] = .1; 
	
	for (int i =0; i < angSpeed.size(); ++i){
		// Construct the bar
		Bar2D bar(coef, angSpeed[i]);
		std::string evolutionFilename = evolutionFileName("BarEvolution/VaryingOmega", angSpeed[i]);
		std::string outFilename = coeffFileName("BarEvolution/VaryingOmega", angSpeed[i]);
		solver.barRotation(bar, outFilename, evolutionFilename, true, true); 
	}
}

void barVaryingKka()
{
	std::vector<int> Kka{4, 5, 6, 7};
	for (int i = 0; i < 4; i++)
	{
		std::cout << "Evolution for Kka: " << Kka[i] << '\n';
		std::string kernelFileName = kernelName("Kernels", "Kalnajs", Kka[i], 10, 2);
		VolterraSolver solver(kernelFileName, 10, 2, 2000, 0.025);
		solver.activeFraction(.25);	
		 
		Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
		coef[0] = .1;
		Bar2D bar(coef, 0);

		std::string outFilename = coeffFileName("BarEvolution/VaryingKka", 0.01*Kka[i]);
		std::string evolutionFilename = evolutionFileName("BarEvolution/VaryingKka", 0.01*Kka[i]);

		solver.barRotation(bar, outFilename, evolutionFilename); 
	}
}

void barVaryingRka()
{
	std::vector<double> Rka{5, 10, 15, 20};
	for (int i = 0; i < 4; i++)
	{
		std::cout << "Evolution for Rka: " << Rka[i] << '\n';
		std::string kernelFileName = kernelName("Kernels", "Kalnajs", 4, Rka[i], 2);
		VolterraSolver solver(kernelFileName, 10, 2, 2000, 0.025);
		solver.activeFraction(.25);	
		 
		Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
		coef[0] = .1;
		Bar2D bar(coef, 0);

		std::string outFilename = coeffFileName("BarEvolution/VaryingRka", 0.01*Rka[i]);
		std::string evolutionFilename = evolutionFileName("BarEvolution/VaryingRka", 0.01*Rka[i]);

		solver.barRotation(bar, outFilename, evolutionFilename); 
	}
}


// Gaussian Bar //

Eigen::VectorXcd gaussianResolving(const double radius)
{
	int nMax{24};
	std::vector<double> params{static_cast<double>(nMax), .5, 15};
	PotentialDensityPairContainer<GaussianLogBasis> basisFunctions(params, nMax, 2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(nMax+1);

	for (int i = 0; i <= nMax; ++i){
		coeff(i) = basisFunctions.potential(radius, i);
	}
	coeff = -(basisFunctions.scriptE()).inverse()*coeff;

	// The following code save the reconstruciton of the density and the potential 
	/*std::vector<double> radii;
	for (double r = .1; r<5; r += .01){radii.push_back(r);}

	std::vector<double> density =  basisFunctions.oneDdensity(radii, coeff);
	std::vector<double> potential =  basisFunctions.oneDpotential(radii, coeff);

	std::ofstream out("Plotting/density.csv");
	for (int i = 1; i < density.size(); ++i){
		out << radii[i] << ',' <<  density[i] <<','<< potential[i] << '\n';
	}

	Eigen::ArrayXXd density = basisFunctions.densityArrayReal(coeff, 500, 5);
	std::ofstream out("Plotting/density2d.csv");
	for (int i = 0; i < density.rows(); ++i)
	{
		for (int j = 0; j < density.cols(); ++j)
		{
			if (j != density.cols()-1){
				out << density(i,j) << ',';
			}
			else {out << density(i,j) << '\n';}
		}
	}*/ 

	return coeff;
}

void gaussianBarEvolution(){
	double barRadius{2}, timeStep{0.01};
	int nMax{24}, numbTimeSteps{2000};
	std::string kernelName{"Kernels/GaussianLog_2.out"};

	Bar2D bar(gaussianResolving(barRadius), 1/barRadius);
	VolterraSolver solver( kernelName, nMax, 2, numbTimeSteps, timeStep);

	solver.barRotation(bar, "Plotting/barCoefficents.csv", "Plotting/barEvolution.csv", false, false);
}



