#include <iostream>
#include <typeinfo>
#include <Eigen/Dense>
#include <vector>

#include "Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "DF_Class/Mestel.h"
#include "Volterra_Solver/VolterraSolver.h"
#include "Bar2D/Bar2D.h"
#include "physics.h"


void gaussianScriptE(int m2)
{
	std::vector<double> params{24, .5, 15};
	PotentialDensityPairContainer<GaussianLogBasis> PD(params, 24,m2);
	PD.scriptE("gaussianScriptE.csv");
}


void generatingKalnajsBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,m2);

	ActionAngleBasisContainer test(10, m2, 5, 101, 20); 
	test.scriptW(PD, DF, "Kalnajs");
}


void generatingSpiralBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{24, .5, 15};
	PotentialDensityPairContainer<GaussianLogBasis> PD(params, 24,m2);

	ActionAngleBasisContainer test(24, m2, 10, 301, 20); 
	test.scriptW(PD, DF, "GaussianLog");
}



void generatingKalnajsKernels(int m2)
{
	ActionAngleBasisContainer test("Kalnajs", 10, m2, 5, 101, 20);
	Mestel DF;


	VolterraSolver solver(10, m2, 2000, 0.01);

	std::string kernel = "Kernels/Kalnajs" +std::to_string(m2) +".out";
	solver.generateKernel(kernel, DF, test);
}

void generatingGaussianKernels(int m2)
{
	ActionAngleBasisContainer test("GaussianLog", 24, m2, 10, 301, 20);
	Mestel DF;


	VolterraSolver solver(24, m2, 2000, 0.01);

	std::string kernel = "Kernels/GaussianLog" +std::to_string(m2) +".out";
	solver.generateKernel(kernel, DF, test);

	VolterraSolver readingIn(kernel, 24, m2, 2000, 0.01);
}




void kalnajsPerturbation()
{
	std::ofstream out("someperturbation");
	out << 10 << '\n';
	for (int i =0; i<2000; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			out << 0.01*sin(M_PI * (i/2000.0)) << " " << 0 << " ";
		}
		for (int j = 3; j<10; ++j)
		{
			out << 0 << " " << 0 << " ";
		}
		out << 0 << " " << 0 << '\n'; 
	}
	out.close();
}


void gaussianPerturbation()
{
	std::ofstream out("someperturbation");
	out << 24 << '\n';
	for (int i =0; i<2000; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			out << 0 << " " << 0 << " ";			
		}
		for (int j = 10; j<13; ++j)
		{
			out << 0.01*sin(M_PI * (i/2000.0)) << " " << 0 << " ";
		}
		for (int j = 13; j < 24; ++j)
		{
			out << 0 << " " << 0 << " ";
		}
		out << 0 << " " << 0 << '\n'; 
	}
	out.close();
}


void testEvolutionKalanajs(int m2)
{
	kalnajsPerturbation();
	VolterraSolver solver("kernelFileName.csv", 10, m2, 2000, 0.01);
	solver.activeFraction(.25);	
	solver.volterraSolver("evolution.csv", "someperturbation", true);
	solver.resetActiveFraction();
}


void testEvolutionGaussian(int m2)
{
	gaussianPerturbation();
	std::string kernel = "Kernels/GaussianLog" +std::to_string(m2) +".out";
	VolterraSolver solver(kernel, 24, m2, 2000, 0.01);
	solver.activeFraction(1/9.5);
	solver.volterraSolver("evolution.csv", "someperturbation", false);
	solver.resetActiveFraction();
}



void kernelFlipped(){
	VolterraSolver solver0("Kernels/Kalnajs_0.out", 10, 0, 2000, 0.01);
	solver0.activeFraction(.25);
	solver0.kernelWrite2fileFlipped("KalnajsBlock0.out");

	VolterraSolver solver1("Kernels/Kalnajs_1.out", 10, 1, 2000, 0.01);
	solver1.activeFraction(.25);
	solver1.kernelWrite2fileFlipped("KalnajsBlock1.out");
	
	VolterraSolver solver2("Kernels/Kalnajs_2.out", 10,2 , 2000, 0.01);
	solver2.activeFraction(.25);
	solver2.kernelWrite2fileFlipped("KalnajsBlock2.out");

}



// Bar Rotating Physics

void kalnajBasisFunctionsVaryingK() // Outputs the values of the density basis functions for plotting
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

// Caclculation of the MOI. 

void barTesting()
{
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,2);

	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	coef[0] = 1;
	Bar2D bar(coef, PD, 0);

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

	VolterraSolver solver("Kernels/Kalnajs_2.out", 10, 2, 2000, 0.01);
	solver.activeFraction(.25);	
	


	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,2);
	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	coef[0] = .01; 
	
	for (int i =0; i < angSpeed.size(); ++i){
		// Construct the bar
		Bar2D bar(coef, PD, angSpeed[i]);
		std::string evolutionFilename = evolutionFileName("BarSlowing", angSpeed[i]);
		std::string outFilename = coeffFileName("BarSlowing", angSpeed[i]);
		solver.barRotation(bar, outFilename, evolutionFilename, true, true); 
	}
}




