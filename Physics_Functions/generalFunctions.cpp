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


void gaussianScriptE(int m2)
{
	std::vector<double> params{24, .5, 15};
	PotentialDensityPairContainer<GaussianLogBasis> PD(params, 24,m2);
	PD.scriptE("gaussianScriptE.csv");
}

// Basis Function Generation //
// ------------------------- // 

// Put in a function to generate the names for the files that the BF are kept it


void generatingKalnajsBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,m2);

	ActionAngleBasisContainer test(10, m2, 5, 101, 20); 
	test.scriptW(PD, DF, "Kalnajs"); // Use file function name here
}

void generatingSpiralBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{24, .5, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params, 24,m2);

	ActionAngleBasisContainer test(24, m2, 10, 301, 20); 
	test.scriptW(PD, DF, "GaussianLog");
}

void kalnajBFVaryingK()
{
	Mestel DF;
	std::vector<double> Kka{4, 5, 6, 7};
	for (int i = 0; i < 4; ++i){
		
		std::cout << "Calculationg BF for Kka: " << Kka[i] << '\n';

		std::vector<double> params{Kka[i], 10};
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2);

		ActionAngleBasisContainer test(10, 2, 5, 101, 10);
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

		ActionAngleBasisContainer test(10, 2, 5, 101, 10);
		std::string file = "Kalnajs/Kalnajs_4_" + std::to_string((int) Rka[i]); // Use file function name here
		test.scriptW(PD, DF, file);
	}
}



// Kernel Generation //
// ----------------- // 

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

void spiralTestEvolution()
{
	VolterraSolver solver0("Kernels/GaussianLog_0.out", 24, 0, 2000, 0.01);
	solver0.activeFraction(.25);
	solver0.volterraSolver("GaussianLogTest0.csv", "someperturbation", false);	

	VolterraSolver solver1("Kernels/GaussianLog_1.out", 24, 1, 2000, 0.01);
	solver1.activeFraction(.25);
	solver1.volterraSolver("GaussianLogTest1.csv", "someperturbation", false);	

	VolterraSolver solver2("Kernels/GaussianLog_2.out", 24, 2, 2000, 0.01);
	solver2.activeFraction(.25);
	solver2.volterraSolver("GaussianLogTest2.csv", "someperturbation", false);	
}


std::string outComplex(std::complex<double> number)
{
	if (imag(number)>=0){
		return std::to_string(real(number)) + '+' + std::to_string(imag(number)) + 'j';
	}
	else{
		return std::to_string(real(number)) + std::to_string(imag(number)) + 'j';
	}
}

void savingArray(Eigen::ArrayXXcd array, std::string filename)
{
	std::ofstream out(filename);
	for (int i = 0 ; i< array.rows(); ++i)
	{
		for (int j = 0; j < array.cols(); ++j)
		{
			if (j == array.cols()-1){
				out << outComplex(array(i,j)) << '\n';
			}
			else {
				out << outComplex(array(i,j)) << ',';
			}
		}
	}

}

void savingPotentialArrays(std::string dir)
{
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(10+1);

	for (int i = 0; i<11; ++i)
	{
		coeff[i] = 1;
		Eigen::ArrayXXcd grid = PD.potentialArray(coeff, 401, 10);
		std::string filename = dir + "/KalnajsGrid_" + std::to_string(i) + ".csv";
		savingArray(grid, filename);
		coeff[i] = 0;
	}
}