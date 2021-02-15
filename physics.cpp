#include <iostream>
#include <typeinfo>
#include <Eigen/Dense>
#include <vector>
#include <complex>

#include "Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "DF_Class/Mestel.h"

#include "Volterra_Solver/VolterraSolver.h"
#include "Volterra_Solver/ExpansionCoeff.h"

#include "Bar2D/Bar2D.h"

#include "physics.h"


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

std::string kernelName(std::string dir, std::string stem, int Kka, int Rka, int m2)
{
	return dir + "/" + stem + "_" + std::to_string(Kka) + "_" + std::to_string(Rka) + "_" + std::to_string(m2) + ".out";
}

std::string kernelName(std::string dir, std::string stem, int Kka, double Rka, int m2) 
{
	return dir + "/" + stem + "_" + std::to_string(Kka) + "_" + std::to_string((int) round(Rka)) + "_" + std::to_string(m2) + ".out";
}


void kalnajsKernelsVaryingK() // Write in to use the new filename func
{
	Mestel DF;
	std::vector<double> Kka{4, 5, 6, 7};
	for (int i = 0; i < 4; ++i){
		
		std::cout << "Calculationg Kernels for Kka: " << Kka[i] << '\n';
		std::string file = "Kalnajs/Kalnajs_" + std::to_string((int) Kka[i]) + "_10";
		

		ActionAngleBasisContainer test(file, 10, 2, 5, 101, 20);
		VolterraSolver solver(10, 2, 2000, 0.025);
		std::string kernel = "Kernels/Kalnajs_" + std::to_string((int) Kka[i]) +"_10_2.out";
		solver.generateKernel(kernel, DF, test);
	}
}


void kalnajsKernelsVaryingR() // Write in to use the new filename func
{
	Mestel DF;
	std::vector<double> Rka{5, 10, 15, 20};
	for (int i = 0; i < 4; ++i){
		
		std::cout << "Calculationg kernel for Rka: " << Rka[i] << '\n';
		std::string file = "Kalnajs/Kalnajs_4_" + std::to_string((int) Rka[i]);
		

		ActionAngleBasisContainer test(file, 10, 2, 5, 101, 20);
		VolterraSolver solver(10, 2, 2000, 0.025);
		std::string kernel = "Kernels/Kalnajs_4_" + std::to_string((int) Rka[i]) +"_2.out";
		solver.generateKernel(kernel, DF, test);
	}
}

void kalnajsKernelsVaryingSigma(int l)
{
	std::vector<double> littleSigma{.25, .35, .50};
	
	std::string file = "Kalnajs/Kalnajs_4_20";
	ActionAngleBasisContainer test(file, 10, l, 5, 101, 20);
	
	VolterraSolver solver(10, l, 2000, 0.025);
	for (int i = 0; i < littleSigma.size(); ++i){
		
		std::cout << "Calculating kernel for littleSigma: " << littleSigma[i] << '\n';
		Mestel DF(1, 1, littleSigma[i]);
		
		std::string kernel = "Disk_Kicking/littleSigma_" + std::to_string((int) round(littleSigma[i]*100)) 
		+ "/Kalnajs"+ "_" + std::to_string(l) +".out";
		
		solver.generateKernel(kernel, DF, test);
	}
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


// Test Evolution //
// -------------- // 
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

// Paper 1 Functions //

std::string perturbationFilename(double littleSigma, double radius, int fourierHarmonic)
{
	return "Disk_Kicking/littleSigma_" + std::to_string((int) round(100*littleSigma)) + 
			"/Perturbation_" + std::to_string((int) round(10*radius)) + "_" +std::to_string(fourierHarmonic) + ".out";
}

template <class T>
void outPutPerturbation(const T & pd, double littleSigma, double radius, int numbTimeSteps = 2000)
{
	Eigen::VectorXcd t0Coeff(pd.maxRadialIndex()+1);
	for (int i = 0; i <= pd.maxRadialIndex(); ++i){
		t0Coeff[i] = pd.potential(radius, i);
	}
	// Do we need a step where we multiply by scriptE? 
	ExpansionCoeff holding(t0Coeff, numbTimeSteps, pd.maxRadialIndex());
	holding.write2File(perturbationFilename(littleSigma, radius, pd.fourierHarmonic()), 1);
}


void barKickingPerturbations()
{
	std::vector<double> radii {.1,.5,1, 2, 3, 5};

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> pd0(params, 10,0), pd1(params, 10, 1), pd2(params, 10, 2);

	for (auto i = radii.begin(); i != radii.end(); ++i){
		outPutPerturbation(pd0, .25, *i);
		outPutPerturbation(pd1, .25, *i);
		outPutPerturbation(pd2, .25, *i);
	}
}


//void SOMEFUNCTION THAT GENERATES ALL THE PERTURBATION FILES THAT WE NEED? 




// Paper 2 Functions // 

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


void barTesting() // I think that we can get rid of this. 
{
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,2);

	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	coef[0] = .1;
	Bar2D bar(coef, PD, 0);

}


// Misc // 
// ---- //

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
















