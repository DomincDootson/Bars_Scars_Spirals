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

#include "../DF_Function/DFfunction.h"


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

	ActionAngleBasisContainer test("Kalnajs", 10, m2, 5, 101, 20); 
	test.scriptW(PD, DF, "Kalnajs/Kalnajs_4_20"); // Use file function name here
}

void generatingSpiralBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{24, .5, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> PD(params, 24,m2);

	ActionAngleBasisContainer test("GaussianLog", 24, m2, 10, 301, 20); 
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



// Kernel Generation //
// ----------------- // 

void generatingKalnajsKernels(int m2)
{
	ActionAngleBasisContainer test("Kalnajs/Kalnajs_4_20", "Kalnajs", 10, m2, 5, 101, 20);
	Mestel DF;


	/*VolterraSolver solver(10, m2, 2000, 0.025);

	std::string kernel = "test2000.out";//"Kernels/Kalnajs" +std::to_string(m2) +".out";
	solver.generateKernel(kernel, DF, test);

	VolterraSolver solver1(10, m2, 1000, 0.05);

	std::string kernel1 = "test1000.out";//"Kernels/Kalnajs" +std::to_string(m2) +".out";
	solver1.generateKernel(kernel1, DF, test);*/

	VolterraSolver solver2(10, m2, 200, 0.1);

	std::string kernel2 = "test200.out";//"Kernels/Kalnajs" +std::to_string(m2) +".out";
	solver2.generateKernel(kernel2, DF, test);
}

void generatingGaussianKernels(int m2)
{
	ActionAngleBasisContainer test("GaussianLog", "GaussianLog", 24, m2, 10, 301, 20);
	Mestel DF;


	VolterraSolver solver(24, m2, 200, 0.1);

	//std::string kernel = "Kernels/GaussianLog" +std::to_string(m2) +".out";
	std::string kernel = "test200.out";
	solver.generateKernel(kernel, DF, test);

	//VolterraSolver readingIn(kernel, 24, m2, 2000, 0.01);
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

// Testing the evolution //
// --------------------- // 

std::complex<double> phase(double time){
	std::complex<double> unitComplex(0,1);
	return exp(-2 * (0.25*M_PI) *  unitComplex);//exp(-2 * time * 0 *  unitComplex); // could it be this function? 
}

void kalnajsPerturbation(int nStep)
{
	std::ofstream out("someperturbation"); //someperturbation
	out << 10 << '\n';
	for (int i =0; i<100; ++i)
	{
		double size{0.01*sin(M_PI * (i/(double) 200))}, time{i*0.1};
		out << size * phase(time).real() << " " << size * phase(time).imag() << " ";
		for (int j = 1; j<10; ++j)
		{
			out << 0 << " " << 0 << " ";
		}
		out << 0 << " " << 0 << '\n'; 
	}
	for (int i =100; i<200; ++i)
	{
		
		double size{0.01}, time{i*0.1};
		out << size * phase(time).real() << " " << size * phase(time).imag() << " ";
		
		for (int j = 1; j<10; ++j)
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
	kalnajsPerturbation(200);
	
	VolterraSolver solver("test200.out", 10, m2, 200, 0.1);
	solver.activeFraction(.25);	
	solver.coefficentEvolution("evolution200.csv", 
		"someperturbation", false);
	solver.resetActiveFraction();
}


void testEvolutionGaussian(int m2)
{
	gaussianPerturbation();
	std::string kernel = "Kernels/GaussianLog" +std::to_string(m2) +".out";
	VolterraSolver solver(kernel, 24, m2, 2000, 0.01);
	solver.activeFraction(1/9.5);
	solver.coefficentEvolution("evolution.csv", "someperturbation", false);
	solver.resetActiveFraction();
}

void spiralTestEvolution()
{
	VolterraSolver solver0("Kernels/GaussianLog_0.out", 24, 0, 2000, 0.01);
	solver0.activeFraction(.25);
	solver0.coefficentEvolution("GaussianLogTest0.csv", "someperturbation", false);	

	VolterraSolver solver1("Kernels/GaussianLog_1.out", 24, 1, 2000, 0.01);
	solver1.activeFraction(.25);
	solver1.coefficentEvolution("GaussianLogTest1.csv", "someperturbation", false);	

	VolterraSolver solver2("Kernels/GaussianLog_2.out", 24, 2, 2000, 0.01);
	solver2.activeFraction(.25);
	solver2.coefficentEvolution("GaussianLogTest2.csv", "someperturbation", false);	
}

// Testing Bar Testing //
// ------------------- //
template <class Tbf>
Eigen::VectorXcd gaussianBar(const Tbf & pd) {
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(pd.maxRadialIndex() + 1);
	for (int i = 0; i < pd.maxRadialIndex() + 1; ++i){
		coeff(i) = 0.01*pd.potential(2.06271, i);
	}
	return - (pd.scriptE()).inverse() * coeff;
}


void barSize(int nStep){
	std::ofstream out("Bar2D/barSize.out"); out << nStep << '\n';
	for (int i =0; i < nStep; ++i) {
		out << i * 0.1 << " " << sin(0.5*M_PI*(i/((double) nStep))) << '\n'; 
	}
	out.close();
}

void barTesting(double omega0){
	barSize(100);
	std::vector<double> params{24, .5, 15};
 	PotentialDensityPairContainer<GaussianLogBasis> pd(params, 24, 2);

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(24+1); coeff =  gaussianBar(pd);
	Bar2D bar(coeff, omega0, "Bar2D/barSize.out");
 

	VolterraSolver solver("test200.out", 24, 2, 200, 0.1);
	solver.activeFraction(.25);
	solver.barRotation(bar, "barCoeff.csv", "barEvolution.csv", false, false, true);
}

void gridTesting(double angD){
	std::complex<double> unitComplex(0,1);
	double ang = (angD/360.0) * 2*M_PI;

	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(11); 
	coeff(0) = .01*exp(-2*ang*unitComplex);

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 2);

	Eigen::ArrayXXd array = PD.potentialArrayReal(coeff, 401, 20);
	std::ofstream out("nBody/pertGrid.csv");
	for (int i =0; i < array.rows(); ++i){
		for (int j = 0; j < array.cols(); ++j){
			if (j == array.cols()-1) {out << array(i,j) <<'\n';}
			else {out << array(i,j) << ',';}
		}
	}
	out.close();

}
