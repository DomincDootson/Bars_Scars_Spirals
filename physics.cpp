#include <iostream>
#include <typeinfo>

#include "Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"

#include "DF_Class/Mestel.h"


#include "Volterra_Solver/VolterraSolver.h"

#include "Bar2D/Bar2D.h"

#include <Eigen/Dense>

#include <vector>
#include "physics.h"



void generatingBF(int m2)
{
	Mestel DF;
	
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,m2);

	ActionAngleBasisContainer test(10, m2, 5, 101, 20); 
	test.scriptW(PD, DF, "Kalnajs");
}



void generatingKernels(int m2)
{
	ActionAngleBasisContainer test("Kalnajs", 10, m2, 5, 101, 20);
	Mestel DF;


	VolterraSolver solver(10, m2, 2000, 0.01);
	solver.generateKernel("kernelFileName.csv", DF, test);

	//VolterraSolver readingIn("kernelFileName.csv", 10, 1, 5, 0.01);
}

void somePerturbation()
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

void testEvolution()
{
	somePerturbation();
	VolterraSolver solver("kernelFileName.csv", 10, 0, 2000, 0.01);
	solver.activeFraction(.25);	
	solver.volterraSolver("evolution.csv", "someperturbation", true);
	solver.resetActiveFraction();
}

void barTesting()
{
	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,2);

	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(11);
	Bar2D bar(coef, PD, 0);

}


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
		PotentialDensityPairContainer<KalnajsBasis> PD(params, 10, 1);
		for (int i = 0; i<999; ++i)
		{
			out << PD.density(i*0.01,0) << ',';
		}
		out << PD.density(999*0.01,0)  << '\n';

	}
	out.close();
}