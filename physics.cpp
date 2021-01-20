#include <iostream>
#include <typeinfo>

#include "Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"

#include "DF_Class/Mestel.h"


#include "Volterra_Solver/VolterraSolver.h"

#include <Eigen/Dense>

#include <vector>
#include "physics.h"



void generatingBF()
{
	Mestel DF;

	std::vector<double> params{4, 20};
	PotentialDensityPairContainer<KalnajsBasis> PD(params, 10,1);

	ActionAngleBasisContainer test(10, 1, 5, 101, 20); 
	test.scriptW(PD, DF, "Kalnajs");
}



void generatingKernels()
{
	ActionAngleBasisContainer test("Kalnajs", 10, 1, 5, 101, 20);
	Mestel DF;


	VolterraSolver solver(10, 1, 5, 0.01);
	solver.generateKernel(DF, test);

	VolterraSolver readingIn("kernelFileName.csv", 10, 1, 5, 0.01);
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