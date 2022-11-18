#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <numeric>

#include "../Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Potential_Density_Pair_Classes/GaussianLogBasis.h"
//#include "../Potential_Density_Pair_Classes/SpiralBasis.h"


#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../Potential_Density_Pair_Classes/TheoreticalPDContainer.h"

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../DF_Class/Mestel.h"

#include "../Volterra_Solver/VolterraSolver.h"
#include "../Volterra_Solver/ExpansionCoeff.h" 

Eigen::VectorXcd guassianRing(double centre, double width) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");

	auto ic = [centre, width] (double radius) {return   1/sqrt(2 * M_PI*width*width)*exp(-0.5 *pow((radius-centre)/width , 2));}; 
	auto background_density = [] (double radius) {return (0.5/(2*M_PI *radius));};

	Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(49);
	double stepSize{0.01};
	for (int n = 0; n < coef.size(); ++n) {
		for (double r= stepSize; r < 15; r += stepSize) {coef(n) += - 2 * M_PI * pd.potential(r, n) * r * ic(r) * r;}
	}

	return stepSize*coef; 
	//for (int n = 0; n < coef.size(); ++n) {coef(n) = - 2 * M_PI*centre*pd.potential(centre, n);}
	//return coef; 
}


void ringEvolution(const std::string & filename, bool isSelfConsistent) {
	PotentialDensityPairContainer<KalnajsNBasis> pd("Potential_Density_Pair_Classes/Kalnajs_Numerical/KalnajsNumerical_0.dat");
	VolterraSolver solver("Kernels/Kalnajs_0.out", 48, 0, 100, 0.35543); 
	solver.activeFraction(0.5);

	solver.setInitalPerturbation(guassianRing(8, 1)); 	
	if (isSelfConsistent) {solver.deltaPerturbationConsistent();}
	else {solver.deltaPerturbationTest();}

	std::vector<double> radii(100);
	std::iota(radii.begin(), radii.end(), 0); 
	std::for_each(radii.begin(), radii.end(), [] (double & n){n =(0.02*n)+4;});

	solver.density1dEvolution(filename, pd); 
}