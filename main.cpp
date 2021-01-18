/*#include <iostream>
#include <typeinfo>

#include "Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"

#include "DF_Class/Mestel.h"


#include "Volterra_Solver/VolterraSolver.h"

#include <Eigen/Dense>

#include <vector>*/ 

#include "physics.h"
	

int main()
{


	//Mestel DF;

	kalnajBasisFunctionsVaryingK();




	//ActionAngleBasisContainer test(10, 0, 5, 101, 20); 
	//test.scriptW(PD, DF, "Kalnajs");
	/*Eigen::VectorXcd coef(25);

	for (int i = 0; i<25; ++i){
		coef(i) = i;
	}
	

	Eigen::MatrixXcd potential = PD.potentialGrid(coef, 801, 20);
	Eigen::MatrixXcd density = PD.densityGrid(coef, 801, 20);

	std::cout << PD.potentialResolving(potential, 20) << '\n' << '\n'; 
	std::cout << PD.densityResolving(density, 20) << '\n'; */ 





	// Things to do tomorrow
	//	 1) get the code to put the basis funcitons into the right file
	//	 2) Get the basis functions to put the code in the right dir

	return 0;
}