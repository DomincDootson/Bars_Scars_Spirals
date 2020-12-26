#include <iostream>
#include <typeinfo>

#include "Potential_Density_Pair_Classes/KalnajsBasis.h"
#include "Potential_Density_Pair_Classes/GaussianLogBasis.h"
#include "Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

#include "Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"

#include "DF_Class/Mestel.h"


#include "Volterra_Solver/VolterraSolver.h"


int main()
{


	Mestel DF;
	PotentialDensityPairContainer<KalnajsBasis> PD(10,0);
	ActionAngleBasisContainer test(10, 0, 5, 101, 20); 
	test.scriptW(PD, DF, "Kalnajs");




	// Things to do tomorrow
	//	 1) get the code to put the basis funcitons into the right file
	//	 2) Get the basis functions to put the code in the right dir
	//Gaussian.scriptE("gaussian1.out");
	return 0;
}