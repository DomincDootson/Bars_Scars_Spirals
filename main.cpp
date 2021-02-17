//#include "physics.h"

#include "Physics_Functions/barEvolutionFunctions.h"
#include "Physics_Functions/densityEvolutionFunctions.h"
#include "Physics_Functions/generalFunctions.h"

int main()
{
	//diskKickingPerturbations();
	diskKicking();
	
	
	
	


	//kalnajsKernelsVaryingK();
	//kalnajsKernelsVaryingR();
	/*barVaryingAngularSpeed(); // Genereates data for different rotational speeds of bars
	barVaryingKka();
	barVaryingRka();*/
	
	
	return 0;
}


// Things to do,
//	1) If the kernel file doesn't exist it will run the error wrong number of timesteps, can we make it check if the kernel exists? 