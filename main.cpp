#include "physics.h"

int main()
{
	//barKickingPerturbations();
	//kalnajsKernelsVaryingSigma(2);
	//generatingKalnajsBF(1);
	kalnajsKernelsVaryingSigma(1);
	
	


	//kalnajsKernelsVaryingK();
	//kalnajsKernelsVaryingR();
	/*barVaryingAngularSpeed(); // Genereates data for different rotational speeds of bars
	barVaryingKka();
	barVaryingRka();*/
	
	
	return 0;
}


// Things to do,
//	1) If the kernel file doesn't exist it will run the error wrong number of timesteps, can we make it check if the kernel exists? 