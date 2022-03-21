#include "nBodyFunctions.h"

#include <iostream>
#include <cmath>

int main(){
	//barEvolutionKalnajs("Evolution_Test_Cold", false, 0.35); 
	//barEvolutionKalnajs("Evolution_Test_Warm", false, 0.45);

	orbitSection(); 
	/* 
	barEvolutionKalnajs("Evolution_Self_Cold", true, 0.35); 
	barEvolutionKalnajs("Evolution_Self_Warm", true, 0.45); */
	
	//spiralTesting();


	return 0;
} 

// Things that need to be put right
// 1) Make the code run the sampling script
// 2) Outputting coef and not calculating jacobi