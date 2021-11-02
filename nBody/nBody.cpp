#include "nBodyFunctions.h"

#include <iostream>
#include <cmath>

double rem(double x, double y) {return x - y*floor(x/y);}

int main(){
	//barEvolution();
	
	//kalanajTest();
	//barEvolutionKalnajs();

	barEvolutionKalnajs(); 
	//individualOrbits();


	//std::cout << rem(102, 2*M_PI) << '\n';
	return 0;
} 