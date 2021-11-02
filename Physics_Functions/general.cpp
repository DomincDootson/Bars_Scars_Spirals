#include "generalFunctions.h"
#include <iostream>
#include <chrono>

#include <vector>
// Things to check,
// Is the rotation rate correct (are we nissing a factor of 2pi somewhere)
// Place the bar at different angles and see if that makes a difference 

int main()
{
	//generatingSpiralBF(2);
	//generatingSpiralBFDiffTemp(2);
	//generatingKalnajsKernels(2);
	//generatingKalnajsBF(2);
	//testingBarTorque();
	testingFitting();
	return 0;
  // Can we put in some tests to make sure that the read in kernel that has the correct params
}

