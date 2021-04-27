#include "generalFunctions.h"
#include <iostream>
#include <chrono>
// Things to check,
// Is the rotation rate correct (are we nissing a factor of 2pi somewhere)
// Place the bar at different angles and see if that makes a difference 

int main()
{
	std::cout << "This is the general cpp\n";
	//generatingKalnajsBF(2);
  //generatingKalnajsKernels(2);
    /*auto t1 = std::chrono::high_resolution_clock::now();
    generatingKalnajsKernels(1);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "f() took "
              << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()
              << " milliseconds\n"; */

    //generatingKalnajsKernels(1);
  
  //generatingGaussianKernels(2);
  barTesting(0.1);
	return 0;
  // Can we put in some tests to make sure that the read in kernel that has the correct params
}

