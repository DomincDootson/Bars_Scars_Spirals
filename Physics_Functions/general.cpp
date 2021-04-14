#include "generalFunctions.h"
#include <iostream>
#include <chrono>
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
   testEvolutionKalanajs(2);
	return 0;
  // Can we put in some tests to make sure that the read in kernel that has the correct params
}

