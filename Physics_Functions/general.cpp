#include "generalFunctions.h"
#include <iostream>
#include <chrono>
int main()
{
	std::cout << "This is the general cpp\n";
	
	auto t1 = std::chrono::high_resolution_clock::now();
    generatingKalnajsKernels(1);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "f() took "
              << std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()
              << " milliseconds\n"; 

    testEvolutionKalanajs(1);
	return 0;
}