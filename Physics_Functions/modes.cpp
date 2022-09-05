#include <iostream>
#include <vector>
#include <string>
#include "modesFunctions.h"


void mutlipleModes() {
	std::vector<int> radii {10, 12, 14, 16, 18, 20, 22, 24, 26, 28};

	for (auto r : radii) {
		uniformSearchScarredDensity(std::to_string(r), "25", "10");
	}
}

int main() {

	//generatingEigenMode("Response_Matrix/test.out", "Plotting/Modes_Data/JB_mode.csv");
	//unstableTimeEvolution();

	//generatingEigenMode("Response_Matrix/test.out", "Plotting/Modes_Data/JB_mode.csv", "Plotting/Modes_Data/JB_mode_Res.csv", "Plotting/Modes_Data/JB_mode_Time.csv");
	//testingKernelMethod();
	//uniformSearchKernel("Plotting/Modes_Data/Unstable_Mode_Search_kernel_nu_8.csv", 4, 0.23);
	//uniformSearchMatrix("Plotting/Modes_Data/Unstable_Mode_Search_nu_4.csv", 4);
	
	// saveResponseMatrix("Response_Matrix/RM/RM_Xi_8_Nu_6.out", 0.827966, 0.126923, 0.8);
	// saveResponseMatrix("Response_Matrix/RM/RM_Xi_7_Nu_6.out", 0.798136, 0.080769, 0.7);
	// saveResponseMatrix("Response_Matrix/RM/RM_Xi_6_Nu_6.out", 0.765932, 0.034615, 0.6);
	// saveResponseMatrix("Response_Matrix/RM/RM_Xi_5_Nu_6.out", 0.74661, -0.001231, 0.5);

	generatingEigenMode("Response_Matrix/RM/RM_Xi_8_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_8_Nu_6.out");
	generatingEigenMode("Response_Matrix/RM/RM_Xi_7_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_7_Nu_6.out");
	generatingEigenMode("Response_Matrix/RM/RM_Xi_6_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_6_Nu_6.out");
	generatingEigenMode("Response_Matrix/RM/RM_Xi_5_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_5_Nu_6.out");

	//uniformSearch("Plotting/Modes_Data/Chi_Search", 4);
	//uniformSearch("Plotting/Modes_Data/Chi_Search", 8);

	//uniformSearchScarredDensity("12", "25", "-05", true);
	//mutlipleModes();

	

	

	return 0; 

} //