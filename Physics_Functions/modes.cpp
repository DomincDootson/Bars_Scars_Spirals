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

	// generatingEigenMode("Response_Matrix/RM/RM_Xi_8_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_8_Nu_6.out");
	// generatingEigenMode("Response_Matrix/RM/RM_Xi_7_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_7_Nu_6.out");
	// generatingEigenMode("Response_Matrix/RM/RM_Xi_6_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_6_Nu_6.out");
	// generatingEigenMode("Response_Matrix/RM/RM_Xi_5_Nu_6.out", "Plotting/Modes_Data/Xi_Modes/RM_Xi_5_Nu_6.out");

	//uniformSearch("Plotting/Modes_Data/Chi_Search", 4);
	//uniformSearch("Plotting/Modes_Data/Chi_Search", 8);

	// uniformSearchScarredDensity("10", "25", "-95_G", false);
	// uniformSearchScarredDensity("11", "25", "-95_G", false);
	// uniformSearchScarredDensity("12", "25", "-95_G", false);
	// uniformSearchScarredDensity("13", "25", "-95_G", false);
	// uniformSearchScarredDensity("14", "25", "-95_G", false);
	// uniformSearchScarredDensity("15", "25", "-95_G", false);
	// uniformSearchScarredDensity("16", "25", "-95_G", false);
	// uniformSearchScarredDensity("17", "25", "-95_G", false);
	// uniformSearchScarredDensity("18", "25", "-95_G", false);
	// uniformSearchScarredDensity("19", "25", "-95_G", false);
	// uniformSearchScarredDensity("20", "25", "-95_G", false);
	//mutlipleModes();

	// saveResponseMatrix("Response_Matrix/RM/Single_Scarred_RM/RM_12_25_-95.csv", "Kernels/Scarred_Kernels/AM_Scarred_Kernel_R_12_W_25_D_-95_G.out", 0.564407, 0.044407, 0.5 * (11.6873/12));
	// saveResponseMatrix("Response_Matrix/RM/Single_Scarred_RM/RM_14_25_-95.csv", "Kernels/Scarred_Kernels/AM_Scarred_Kernel_R_14_W_25_D_-95_G.out", 0.510169, 0.044916, 0.5 * (11.6873/12));
	// saveResponseMatrix("Response_Matrix/RM/Single_Scarred_RM/RM_16_25_-95.csv", "Kernels/Scarred_Kernels/AM_Scarred_Kernel_R_16_W_25_D_-95_G.out", 0.469492, 0.042373, 0.5 * (11.6873/12));
	// saveResponseMatrix("Response_Matrix/RM/Single_Scarred_RM/RM_18_25_-95.csv", "Kernels/Scarred_Kernels/AM_Scarred_Kernel_R_18_W_25_D_-95_G.out", 0.428814, 0.038305, 0.5 * (11.6873/12));
	// saveResponseMatrix("Response_Matrix/RM/Single_Scarred_RM/RM_20_25_-95.csv", "Kernels/Scarred_Kernels/AM_Scarred_Kernel_R_20_W_25_D_-95_G.out", 0.394915, 0.034746, 0.5 * (11.6873/12));


	eigenvector2Density("Response_Matrix/RM/Single_Scarred_RM/EV_20_25_-95.out", "Plotting/Modes_Data/Scarred_Density/SD_20_25_-95.csv");
	eigenvector2Density("Response_Matrix/RM/Single_Scarred_RM/EV_12_25_-95.out", "Plotting/Modes_Data/Scarred_Density/SD_12_25_-95.csv");
	eigenvector2Density("Response_Matrix/RM/Single_Scarred_RM/EV_14_25_-95.out", "Plotting/Modes_Data/Scarred_Density/SD_14_25_-95.csv");
	eigenvector2Density("Response_Matrix/RM/Single_Scarred_RM/EV_16_25_-95.out", "Plotting/Modes_Data/Scarred_Density/SD_16_25_-95.csv");
	eigenvector2Density("Response_Matrix/RM/Single_Scarred_RM/EV_18_25_-95.out", "Plotting/Modes_Data/Scarred_Density/SD_18_25_-95.csv");

	

	return 0; 

} //