#include "modesFunctions.h"

int main() {

	//generatingEigenMode("Response_Matrix/test.out", "Plotting/Modes_Data/JB_mode.csv");
	//unstableTimeEvolution();

	//generatingEigenMode("Response_Matrix/test.out", "Plotting/Modes_Data/JB_mode.csv", "Plotting/Modes_Data/JB_mode_Res.csv", "Plotting/Modes_Data/JB_mode_Time.csv");
	//testingKernelMethod();
	
	unstableSearch("Plotting/Modes_Data/Unstable_Mode_Search.csv");
	//uniformSeach("Plotting/Modes_Data/Unstable_Mode_Search_nu_4.csv");
	return 0; 

} 