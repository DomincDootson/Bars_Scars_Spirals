#include "scarredEvolutionFunctions.h"

int main() {

	//circularCutThrought(); 
	//scarredPotential(); 
	//backgroundDF();

	//testingScaredMestelEvolution(); 
	//scarredModes();

	//circularInfall(); 
	//angularMomentumScar();

	//scarredDensity(0.366, "test.csv");
	//savingEvolutionKernel(1.0, 400, 0.5, "AM_Scarred_Kernel_R_10_W_25_D_10.out", false);
	// savingEvolutionKernel(1.0, 800, 0.5, "AM_Scarred_Kernel_R_10_W_25_D_-05_G.out", false);	
	// savingEvolutionKernel(1.1, 800, 0.5, "AM_Scarred_Kernel_R_11_W_25_D_-05_G.out", false);	
	// savingEvolutionKernel(1.2, 800, 0.5, "AM_Scarred_Kernel_R_12_W_25_D_-05_G.out", false);	
	// savingEvolutionKernel(1.3, 800, 0.5, "AM_Scarred_Kernel_R_13_W_25_D_-05_G.out", false);	
	// savingEvolutionKernel(1.4, 800, 0.5, "AM_Scarred_Kernel_R_14_W_25_D_-05_G.out", false);	
	// savingEvolutionKernel(1.5, 800, 0.5, "AM_Scarred_Kernel_R_15_W_25_D_-95_G.out", false);	
	// savingEvolutionKernel(1.6, 800, 0.5, "AM_Scarred_Kernel_R_16_W_25_D_-95_G.out", false);	
	// savingEvolutionKernel(1.7, 800, 0.5, "AM_Scarred_Kernel_R_17_W_25_D_-95_G.out", false);	
	// savingEvolutionKernel(1.8, 800, 0.5, "AM_Scarred_Kernel_R_18_W_25_D_-95_G.out", false);	
	// savingEvolutionKernel(1.9, 800, 0.5, "AM_Scarred_Kernel_R_19_W_25_D_-95_G.out", false);	
	//savingEvolutionKernel(0.15, 800, 0.95, "AM_Scarred_Kernel_R_20_W_15_D_-95_G_inner.out", false);	
	//savingEvolutionKernel(0.20, 800, 0.50, "AM_Scarred_Kernel_R_20_W_20_D_-50_G_inner.out", false);	

	

	savingDensityFile(0, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-0_G.csv");
	savingDensityFile(0.1, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-10_G.csv");
	savingDensityFile(0.2, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-20_G.csv");
	savingDensityFile(0.3, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-30_G.csv");
	savingDensityFile(0.4, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-40_G.csv");
	savingDensityFile(0.5, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-50_G.csv");
	savingDensityFile(0.6, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-60_G.csv");
	savingDensityFile(0.7, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-70_G.csv");
	savingDensityFile(0.8, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-80_G.csv");
	savingDensityFile(0.9, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-90_G.csv");
	savingDensityFile(1.0, "WKB_Disc/Disc_Density/Tapered_R_10_W_25_D_-100_G.csv");

	return 0;
}