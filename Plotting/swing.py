from math import *
from Density_Classes.OneDdensity import *

def axisymmetric_plots(files = ["Swing_Data/Axisymmetric/TP_sheet_axisymmetric_G_pert.csv", "Swing_Data/Axisymmetric/TP_Linear_axisymmetric_G_pert.csv", "Swing_Data/Axisymmetric/SC_sheet_axisymmetric_G_pert.csv", "Swing_Data/Axisymmetric/SC_Linear_axisymmetric_G_pert.csv"], indices = [5,20, 35, 50], centre = 8):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	densities = [OneDdensity(file) for file in files]

	fig, axs = plt.subplots(nrows = len(indices), ncols = 2, sharex = True, sharey = True)

	for i in range(len(indices)):

		axs[i, 0].plot(densities[0].radii-centre, densities[0][indices[i]], color = 'royalblue', label = "Sheet")
		axs[i, 0].plot(densities[1].radii-centre, (1/(2*pi))*densities[1].radii*densities[1][indices[i]], color = 'firebrick', label = "Disc")

		axs[i, 1].plot(densities[2].radii-centre, densities[2][indices[i]], color = 'royalblue')
		axs[i, 1].plot(densities[3].radii-centre, (1/(2*pi))*densities[3].radii*densities[3][indices[i]], color = 'firebrick')

		axs[i,1].set_ylabel(r"$t=$%.2f$\times\frac{2\times \pi}{\kappa}$" % (indices[i]*0.01), fontsize = 15)
		axs[i,1].yaxis.set_label_position("right")

	fig.text(0.01, 0.5, r"$\Sigma/\hat{\Sigma}_{e}$", va='center', rotation='vertical', fontsize = 20)	
	axs[-1,0].set_xlabel(r"$x/\Delta_{x}$", fontsize = 15)
	axs[-1,1].set_xlabel(r"$x/\Delta_{x}$", fontsize = 15)

	axs[0,0].set_title("Test Particle", fontsize = 20)
	axs[0,1].set_title("Self Consistent", fontsize = 20)
	axs[-1,0].legend(fontsize = 12)

	plt.show()


axisymmetric_plots()
