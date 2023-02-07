from mpl_toolkits.axes_grid1 import make_axes_locatable

from math import *
from Density_Classes.OneDdensity import *
from Density_Classes.TwoDdensity import *

def axisymmetric_comparison_plots(files = ["Swing_Data/Axisymmetric/TP_sheet_axisymmetric_G_pert.csv", "Swing_Data/Axisymmetric/TP_Linear_axisymmetric_G_pert.csv", "Swing_Data/Axisymmetric/SC_sheet_axisymmetric_G_pert.csv", "Swing_Data/Axisymmetric/SC_Linear_axisymmetric_G_pert.csv"], indices = [5,20, 35, 50], centre = 8, delta = 1):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	densities = [OneDdensity(file) for file in files]

	fig, axs = plt.subplots(nrows = len(indices), ncols = 2, sharex = True, sharey = True)

	for i in range(len(indices)):

		axs[i, 0].plot((densities[0].radii-centre)/delta, densities[0][indices[i]], color = 'royalblue', label = "Sheet")
		axs[i, 0].plot((densities[1].radii-centre)/delta, densities[1][indices[i]], color = 'firebrick', label = "Disc")

		axs[i, 1].plot((densities[2].radii-centre)/delta, densities[2][indices[i]], color = 'royalblue')
		axs[i, 1].plot((densities[3].radii-centre)/delta, densities[3][indices[i]], color = 'firebrick')

		axs[i,1].set_ylabel(r"$t=$%.2f$\times\frac{2\times \pi}{\kappa}$" % (indices[i]*0.01), fontsize = 15)
		axs[i,1].yaxis.set_label_position("right")

	fig.text(0.01, 0.5, r"$\Sigma/\hat{\Sigma}_{e}$", va='center', rotation='vertical', fontsize = 20)	
	axs[-1,0].set_xlabel(r"$x/\Delta_{x}$", fontsize = 15)
	axs[-1,1].set_xlabel(r"$x/\Delta_{x}$", fontsize = 15)

	axs[0,0].set_title("Test Particle", fontsize = 20)
	axs[0,1].set_title("Self Consistent", fontsize = 20)
	axs[-1,0].legend(fontsize = 12)

	plt.show()

def amplificationPlots(files = ["Swing_Data/Amplification_Plots/Binney_Amplification.csv", "Swing_Data/Amplification_Plots/Toomre_Amplification.csv"]):
	amplification = [OneDdensity(file) for file in files]

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(ncols = 2, sharex = True)
	axs[0].set_yscale('log')
	axs[1].set_yscale('log')

	colors, Q = ["navy","cornflowerblue","lightcoral","firebrick"], [1.2, 1.3, 1.5, 2.0]

	axs[0].plot(amplification[0].radii, amplification[0][0], color = colors[0], label = str(Q[0]))
	axs[1].plot(amplification[1].radii, amplification[1][0], color = colors[0])

	axs[0].plot(amplification[0].radii, amplification[0][1], color = colors[1], label = str(Q[1]))
	axs[1].plot(amplification[1].radii, amplification[1][1], color = colors[1])

	axs[0].plot(amplification[0].radii, amplification[0][2], color = colors[2], label = str(Q[2]))
	axs[1].plot(amplification[1].radii, amplification[1][2], color = colors[2])

	axs[0].plot(amplification[0].radii, amplification[0][3], color = colors[3], label = str(Q[3]))
	axs[1].plot(amplification[1].radii, amplification[1][3], color = colors[3])


	axs[0].legend(title = "Q", fontsize = 15, title_fontsize = 15)

	axs[0].set_xlabel(r"$\lambda/\lambda_{\mathrm{crit}}$", fontsize = 15)
	axs[1].set_xlabel(r"$\lambda/\lambda_{\mathrm{crit}}$", fontsize = 15)

	axs[0].set_ylabel(r"$A_{\mathrm{max}}$", fontsize = 15)
	axs[1].set_ylabel(r"$T_{\mathrm{max}}$", fontsize = 15)



	plt.show()

def max_density_Response():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	CC =  readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentClockwise.csv")
	
	sheet_cold, sheet_warm = readingInRealCSV("swing_test_12.csv"), readingInRealCSV("swing_test_15.csv")
	plt.scatter((4/CC[1:7,0]), CC[1:7,1], color = 'firebrick')

	
	
	plt.plot(4/sheet_cold[:,0], sheet_cold[:,1], color = 'royalblue', label = "1.3")
	plt.plot(4/sheet_warm[:,0], sheet_warm[:,1], color = 'firebrick', label = "1.5")
	


	plt.legend(title = "Q", fontsize = 15, title_fontsize  =15)
	plt.ylabel(r"$\rho_{Max}/M$", fontsize = 15)
	plt.xlabel(r"$\lambda/\lambda_{crit}$", fontsize = 15)

	plt.yscale('log')

	plt.show()

def selfgrav_Amplification():
	CC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentClockwise.csv")
	TC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/TestClockwise.csv")
	amp = OneDdensity("Swing_Data/Amplification_Plots/Binney_Amplification.csv")
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	
	plt.plot(amp.radii, amp[1], color = 'royalblue', label = 1.3)
	
	plt.scatter(4/CC[1:6,0], CC[1:6,1]/TC[1:6,1], color = 'firebrick')
	
	plt.scatter(4/CC[10::2,0], CC[10::2, 1]/TC[10::2,1], color = 'firebrick')
	plt.plot(amp.radii, amp[2], color = 'firebrick', label = "1.5")


	
	

	plt.legend(title_fontsize = 15, title = "Q")
	plt.ylabel(r"$A_{\mathrm{max}}$", fontsize = 15)
	plt.xlabel(r"$\lambda/\lambda_{crit}$", fontsize = 15)

	 
	plt.yscale('log')
	plt.show()

def sum_over_l(filestem, max_l, nrows):
	densities = [TwoDdensity(filestem + str(l) +".csv") for l in range(max_l)]
	arrays = []
	for i in range(1, nrows+1):
		arrays.append(np.zeros_like(densities[0][0]))

		for den in densities:
			arrays[-1] += den[i]

	return arrays

def max_values(densities):
	
	max_vals = [np.amax(den) for den in densities]
	return max(max_vals)

def min_values(densities):
	min_vals = [np.amin(den) for den in densities]
	return min(min_vals)

def swing_plots(filestem = ["Swing_Data/Amplification_Fixed_5/CC_Evolution_", "Swing_Data/Amplification_Fixed_5/TC_Evolution_", "Swing_Data/Amplification_Fixed_5/CA_Evolution_", "Swing_Data/Amplification_Fixed_5/TA_Evolution_"] , max_l = 10, nrows = 4):
	densities = [sum_over_l(file, max_l, nrows) for file in filestem] 
	max_vals, min_vals = [max_values(dens) for dens in densities], [min_values(dens) for dens in densities]
	colorbar_scale = [max(abs(mi), abs(mx)) for mi, mx, in zip(min_vals, max_vals)]
	print(colorbar_scale)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

		
	
	fig, axs = plt.subplots(ncols = len(filestem), nrows = nrows, sharex = True, sharey = True)
	
	# Can we set them to be square 

	x, y  = np.linspace(-8, 8, np.shape(densities[0][0])[0]), np.linspace(-8, 8, np.shape(densities[0][0])[1])
	XX, YY = np.meshgrid(x, y)
	xCir, yCir = [5 * cos(t) for t in np.linspace(0, 2*pi)], [5 * sin(t) for t in np.linspace(0, 2*pi)]

	for i in range(nrows):
		for j in range(len(filestem)):
			max_d = np.amax((densities[j])[i])
			im = axs[i,j].contourf(XX, YY, (densities[j])[i], levels = np.linspace(-colorbar_scale[j], colorbar_scale[j]))
			axs[i,j].set_aspect('equal', 'box')
			axs[i,j].plot(xCir, yCir, color = 'firebrick', linestyle = '--', alpha = 0.6)
			# if i == nrows-1:
				
			# 	fig.colorbar(im, ax = axs[i,j], orientation = 'horizontal', use_gridspec = True)
			#axs[i,j].contour(XX, YY, (densities[j])[i], levels = [0.25*max_d, 0.5*max_d, 0.75*max_d], colors = "firebrick")

	axs[0,0].set_title("Self Consistent", fontsize = 15)
	axs[0,2].set_title("Self Consistent", fontsize = 15)
	axs[0,1].set_title("Test Particle", fontsize = 15)
	axs[0,3].set_title("Test Particle", fontsize = 15)

	
	axs[0,-1].set_ylabel(r"t=1.66", fontsize = 15)
	axs[0,-1].yaxis.set_label_position("right")
	axs[1,-1].set_ylabel(r"t=3.33", fontsize = 15)
	axs[1,-1].yaxis.set_label_position("right")
	axs[2,-1].set_ylabel(r"t=5.00", fontsize = 15)
	axs[2,-1].yaxis.set_label_position("right")
	axs[3,-1].set_ylabel(r"t=6.66", fontsize = 15)
	axs[3,-1].yaxis.set_label_position("right")

	xCir, yCir = [2.5 * cos(t) for t in np.linspace(0, 2*pi)], [2.5 * sin(t) for t in np.linspace(0, 2*pi)]
	axs[-1,0].plot(xCir, yCir, color = 'firebrick', linestyle = '--', alpha = 0.6)
	axs[-1,2].plot(xCir, yCir, color = 'firebrick', linestyle = '--', alpha = 0.6)

	plt.figtext(0.31,0.93, "Leading", va="center", ha="center", fontsize=17)
	plt.figtext(0.72,0.93,"Trailing", va="center", ha="center", fontsize=17)

	plt.show()

	# Create the plots 
	


#axisymmetric_comparison_plots()
#axisymmetric_comparison_plots(files = ["Swing_Data/NonAxisymmetric_Comparisons/TP_sheet_Nonaxisymmetric_G_pert.csv", "Swing_Data/NonAxisymmetric_Comparisons/TP_Linear_Nonaxisymmetric_G_pert.csv", "Swing_Data/NonAxisymmetric_Comparisons/SC_sheet_Nonaxisymmetric_G_pert.csv", "Swing_Data/NonAxisymmetric_Comparisons/SC_Linear_Nonaxisymmetric_G_pert.csv"], indices = [10,20, 30, 40], centre = 8, delta = 0.5)

#amplificationPlots()

swing_plots()

#selfgrav_Amplification()

#max_density_Response()

# den = TwoDdensity("Swing_Data/Amplification_Fixed_5/CC_Evolution_0.csv")
# plt.imshow(den[1]- den[2])
# plt.show()