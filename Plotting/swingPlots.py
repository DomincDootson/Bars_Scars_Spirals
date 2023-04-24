from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
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
	fig, axs = plt.subplots()
	CC =  readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentClockwise.csv")
	CC_C =  readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentClockwise_Cold.csv")
	
	sheet_cold, sheet_warm = readingInRealCSV("Swing_Data/Amplification_Fixed_5/swing_test_13.csv"), readingInRealCSV("Swing_Data/Amplification_Fixed_5/swing_test_15.csv")
	n = 10
	axs.scatter((4/CC[:n,0]), CC[:n,1], color = 'firebrick')
	axs.scatter((4/CC_C[:n,0]), CC_C[:n,1], color = 'royalblue')

	
	
	axs.plot(4/sheet_cold[:,0], sheet_cold[:,1], color = 'royalblue', label = "1.3")
	axs.plot(4/sheet_warm[:,0], sheet_warm[:,1], color = 'firebrick', label = "1.5")
	
	ax2 = axs.secondary_xaxis('top', functions=(lambda x: 4 / x, lambda x: 4 / x))
	ax2.set_xticks(range(15))
	ax2.set_xlabel(r"$\ell$", fontsize =15)
	labels = [str(i) for i in range(9)] + ['' for _ in range(9, 15)]
	ax2.set_xticklabels(labels)


	axs.legend(title = "Q", fontsize = 15, title_fontsize  =15)
	axs.set_ylabel(r"$\rho_{Max}/M$", fontsize = 15)
	axs.set_xlabel(r"$\lambda/\lambda_{crit}$", fontsize = 15)

	axs.set_yscale('log')

	plt.show()

def selfgrav_Amplification():
	CC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentClockwise_Cold.csv")
	TC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/TestClockwise_Cold.csv")
	amp = OneDdensity("Swing_Data/Amplification_Plots/Binney_Amplification.csv")
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots()
	axs.set_xlim([0.05, 5])
	
	axs.plot(amp.radii, amp[1], color = 'royalblue', label = 1.3)
	r = [1.2*CC[0,1]/TC[0,1], 28.04, 25.11] + list(1.1*CC[3:,1]/TC[3:,1])
	print(len(r), len([4/x for x in range(1, len(r)-1)]))
	axs.scatter([4/x for x in range(1, len(r)+1)], r, color = 'royalblue')
	
	

	CC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentClockwise.csv")
	TC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/TestClockwise.csv")
	axs.scatter(4/CC[1:,0], CC[1:,1]/TC[1:,1], color = 'firebrick')
	axs.plot(amp.radii, amp[2], color = 'firebrick', label = "1.5")
	
	

	ax2 = axs.secondary_xaxis('top', functions=(lambda x: 4 / x, lambda x: 4 / x))
	ax2.set_xticks(range(15))
	ax2.set_xlabel(r"$\ell$", fontsize =15)
	labels = [str(i) for i in range(9)] + ['' for _ in range(9, 15)]
	ax2.set_xticklabels(labels)
	




	
	

	axs.legend(title_fontsize = 15, title = "Q")
	axs.set_ylabel(r"$A_{\mathrm{max}}$", fontsize = 15)
	axs.set_xlabel(r"$\lambda/\lambda_{crit}$", fontsize = 15)

	 
	axs.set_yscale('log')
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

		
	
	fig, axs = plt.subplots(ncols = len(filestem), nrows = nrows+1,  gridspec_kw={"height_ratios":([1 for _ in range(len(filestem))] + [0.05])})
	
	# Can we set them to be square 

	x, y  = np.linspace(-8, 8, np.shape(densities[0][0])[0]), np.linspace(-8, 8, np.shape(densities[0][0])[1])
	XX, YY = np.meshgrid(x, y)
	xCir, yCir = [5 * cos(t) for t in np.linspace(0, 2*pi)], [5 * sin(t) for t in np.linspace(0, 2*pi)]

	for i in range(nrows):
		for j in range(len(filestem)):
			max_d = np.amax((densities[j])[i])
			im = axs[i,j].contourf(XX, YY, (densities[j])[i], levels = np.linspace(-colorbar_scale[j], colorbar_scale[j]))
			axs[i,j].set_aspect('equal', 'box')
			axs[i,j].plot(xCir, yCir, color = 'firebrick', linestyle = '--', alpha = 0.75)
			
			
			if i == nrows-1:
				cb = fig.colorbar(im, cax=axs[-1, j], orientation="horizontal", ticks = [-colorbar_scale[j], 0, colorbar_scale[j]])
				
				cb.ax.set_xticklabels([f"{-colorbar_scale[j]:.1f}", "0.0", f"{colorbar_scale[j]:.1f}"])
				cb.ax.set_xlabel(r"$\rho(R, \phi)$", fontsize = 15)
			

	for i in range(nrows-1):
		for j in range(len(filestem)):
			axs[i,j].xaxis.set_ticks([])
	for i in range(nrows):
		for j in range(1,len(filestem)):
			axs[i,j].yaxis.set_ticks([])

	axs[0,0].set_title("Self Consistent", fontsize = 15)
	axs[0,2].set_title("Self Consistent", fontsize = 15)
	axs[0,1].set_title("Test Particle", fontsize = 15)
	axs[0,3].set_title("Test Particle", fontsize = 15)

	
	axs[0,-1].set_ylabel(r"$\Omega t=1.66$", fontsize = 15)
	axs[0,-1].yaxis.set_label_position("right")
	axs[1,-1].set_ylabel(r"$\Omega t=3.33$", fontsize = 15)
	axs[1,-1].yaxis.set_label_position("right")
	axs[2,-1].set_ylabel(r"$\Omega t=5.00$", fontsize = 15)
	axs[2,-1].yaxis.set_label_position("right")
	axs[3,-1].set_ylabel(r"$\Omega t=6.66$", fontsize = 15)
	axs[3,-1].yaxis.set_label_position("right")

	xCir, yCir = [2.5 * cos(t) for t in np.linspace(0, 2*pi)], [2.5 * sin(t) for t in np.linspace(0, 2*pi)]
	axs[-2,0].plot(xCir, yCir, color = 'firebrick', linestyle = '--', alpha = 0.75)
	axs[-2,2].plot(xCir, yCir, color = 'firebrick', linestyle = '--', alpha = 0.75)

	plt.figtext(0.31,0.93, "Leading Perturbation", va="center", ha="center", fontsize=17)
	plt.figtext(0.72,0.93,"Trailing Perturbation", va="center", ha="center", fontsize=17)

	plt.show()

	# Create the plots 

## Cold Limit Plot ##


def cold_limit_Plot():
	sheet_13 = readingInRealCSV("Swing_Data/Cold_Limit/Cold_Limit_Sheet_13.csv")

	sheet_13_c_2, sheet_13_c_4, sheet_13_c_8 = readingInRealCSV("Swing_Data/Cold_Limit/Disc_Consistent_13_2.csv"), readingInRealCSV("Swing_Data/Cold_Limit/Disc_Consistent_13_4.csv"), readingInRealCSV("Swing_Data/Cold_Limit/Disc_Consistent_13_8.csv")
	sheet_13_t_2, sheet_13_t_4, sheet_13_t_8 = readingInRealCSV("Swing_Data/Cold_Limit/Disc_Test_13_2.csv"), readingInRealCSV("Swing_Data/Cold_Limit/Disc_Test_13_4.csv"), readingInRealCSV("Swing_Data/Cold_Limit/Disc_Test_13_8.csv")

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig,axs = plt.subplots()

	#axs.plot((sheet_13[0,:]*0.5*2), sheet_13[1,:])
	# axs.plot(sheet_13[0,:]*0.5*4, sheet_13[2,:])
	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-0.00, vmax=0.30))
	axs.legend(title = r"$\ell$", title_fontsize = 15, fontsize = 12)


	axs.plot((sheet_13[0,:]*8), sheet_13[1,:], color = "firebrick", label = "Shearing Sheet")

	for i in range(np.shape(sheet_13_c_2)[0]):
		axs.scatter(sheet_13_c_2[i,0] *2, sheet_13_c_2[i,1]/sheet_13_t_2[i,1], marker = 'o', color = cmap.to_rgba((i+1)*0.25/(np.shape(sheet_13_c_2)[0]-1)))
		axs.scatter(sheet_13_c_2[i,0] *4, 1.5*sheet_13_c_4[i,1]/sheet_13_t_4[i,1], marker = 'x',color = cmap.to_rgba((i+1)*0.25/(np.shape(sheet_13_c_2)[0]-1)))
		axs.scatter(sheet_13_c_2[i,0] *8, 2.5*sheet_13_c_8[i,1]/sheet_13_t_8[i,1], marker = '*',color = cmap.to_rgba((i+1)*0.25/(np.shape(sheet_13_c_2)[0]-1)))
		if (i ==0):
			axs.scatter(sheet_13_c_2[i,0] *2, sheet_13_c_2[i,1]/sheet_13_t_2[i,1], label = r'$\ell = 2$', marker = 'o', color = cmap.to_rgba((i+1)*0.25/(np.shape(sheet_13_c_2)[0]-1)))
			axs.scatter(sheet_13_c_2[i,0] *4, 1.5*sheet_13_c_4[i,1]/sheet_13_t_4[i,1], label = r'$\ell = 4$', marker = 'x', color = cmap.to_rgba((i+1)*0.25/(np.shape(sheet_13_c_2)[0]-1)))
			axs.scatter(sheet_13_c_2[i,0] *8, 2.5*sheet_13_c_8[i,1]/sheet_13_t_8[i,1], label = r'$\ell = 8$', marker = '*', color = cmap.to_rgba((i+1)*0.25/(np.shape(sheet_13_c_2)[0]-1)))

	axs.set_xlabel(r"$\xi \ell$", fontsize = 15)
	axs.set_ylabel(r"$A_{max}$", fontsize = 15)

	cbar = fig.colorbar(cmap)
	cbar.set_label( r"$\sigma_{R}$", fontsize = 15) 
	axs.legend(fontsize = 12)

	plt.show()



#axisymmetric_comparison_plots()
#axisymmetric_comparison_plots(files = ["Swing_Data/NonAxisymmetric_Comparisons/TP_sheet_Nonaxisymmetric_G_pert.csv", "Swing_Data/NonAxisymmetric_Comparisons/TP_Linear_Nonaxisymmetric_G_pert.csv", "Swing_Data/NonAxisymmetric_Comparisons/SC_sheet_Nonaxisymmetric_G_pert.csv", "Swing_Data/NonAxisymmetric_Comparisons/SC_Linear_Nonaxisymmetric_G_pert.csv"], indices = [10,20, 30, 40], centre = 8, delta = 0.5)

#amplificationPlots()

#swing_plots()
#cold_limit_Plot()

#selfgrav_Amplification()

max_density_Response()

# den = TwoDdensity("Swing_Data/Amplification_Fixed_5/CC_Evolution_0.csv")
# plt.imshow(den[1]- den[2])
# plt.show()