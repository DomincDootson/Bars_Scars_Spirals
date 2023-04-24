from generalPlottingFunctions import *
from Density_Classes.TwoDdensity import *
import matplotlib.ticker as ticker

def circle(radius):
	return [radius * cos(theta) for theta in np.linspace(0,2*3.14)], [radius * sin(theta) for theta in np.linspace(0,2*3.14)]

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if x ==0:
    	return "0"
    return r'${} \times 10^{{{}}}$'.format(a, b)

def greensFunctionEvolution(filename = "../Disk_Kicking/littleSigma_35/Density2D20_2.csv"): 
	density2D = TwoDdensity(filename)

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(1, 4, sharex = True, sharey = True) 
	spacing = 10/(density2D.nCols-1)
	centre = (density2D.nCols-1)*0.5
					
	x = np.arange(-5,5+spacing, spacing)
	y = np.arange(-5,5+spacing, spacing)
	XX, YY = np.meshgrid(x, y)
	time = [1, 24, 44, 64]

	xCir, yCir = circle(1.95)
	absMaxValue = density2D.maxDensityEvolution()[5]


	for j in range(4):
			
		#axs[i,j].imshow(density2D.densityAtTime(time[i*4 + j]))
		if j==0:
			contourFilled = axs[j].contourf(XX, YY, density2D.densityAtTime(time[j]), levels = 100, vmin = -absMaxValue, vmax = absMaxValue)
		axs[j].contourf(XX, YY, density2D.densityAtTime(time[j]), levels = 100, vmin = -absMaxValue, vmax = absMaxValue)
		axs[j].set(aspect = 1)
		if j==0:
			axs[j].set_title(r"$t=0$")
		else:
			axs[j].set_title(r"$t=$ " + str(time[j] *0.25) + r"$\times 2\pi$")
	
	axs[0].plot(xCir, yCir, color = 'firebrick', alpha = 0.8, linestyle = '--')
			


		
	
	fig.tight_layout()
	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	cbar = fig.colorbar(contourFilled, cax=cbar_ax, format=ticker.FuncFormatter(fmt))
	cbar.set_ticks([-1.5*(10**-3), 0, 1.5*(10**-3)])

	plt.show()

def kalnajsPlots():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	den, pot = readingInRealCSV("KalnajsTesting/kalnajsDensity.csv"), readingInRealCSV("KalnajsTesting/kalnajsPotential.csv")
	r = np.linspace(0,1, np.shape(den)[0])

	labelP, labelD = [r"$\mathcal{U}^{2}_{0}(r)$", r"$\mathcal{U}^{2}_{1}(r)$", r"$\mathcal{U}^{2}_{2}(r)$"], [r"$\mathcal{D}^{2}_{0}(r)$", r"$\mathcal{D}^{2}_{1}(r)$", r"$\mathcal{D}^{2}_{2}(r)$"]
	colorP, colorD = ["cornflowerblue", "#1f77b4", "midnightblue"], ["firebrick", "indianred","lightcoral"]

	fig, axs = plt.subplots(nrows = 2, ncols = 1, sharex = True)

	for i in range(3):
		axs[0].plot(r, pot[:,i], label = labelP[i], color = colorP[i])
		axs[1].plot(r, den[:,i], label = labelD[i], color = colorD[i])
	


	axs[0].legend(loc = "upper left", fontsize = 14)
	axs[1].legend(loc = "lower left", fontsize = 14)

	axs[0].set_ylabel("Potential", fontsize = 14)
	axs[1].set_ylabel("Density", fontsize = 14)

	axs[1].set_xlabel("Radius", fontsize = 14)
	
	axs[1].set_xticks([0, 0.25, 0.5, 0.75, 1])
	axs[1].set_xticklabels(["0", r"$r_{ka}/4$", r"$r_{ka}/2$", r"$3r_{ka}/4$", r"$r_{ka}$"]) 

	plt.show()


def greenComparison():
	fig, axs = plt.subplots(nrows = 4, ncols = 2, sharex = True, sharey = 'col')
	t_lst = [0, 2,4,6]

	data = [TwoDdensity("Greens_Data/Kalnajs_Green_2.csv"), TwoDdensity("Greens_Data/Kalnajs_Green_2.csv"), TwoDdensity("Greens_Data/Kalnajs_Test_Green_2.csv"), TwoDdensity("Greens_Data/Kalnajs_Test_Green_2.csv")]
	r = np.linspace(0, 5, data[0][0].shape[0])
	
	for i in range(4):
		t = t_lst[i]
		r, amp, phase = data[0].fourierCoeffT(2, t, 10)
		#axs[i,0].imshow(data[0][t])
		axs[i, 0].plot(r, amp, color = 'firebrick', linestyle = '--')
		axs[i, 0].plot(r, amp*np.cos(phase), color = 'firebrick', label = 'Kalnajs')

		r, amp, phase = data[1].fourierCoeffT(2, t, 10)
		# axs[i,1].imshow(data[1][t])
		axs[i, 0].plot(r, amp, color = 'royalblue', linestyle = '--')
		axs[i, 0].plot(r, amp*np.cos(phase), color = 'royalblue', label = 'Gaussian')

		r, amp, phase = data[2].fourierCoeffT(2, t, 10)
		#axs[i,2].imshow(data[2][t])
		axs[i, 1].plot(r, amp, color = 'firebrick', linestyle = '--')
		axs[i, 1].plot(r, amp*np.cos(phase), color = 'firebrick')

		r, amp, phase = data[3].fourierCoeffT(2, t, 10)
		# axs[i,3].imshow(data[3][t])
		axs[i, 1].plot(r, amp, color = 'royalblue', linestyle = '--')
		axs[i, 1].plot(r, amp*np.cos(phase), color = 'royalblue')

	axs[0,0].set_xlim([0,5])
	
	axs[-1,0].set_xlabel(r"$R$", fontsize = 15)
	axs[-1,1].set_xlabel(r"$R$", fontsize = 15)

	axs[0,0].set_ylabel(r"$t = 0$", fontsize = 15)
	axs[1,0].set_ylabel(r"$t = 4.1\times 2\pi$", fontsize = 15)
	axs[2,0].set_ylabel(r"$t = 8.3\times 2\pi$", fontsize = 15)
	axs[3,0].set_ylabel(r"$t = 12.4\times 2\pi$", fontsize = 15)

	axs[0,1].set_title(r"Test Particle", fontsize = 15)
	axs[0,0].set_title(r"Self Consistent", fontsize = 15)

	axs[0,0].legend(fontsize = 12)



	plt.show()



## Decomposition of the Greens Functions ##
## ------------------------------------- ##

def poissonKernel(r, rp, phi, ep=0):
	return -1/sqrt(r**2 +rp**2 - 2*r*rp*cos(phi) +ep**2)

def poissonKernelHarmonic(r, rp, l =2, ep = 0):
	phi = np.linspace(0, 2*pi, 100)
	integrand = [e**(-1j * l * p)*poissonKernel(r, rp, p, ep) for p in phi]
	
	return np.real(sum(integrand)) * phi[1]/(2*pi)


def interactionPotential():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 2, sharex = True, sharey = True)

	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-0.0, vmax=10))
	data_72, data_48  = readingInRealCSV("Greens_Data/Kalnajs_Interaction_Potential_72.csv"), readingInRealCSV("Greens_Data/Kalnajs_Interaction_Potential_48.csv")
	r_pd = np.linspace(0, 15, 502)
	r_lst = np.linspace(0, 15, 500)
	for r in reversed(range(1,8,2)):
		ker_48 = [poissonKernelHarmonic(r_har, r, ep = 0.1) for r_har in r_lst]
		
		axs[0].plot(r_lst, ker_48, linestyle ='--', color = cmap.to_rgba(r))
		axs[0].plot(r_pd, data_48[int(r*33), :], color = cmap.to_rgba(r))
		axs[0].set_xlim([0,10])

		ker_72 = [poissonKernelHarmonic(r_har, r, ep = 0.08) for r_har in r_lst]
		axs[1].plot(r_lst, ker_72, linestyle ='--', color = cmap.to_rgba(r))
				
		axs[1].plot(r_pd, data_72[int(r*33), :], color = cmap.to_rgba(r))
		axs[1].set_xlim([0,10])

	axs[0].plot(r_pd, data_48[int(r*33), :], color = cmap.to_rgba(r), label = 'Kalnajs Interaction')
	axs[0].plot(r_lst, ker_48, linestyle ='--', color = cmap.to_rgba(r), label = "Poisson Kernel")

	axs[0].legend(fontsize = 12)

	axs[0].set_xlabel(r"$R'$", fontsize = 15)
	axs[1].set_xlabel(r"$R'$", fontsize = 15)

	axs[0].set_title(r"$N=48$", fontsize = 15)
	axs[1].set_title(r"$N=72$", fontsize = 15)

	axs[0].set_ylabel(r"$U(R,R')$", fontsize = 15)

	cbar = fig.colorbar(cmap, ax=axs[:])
	cbar.ax.set_ylabel(r"$R$", rotation =270, fontsize = 15)

	plt.show()

interactionPotential()
#greenComparison()
#greensFunctionEvolution()
#greenComparison()