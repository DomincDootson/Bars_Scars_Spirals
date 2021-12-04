from generalPlottingFunctions import *

def circle(radius):
	return [radius * cos(theta) for theta in np.linspace(0,2*3.14)], [radius * sin(theta) for theta in np.linspace(0,2*3.14)]

def greensFunctionEvolution(): 
	flatternedDensity = list(readingInRealCSV("../Disk_Kicking/littleSigma_35/Density2D20_2.csv"))
	nRows, nCols = int(sqrt(np.size(flatternedDensity[0]))), int(sqrt(np.size(flatternedDensity[0]))) # Assume square
	density2D = [np.reshape(array[:nRows*nCols], (nRows, nCols,)) for array in flatternedDensity] 	

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(2, 4, sharex = True, sharey = True) 
	spacing = 10/(nCols-1)
	centre = (nCols-1)*0.5

	time = [15, 30, 45, 60, 75, 90, 105, 120]

	x = np.arange(-5,5+spacing, spacing)
	y = np.arange(-5,5+spacing, spacing)
	XX, YY = np.meshgrid(x, y)

	xCir, yCir = circle(1.95)
	for i in range(2):
		for j in range(4):
			axs[i,j].contourf(XX, YY, density2D[time[j + i * 4]], 100)
			if (i==0 and j ==2):
				contourFilled = axs[i,j].contourf(XX, YY, density2D[time[j + i * 4]], 100) 
			axs[i,j].plot(xCir, yCir, color = 'firebrick', alpha = 0.8, linestyle = '--')
			axs[i,j].text(-4,4, str(0.25*time[j + i * 4]) +r"$t_{0}$", color = 'black')
	
	fig.tight_layout()
	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	fig.colorbar(contourFilled, cax=cbar_ax)

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


kalnajsPlots()
#greensFunctionEvolution()