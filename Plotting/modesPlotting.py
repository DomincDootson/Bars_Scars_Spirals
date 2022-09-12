from Density_Classes.TwoDdensity import *
from generalPlottingFunctions import *
from CoefficientClass import *
from ModeFinder import *
from scipy.stats import linregress
import csv


def resonance(omega, m1, m2):
	resRadius = (m2 + sqrt(2)*m1) /omega # assume vc = 1

	return [[resRadius * cos(theta) for theta in np.linspace(0, 2*3.14)], [resRadius * sin(theta) for theta in np.linspace(0, 2*3.14)]] 


def eigenModeComparison(files, omega, rMax = 5, rInner = 1):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig,axs = plt.subplots(ncols = len(files), sharey = True, sharex = True)

	for i, file in enumerate(files):
		mode = TwoDdensity(file)
	

		spacing = 2*rMax/(mode.nCols-1)
		centre = mode.centre

		x = np.arange(-rMax,rMax+spacing, spacing)
		y = np.arange(-rMax,rMax+spacing, spacing)
		XX, YY = np.meshgrid(x, y)

		maxValue = mode.maxDensityEvolution()[-1]
		contours = [frac * maxValue for frac in [.2, .5, .8]]

		ILR, CR, OLR = resonance(omega, -1, 2), resonance(omega, 0, 2), resonance(omega, 1, 2)
		rI = [[rInner * cos(theta) for theta in np.linspace(0, 2*3.14)], [rInner * sin(theta) for theta in np.linspace(0, 2*3.14)]] 
		
		axs[i].contourf(XX, YY, mode.densityAtTime(-1), levels = 50)
		axs[i].contour(XX, YY, mode.densityAtTime(-1), levels = contours, colors = 'k')
		axs[i].plot(ILR[0], ILR[1], color = 'firebrick', linestyle ='--')
		axs[i].plot(CR[0], CR[1], color = 'firebrick', linestyle ='--')
		axs[i].plot(OLR[0], OLR[1], color = 'firebrick', linestyle ='--')
		axs[i].plot(rI[0], rI[1], color = 'white', linestyle ='-')
		axs[i].set_aspect('equal')


	axs[0].set_title("Modal Calculation", fontsize = 15)
	axs[1].set_title("Temporal Calculation", fontsize = 15)
	#axs[2].set_title("Residuals", fontsize = 15)

	plt.show()


def comparingXi(files, omega, rMax = 5):
	fig,axs = plt.subplots(nrows = len(files), sharey = True, sharex = True)

	for i, file in enumerate(files):
		mode = TwoDdensity(file)
	

		spacing = 2*rMax/(mode.nCols-1)
		centre = mode.centre

		x = np.arange(-rMax,rMax+spacing, spacing)
		y = np.arange(-rMax,rMax+spacing, spacing)
		XX, YY = np.meshgrid(x, y)

		maxValue = mode.maxDensityEvolution()[-1]
		contours = [frac * maxValue for frac in [.2, .5, .8]]

		ILR, CR, OLR = resonance(omega[i], -1, 2), resonance(omega[i], 0, 2), resonance(omega[i], 1, 2)
		
		axs[i].contourf(XX, YY, mode.densityAtTime(-1), levels = 50)
		axs[i].contour(XX, YY, mode.densityAtTime(-1), levels = contours, colors = 'k')
		axs[i].plot(ILR[0], ILR[1], color = 'firebrick', linestyle ='--')
		axs[i].plot(CR[0], CR[1], color = 'firebrick', linestyle ='--')
		axs[i].plot(OLR[0], OLR[1], color = 'firebrick', linestyle ='--')
		axs[i].set_aspect('equal')

	plt.show()


def expoFit(time, power): # Put this function into the complex class 
	logPower = np.log(power)

	fit = linregress(time, logPower)

	return time, np.exp(fit.intercept + fit.slope * time), fit.slope/2

def calculatePeriods(files, timeEnd):
	coeffs = [CoefficientClass(file) for file in files]
	time = np.linspace(0, timeEnd, coeffs[0].nTime)
	print([coeff.periodFitter(time) for coeff in coeffs])
	

def powerPlots(files, timeEnd):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots()

	coeffs = [CoefficientClass(file) for file in files]
	axs.set_yscale('log')

	for i, coeff in enumerate(coeffs):
		time = np.linspace(0, timeEnd, coeff.nTime)
		axs.plot(time, coeff.totalModePower(), color = "royalblue")

		fitTime, fit, fitSlope = coeff.expoFit(time)
		axs.plot(fitTime, fit, linestyle = '--', label = f"{fitSlope:.2f}", color = "firebrick")


	axs.legend(title = r"$\eta$")
	axs.set_ylabel(r"Response Power, $B_{q}B^{q}$", fontsize = 15)
	pos = [0, 4*(2*pi), 8*(2*pi), 12*(2*pi), 16*(2*pi)]
	label = [0, r"4$\times 2\pi$", r"8$\times 2\pi$", r"12$\times 2\pi$", r"16$\times 2\pi$"]

	plt.xticks(pos, label, fontsize = 12)
	axs.set_xlabel("Time", fontsize = 15)
	plt.show()


def calculatingOmega0Eta(filestem, step):
	files = [filestem + f"{100-i*5}.csv" for i in range(0, 18)]
	eigenFreq = [ModeFinder(file).mostUnstableMode() for file in files]
	chi, omega0, eta = [1 - i *step for i in range(0, len(files))], [omega.real for omega in eigenFreq], [omega.imag for omega in eigenFreq]

	return chi, omega0, eta

def modesVaryingChi(step = 0.05): # save the files here 
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 1)
	colors = ["royalblue", "darkviolet", "firebrick"]
	axs.tick_params(size = 12)
	
	for nu, color in zip([6], colors):
		chi, omega0, eta = calculatingOmega0Eta(f"Modes_Data/Chi_Search/Chi_SearchKernel_Mode_Searching_{nu}_" , step)
		for c, o, e in zip(chi, omega0, eta):
			print(f"{o}, {e}, {c}")
		axs.plot(chi, omega0, label = str(nu), color = color)
	

	plt.rcParams['legend.title_fontsize'] = 15
	axs.set_ylabel(r"$\eta$", fontsize = 15)
	axs.set_xlabel(r"$\xi$", fontsize = 15)
	axs.legend(title = r"$\nu_{t}$", fontsize = 15)

	axs.axhline(color = 'k', linestyle = '--')
	


	plt.show()

def varyingScarPosition(R = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 28]):

	fig, axs = plt.subplots(nrows = 2, ncols = 5)

	for (i,j), ax in np.ndenumerate(axs):
		file = f"Modes_Data/Kernel_Search/Single_Scar/AM_Scarred_Kernel_R_{R[i*5 + j]}_W_25_D_10.csv"
		det = ModeFinder(file)
		ax = det.contourPlot(axs = ax)

	plt.show()


def savingModes(radii, filename = 'VaryingInnerPositionResponse.csv', readinFileStem = "Modes_Data/Kernel_Search/Single_Scar/AM_Scarred_Kernel_R_" ):
	with open(filename, 'w') as f:
		files = [readinFileStem + str(int(radius*10)) + "_W_25_D_-95_G.csv" for radius in radii]
		for file, radius in zip(files, radii):
			print(file)
			mode = ModeFinder(file)
			omegas = mode.modes()
			print(omegas)
			for omega in omegas:
				f.write(f"{radius},{omega.real},{omega.imag}\n")
		

## Plotting Scarred Eigenmodes ## 
## --------------------------- ##

def scarredModesDensity(files, radius, omega0s, rMax = 6):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = len(files), sharex = True, sharey=True)
	
	for file, ax, r, om in zip(files, axs, radius, omega0s):
		mode = TwoDdensity(file)
	
		spacing = 2*rMax/(mode.nCols-1)
		centre = mode.centre

		x = np.arange(-rMax,rMax+spacing, spacing)
		y = np.arange(-rMax,rMax+spacing, spacing)
		XX, YY = np.meshgrid(x, y)

		maxValue = mode.maxDensityEvolution()[-1]
		contours = [frac * maxValue for frac in [.2, .5, .8]]

		ILR, CR, OLR = resonance(om, -1, 2), resonance(om, 0, 2), resonance(om, 1, 2)
		
		ax.contourf(XX, YY, mode.densityAtTime(-1), levels = 100)
		ax.contour(XX, YY, mode.densityAtTime(-1), levels = contours, colors = 'k')
		
		ax.plot(ILR[0], ILR[1], color = 'firebrick', linestyle ='--')
		ax.plot(CR[0], CR[1], color = 'firebrick', linestyle ='--')
		#ax.plot(OLR[0], OLR[1], color = 'firebrick', linestyle ='--')
		ax.plot([r *cos(th) for th in np.linspace(0, 2 * pi)], [r * sin(th) for th in np.linspace(0, 2 * pi)], linestyle = ':', color = "white")
	ax.set_title(f"$R_{{i}} = {r}$", fontsize = 12)


		ax.set_aspect('equal')

	fig.tight_layout()

	plt.show()




#eigenModeComparison(["Modes_Data/JB_mode.csv", "Modes_Data/JB_mode_Time.csv"], 0.88)

#Modes_Data/JB_mode_Coeff.csv

#powerPlots([f"Modes_Data/Coefficent_JB_Mode_{i}.csv" for i in [5, 7, 9, 10]], 100)
#Coefficent_Taper_4
#powerPlots([f"Modes_Data/Coefficent_Taper_{i}.csv" for i in [8,6,4]], 100)

#calculatePeriods([f"Modes_Data/Coefficent_Taper_{i}.csv" for i in [8,6,4]], 100)

#modesVaryingChi([f"Modes_Data/Chi_SearchKernel_Mode_Searching_8_{100-i*5}.csv" for i in range(0, 18)])


#modesVaryingChi()

#varyingScarPosition()
#modesVaryingChi()

#comparingXi([f"Modes_Data/Xi_Modes/RM_Xi_{x}_Nu_6.out" for x in ['5','6','7','8']], [0.74661, 0.77, 0.8, 0.83])

scarredModesDensity([f"Modes_Data/Scarred_Density/SD_{r}_25_-95.csv" for r in ['12', '14','16','18','20']], [1.2, 1.4, 1.6, 1.8, 2.0], [0.564407, 0.510169, 0.469492, 0.428814, 0.394915])
#scarredModesDensity([f"Modes_Data/Scarred_Density/SD_{r}_25_-95.csv" for r in ['14','16', '18']], [1.4, 1.6, 1.8], [0.510169, 0.469492, 0.428814])





