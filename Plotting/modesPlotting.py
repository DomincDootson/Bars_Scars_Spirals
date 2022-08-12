from Density_Classes.TwoDdensity import *
from generalPlottingFunctions import *
from CoefficientClass import *
from scipy.stats import linregress

def resonance(omega, m1, m2):
	resRadius = (m2 + sqrt(2)*m1) /omega # assume vc = 1

	return[[resRadius * cos(theta) for theta in np.linspace(0, 2*3.14)], [resRadius * sin(theta) for theta in np.linspace(0, 2*3.14)]] 


def eigenModeComparison(files, omega, rMax = 5):
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
		
		axs[i].contourf(XX, YY, mode.densityAtTime(-1), levels = 50)
		axs[i].contour(XX, YY, mode.densityAtTime(-1), levels = contours, colors = 'k')
		axs[i].plot(ILR[0], ILR[1], color = 'firebrick', linestyle ='--')
		axs[i].plot(CR[0], CR[1], color = 'firebrick', linestyle ='-')
		axs[i].plot(OLR[0], OLR[1], color = 'firebrick', linestyle ='--')
		axs[i].set_aspect('equal')


	axs[0].set_title("Modal Calculation", fontsize = 15)
	axs[1].set_title("Temporal Calculation", fontsize = 15)
	#axs[2].set_title("Residuals", fontsize = 15)

	plt.show()

def expoFit(time, power):
	logPower = np.log(power)

	fit = linregress(time, logPower)

	return time, np.exp(fit.intercept + fit.slope * time), fit.slope/2 

def powerPlots(files, timeEnd):
	fig, axs = plt.subplots()

	coeffs = [CoefficientClass(file) for file in files]
	axs.set_yscale('log')

	for i, coeff in enumerate(coeffs):
		time = np.linspace(0, timeEnd, coeff.nTime)
		axs.plot(time, coeff.totalModePower())

		fitTime, fit, fitSlope = expoFit(time[coeff.nTime//3:], coeff.totalModePower()[coeff.nTime//3:])
		axs.plot(fitTime, fit, linestyle = '--', label = f"{fitSlope:.2f}")


	axs.legend(title = r"$\eta$")
	plt.show()


#eigenModeComparison(["Modes_Data/JB_mode.csv", "Modes_Data/JB_mode_Time.csv"], 0.88)

#Modes_Data/JB_mode_Coeff.csv

powerPlots([f"Modes_Data/Coefficent_JB_Mode_{i}.csv" for i in [5, 7, 9, 10]], 100)