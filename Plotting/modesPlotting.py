from Density_Classes.TwoDdensity import *
from generalPlottingFunctions import *
from matplotlib.cm import *
from CoefficientClass import *
from ModeFinder import *
from scipy.stats import linregress
import csv
from Density_Classes.WaveFitter import * 

import matplotlib.animation as animation


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
	
	for nu, color in zip([4,6,8], colors):
		chi, omega0, eta = calculatingOmega0Eta(f"Modes_Data/Chi_Search/Chi_SearchKernel_Mode_Searching_{nu}_" , step)
		for c, o, e in zip(chi, omega0, eta):
			print(f"{o}, {e}, {c}")
		axs.plot(chi, eta, label = str(nu), color = color)
	

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

def scarredModesDensityRadii(radius, rMax = 6):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	colors = ['firebrick','royalblue','navy']
	fig, axs = plt.subplots(ncols = 2, sharex = True, sharey=False)

	cmap = ScalarMappable(Normalize(0.9,2.2), cmap = "plasma")


	for r in radius:
		file = f"Modes_Data/Scarred_Density/SD_{10*r:.0f}_25_-95.csv"
		mode = TwoDdensity(file)

		rad, amp, phs = mode.fourierCoeffT(2, -1, rMax) 
		axs[0].plot((1/r)*rad, amp*(rad/(2*pi)), color = cmap.to_rgba(r), label = f"{r:.1f}")
		axs[1].plot((1/r)*rad, amp*(rad/(2*pi))*np.cos(phs), color = cmap.to_rgba(r), label = f"{r:.1f}")
		

	axs[0].set_xlim([0.5,3])
	
	axs[1].set_xlabel(r"$R/R_{in}$", fontsize = 15)
	axs[0].set_xlabel(r"$R/R_{in}$", fontsize = 15)

	axs[0].set_ylabel(r"$\rho(R)/\rho_{0}(R)$", fontsize = 15)
	axs[1].set_ylabel(r"$\rho(R,\phi = 0)/\rho_{0}(R)$", fontsize = 15)

	axs[0].legend(title = r"$R_{in}$", fontsize = 12, title_fontsize = 15)

	plt.show()


def modeEvolutionComplexComponents(filestem):

	pr, pi, nr, ni = TwoDdensity(filestem + "/pos_real.csv"), TwoDdensity(filestem + "/pos_imag.csv"), TwoDdensity(filestem + "/neg_real.csv"), TwoDdensity(filestem + "/neg_imag.csv")
	nr.fourierAnimations(2,rMax = 5)

def modeAnimationVaryingChi(filestem, to_add, filename = None):
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=1, metadata=dict(artist='Me'))

	

	fig, axs = plt.subplots(1,1)

	ims = []
	for xi in to_add:
		modes = ModeFinder(filestem + str(int(100*xi)) +".csv")
		XX, YY, logDet = np.real(modes.omegas), np.imag(modes.omegas), np.log(np.absolute(modes.det))
		xMin, xMax = XX[0,0], XX[-1,-1]
		yMin, yMax = YY[0,0], YY[-1,-1]
		#contourFilled = axs.contourf(XX, YY, logDet, levels = 100, animated = True)
		contourFilled = axs.imshow(logDet,origin = 'lower', extent = [xMin, xMax, yMin, yMax])


		title = fig.text(.4,.9,(f"Active Fraction: {xi}"))
		
		ims.append([contourFilled, title])


	ani = animation.ArtistAnimation(fig, ims, interval=30)
	if (filename):
		print("Animation saved to: " + filename)
		ani.save(filename, writer = writer)
	else:
		plt.show()
	
def coefficentEvolution(file = "Modes_Data/fitting_test.csv"): # Maybe Could include a shaded region for the amplitude? 
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	linear = readingInComplexCSV(file)

	fig,axs = plt.subplots(ncols = 2)
	lst_n = range(0,15,3)
	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-1, vmax=17))

	for n in lst_n:
		axs[0].plot(np.linspace(0, 150, np.shape(linear[:800,n])[0]), np.absolute(linear[:800,n]), color = cmap.to_rgba(n))
		deriv = getGradient(np.linspace(0, 150, np.shape(linear[:800,n])[0]), np.angle(linear[:800,n]), multiple = 0.55, includeMinus = 1)

		axs[1].plot(*deriv, label = f"{n}", color = cmap.to_rgba(n))

	axs[0].set_ylabel(r"$|B_{p}(t)|$", fontsize = 15)
	axs[1].set_ylabel(r"$\partial\arg \left(B_{p}\right)/\partial t$", fontsize = 15)

	axs[0].set_yscale('log')
	t = np.linspace(50,150)
	axs[0].plot(t, 0.025*np.exp(t * 0.22), label = r"$\exp{(\eta t)}$", color = 'k', linestyle = '--')

	
	axs[1].axhline(-0.88, color = 'k', linestyle = '--', label  = r'$-\omega_{0}$')
	axs[1].legend(title = r"Index, $p$", title_fontsize = 15, fontsize =12)
	axs[0].legend(fontsize = 12)

	axs[0].set_title("Amplitude", fontsize = 15)
	axs[1].set_title("Phase", fontsize = 15)

	xlabelPositions, xlabels = [2*pi*i for i in range(0,25, 6)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,25, 6)]
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)
	axs[1].set_xticks(xlabelPositions)
	axs[1].set_xticklabels(xlabels)
	axs[0].set_xlabel("Time", fontsize = 15)
	axs[1].set_xlabel("Time", fontsize = 15)


	
	
	#axs[0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	plt.show()



def fitCoefficents2File(file = "Modes_Data/Chi_Coeff/UnstableModes.csv"):
	filename = lambda nu, chi :  f"Modes_Data/Chi_Coeff/coeff_{nu}_{100*chi:.0f}.csv"
	chi_lst, nu_4, nu_6, nu_8 = [], [], [], []
	for chi in np.arange(0.2, 1., 0.01):
		coeff_4, coeff_6, coeff_8 = CoefficientClass(filename(4, chi), 400), CoefficientClass(filename(6, chi), 400), CoefficientClass(filename(8, chi), 400)
		chi_lst.append(chi)
		nu_4.append(coeff_4.fitUnstableMode())
		nu_6.append(coeff_6.fitUnstableMode())
		nu_8.append(coeff_8.fitUnstableMode())

	# Now write to file 
	with open(file, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f)
		for c, n_4, n_6, n_8 in zip(chi_lst, nu_4, nu_6, nu_8):
			writer.writerow([c, n_4.real, n_4.imag, n_6.real, n_6.imag, n_8.real, n_8.imag])



def kernelMode2File(file = "Modes_Data/Chi_Coeff/UnstableModesKernel.csv"):
	
	chi_4, omega0_4, eta_4 = calculatingOmega0Eta(f"Modes_Data/Chi_Search/Chi_SearchKernel_Mode_Searching_4_" , 0.05)
	chi_6, omega0_6, eta_6 = calculatingOmega0Eta(f"Modes_Data/Chi_Search/Chi_SearchKernel_Mode_Searching_6_" , 0.05)
	chi_8, omega0_8, eta_8 = calculatingOmega0Eta(f"Modes_Data/Chi_Search/Chi_SearchKernel_Mode_Searching_8_" , 0.05)

	nu_4 = [o + 1j * n for o, n in zip(omega0_4, eta_4)]
	nu_6 = [o + 1j * n for o, n in zip(omega0_6, eta_6)]
	nu_8 = [o + 1j * n for o, n in zip(omega0_8, eta_8)]
	print(chi_4)

	with open(file, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f)
		for c, n_4, n_6, n_8 in zip(chi_4, nu_4, nu_6, nu_8):
			writer.writerow([c, n_4.real, n_4.imag, n_6.real, n_6.imag, n_8.real, n_8.imag])

	print(f"Data saved to: {file}")
		
		
	

def varyingChiPlot():
	data = readingInRealCSV("Modes_Data/Chi_Coeff/UnstableModes.csv")
	data_scatter = readingInRealCSV("Modes_Data/Chi_Coeff/UnstableModesKernel.csv")

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 1, sharex = True)
	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=0.4, vmax=1.2))
	
	n = 47
	for x0, x1, y0, y1, c in zip(data[n:,1], data[n+1:,1], data[n:,2], data[n+1:,2], data[n:,0]):
		axs.plot([x0, x1],[y0,y1], color = cmap.to_rgba(c), linestyle = '--')
	axs.plot([x0, x1],[y0,y1], color = cmap.to_rgba(c), linestyle = '--', label = 'Kernel')

	
	for x0, y0, c in zip(data_scatter[:7,1], data_scatter[:7,2], data_scatter[:7,0]):
		axs.scatter(x0, y0+0.01, color = cmap.to_rgba(c))
	axs.scatter(x0, y0+0.01, color = cmap.to_rgba(c), label ='Coefficient')
	
	# for x0, x1, y0, y1, c in zip(data[35:,3], data[36:,3], data[35:,4], data[36:,4], data[35:,0]):
	# 	axs.plot([x0, x1],[y0,y1], color = cmap.to_rgba(c))
	n = 32
	for x0, x1, y0, y1, c in zip(data[n:,5], data[n+1:,5], data[n:,6], data[n+1:,6], data[n+1:,0]):
		axs.plot([x0, x1],[y0,y1], color = cmap.to_rgba(c))
	axs.plot([x0, x1],[y0,y1], color = cmap.to_rgba(c), label = 'Kernel')
	
	for x0, y0, c in zip(data_scatter[:10,5], data_scatter[:10,6], data_scatter[:10,0]):
		axs.scatter(x0, y0, color = cmap.to_rgba(c), marker = 'x')
	axs.scatter(x0, y0, color = cmap.to_rgba(c), marker = 'x', label ='Coefficient')
	
	

	# axs[0].plot(data[47:,0], data[47:,1])
	# axs[0].plot(data[35:,0], data[35:,3])
	# axs[0].plot(data[27:,0], data[27:,5])

	# axs[1].plot(data[47:,0], data[47:,2])
	# axs[1].plot(data[35:,0], data[35:,4])
	# axs[1].plot(data[27:,0], data[27:,6])


	axs.axhline(0, linestyle = '--', color = 'k')

	
	axs.set_ylabel(r"$\eta$", fontsize = 15)
	axs.set_xlabel(r"$\omega_{0}$", fontsize = 15)
	cbar = fig.colorbar(cmap, ax=axs).set_label(label = r"$\xi$", fontsize = 15)

	h, l = axs.get_legend_handles_labels()
	ph = [plt.plot([],marker="", ls="")[0]]*2
	handles = ph + h
	labels = [r"$\nu_{t}$ = 4:", r"$\nu_{t}$ = 8:"] + l
	plt.legend(handles, labels, ncol=3, fontsize =12)
	

	plt.show()


## In-going and Out-going ##

def seperatingWaves():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(ncols = 2, nrows = 2, sharey = 'row', sharex = True)
	wf = WaveFitter("Modes_Data/Mode_Evolution_R_20_W_25_D_-95_G.csv", timeEnd = 25)
	wf.split_waves()


	axs[0,0].plot(wf.radii, np.absolute(wf.outgoing_T(50)), linestyle = '--', label = r"$\rho(R)$", color = 'firebrick')
	axs[0,0].plot(wf.radii, np.absolute(wf.outgoing_T(50))*np.cos(np.angle(wf.outgoing_T(50))), color = 'royalblue')
	
	axs[0,1].plot(wf.radii, np.absolute(wf.ingoing_T(50)), linestyle = '--', color = 'firebrick', label = r"$\rho(R)$")
	axs[0,1].plot(wf.radii, np.absolute(wf.ingoing_T(50))*np.cos(np.angle(wf.ingoing_T(50))), color = 'royalblue', label = r"$\rho(R, \phi = 0)$")

	axs[1,0].plot(*getGradient(wf.radii, np.angle(wf.outgoing_T(50)), includeMinus = 1), color = 'royalblue')
	axs[1,1].plot(*getGradient(wf.radii, np.angle(wf.ingoing_T(50)), includeMinus = 1), color = 'royalblue')
	
	axs[1,0].axhline(0, linestyle = '--', color = 'k')
	axs[1,1].axhline(0, linestyle = '--', color = 'k')

	axs[0,0].axvline(2, linestyle = ':', color = 'gray')
	axs[0,1].axvline(2, linestyle = ':', color = 'gray')
	axs[1,0].axvline(2, linestyle = ':', color = 'gray')
	axs[1,1].axvline(2, linestyle = ':', color = 'gray')

	axs[0,0].axvline(2.71, linestyle = '--', color = 'gray')
	axs[0,1].axvline(2.71, linestyle = '--', color = 'gray')
	axs[1,0].axvline(2.71, linestyle = '--', color = 'gray')
	axs[1,1].axvline(2.71, linestyle = '--', color = 'gray')

	axs[0,0].axvline(5.5, linestyle = ':', color = 'gray')
	axs[0,1].axvline(5.5, linestyle = ':', color = 'gray')
	axs[1,0].axvline(5.5, linestyle = ':', color = 'gray')
	axs[1,1].axvline(5.5, linestyle = ':', color = 'gray')

	#axs[0,0].set_xlim([1., 5.6])
	axs[0,0].set_ylim([-0.1, 0.2])

	axs[0,0].set_title("Out-going", fontsize = 15)
	axs[0,1].set_title("In-going", fontsize = 15)

	axs[0,0].set_ylabel(r"Density", fontsize = 15)
	axs[1,0].set_ylabel(r"$k(R)$", fontsize = 15)
	axs[0,1].legend(fontsize = 12)

	axs[1,0].set_xlabel(r"Radius", fontsize = 15)
	axs[1,1].set_xlabel(r"Radius", fontsize = 15)

	plt.show()


















#eigenModeComparison(["Modes_Data/JB_mode.csv", "Modes_Data/JB_mode_Time.csv"], 0.88)

#Modes_Data/JB_mode_Coeff.csv

#powerPlots([f"Modes_Data/Coefficent_JB_Mode_{i}.csv" for i in [5, 7, 9, 10]], 100)
#Coefficent_Taper_4
#powerPlots([f"Modes_Data/Coefficent_Taper_{i}.csv" for i in [8,6,4]], 100)

#calculatePeriods([f"Modes_Data/Coefficent_Taper_{i}.csv" for i in [8,6,4]], 100)

#modesVaryingChi([f"Modes_Data/Chi_SearchKernel_Mode_Searching_8_{100-i*5}.csv" for i in range(0, 18)])


#modesVaryingChi()
#modes = ModeFinder("Modes_Data/Chi_Search/Chi_SearchKernel_Mode_Searching_4_100.csv")
#modes.contourPlotShow()

#varyingScarPosition()
#modesVaryingChi()

#comparingXi([f"Modes_Data/Xi_Modes/RM_Xi_{x}_Nu_6.out" for x in ['5','6','7','8']], [0.74661, 0.77, 0.8, 0.83])

#scarredModesDensity([f"Modes_Data/Scarred_Density/SD_{r}_25_-95.csv" for r in ['12', '14','16','18','20']], [1.2, 1.4, 1.6, 1.8, 2.0], [0.564407, 0.510169, 0.469492, 0.428814, 0.394915])
#scarredModesDensity([f"Modes_Data/Scarred_Density/SD_{r}_25_-95.csv" for r in ['14','16', '18']], [1.4, 1.6, 1.8], [0.510169, 0.469492, 0.428814])

#scarredModesDensityRadii([1.2, 1.4, 1.6, 1.8, 2.0])
#scarredModesDensityRadii([f"Modes_Data/Scarred_Density/SD_{r}_25_-95.csv" for r in ['12','18','20']], [1.2, 1.8, 2.0], [0.564, 0.428814, 0.394915])


#modeEvolutionComplexComponents("Modes_Data/Real_Imag_Modes")

#modes = ModeFinder("Modes_Data/Chi_Search/VideoKernel_Mode_Searching_4_100.csv")
#modes.contourPlotShow()
#modeAnimationVaryingChi("Modes_Data/Chi_Search/VideoKernel_Mode_Searching_4_", [0.92, 0.94, 0.96, 0.98, 1])
#coefficentEvolution()
#kernelMode2File()
#fitCoefficents2File()
#varyingChiPlot()
#coefficentEvolution("Modes_Data/Chi_Coeff/coeff_4_69.csv")

seperatingWaves()