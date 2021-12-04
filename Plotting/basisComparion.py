from generalPlottingFunctions import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable


from scipy.stats import linregress
from scipy.optimize import minimize 

def plottingPerturbation():
	fig, axs = plt.subplots(nrows=1, ncols=2)

	radii = readingInRealCSV("BF_Comparison/radii.csv")
	pertK, pertG = readingInRealCSV("BF_Comparison/kalnajsPerturbationK.csv"), readingInRealCSV("BF_Comparison/kalnajsPerturbationG.csv")
	

	axs[0].plot(radii, pertK, color = 'royalblue', linestyle = '--', label = 'Kalnajs Perturbation')
	axs[0].plot(radii, pertG, color = 'firebrick', label = 'Gaussian Fit')
	axs[0].set_xlabel(r"Radius $[r_{0}]$")
	axs[0].set_ylabel(r"Potential")
	axs[0].set_title(r"Kalnajs Perturbation")
	axs[0].legend(loc = 'upper right')
	axs[0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	axs[0].axvline(15)

	pertK, pertG = readingInRealCSV("BF_Comparison/gaussianPerturbationK.csv"), readingInRealCSV("BF_Comparison/gaussianPerturbationG.csv")
	axs[1].plot(radii, pertK, color = 'royalblue', label = 'Kalnajs Fit')
	axs[1].plot(radii, pertG, color = 'firebrick', linestyle = '--', label = "Gaussian Perturbation")
	axs[1].set_xlabel(r"Radius $[r_{0}]$")
	axs[1].set_title(r"Gaussian Perturbation")
	axs[1].legend(loc = 'lower right')
	axs[1].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))

	plt.show()


def individualGaussian():
	data = readingInRealCSV("BF_Comparison/individualGaussian.csv")

	radii = readingInRealCSV("BF_Comparison/radii.csv")
	pertK, pertG = readingInRealCSV("BF_Comparison/kalnajsPerturbationK.csv"), readingInRealCSV("BF_Comparison/kalnajsPerturbationG.csv")


	lst = range(1,np.shape(data)[1])
	for n in lst:
		plt.plot(data[:,0], data[:, n])


	plt.plot(radii, pertK, color = 'black', linestyle='--', label = 'Kalnajs')
	plt.plot(radii, pertG, color = 'black', label = 'Fit')

	plt.legend()
	plt.show()


def bfComparison(filenames):
	kalnajs, gaussian = twoDdensity(filenames[0]), twoDdensity(filenames[1])

	times = [50, 100, 150,-1]
	nrows = len(times)

	fig, axs = plt.subplots(nrows = nrows, ncols = 2)
	
	spacing = 20/(kalnajs.nCols-1) # 
	centre = (kalnajs.nCols-1)*0.5

	x = np.arange(-10,10+spacing, spacing)
	y = np.arange(-10,10+spacing, spacing)
	XX, YY = np.meshgrid(x, y)

	for i in range(nrows):
		pcm1 = axs[i, 0].contourf(XX, YY, kalnajs.densityAtTime(times[i]), levels = 100)
		pcm2 = axs[i, 1].contourf(XX, YY, gaussian.densityAtTime(times[i]), levels = 100)

	axs[0,0].set_title("Kalnajs")
	axs[0,1].set_title("Gaussian")

	axins1 = inset_axes(axs[-1, 0], width="100%", height="5%", loc='lower center',borderpad=-5)
	axins2 = inset_axes(axs[-1, 1], width="100%", height="5%", loc='lower center',borderpad=-5) 

	fig.colorbar(pcm1, ax = axs[-1,0], cax=axins1, orientation="horizontal", format='sci')
	fig.colorbar(pcm2, ax = axs[-1,1], cax=axins2, orientation="horizontal", format='%.0e')




	plt.show()


def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if x ==0:
    	return "0"
    return r'${} \times 10^{{{}}}$'.format(a, b)

def densityComparison(filenames): 
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	kalnajsK, gaussianK, kalnajsG, gaussianG = twoDdensity(filenames[0]), twoDdensity(filenames[1]), twoDdensity(filenames[2]), twoDdensity(filenames[3])
	times = [100, 150,-1]


	fig = plt.figure()
	axs = [fig.add_subplot(len(times)+1, 2, 1), fig.add_subplot(len(times)+1, 2, 2)]

	for i in range(len(times)):
		for j in range(len(filenames)):
			axs.append(fig.add_subplot(len(times)+1, len(filenames), 5 + j + i * (len(filenames))))


	radii = readingInRealCSV("BF_Comparison/radii.csv")
	pertK, pertG = readingInRealCSV("BF_Comparison/kalnajsPerturbationK.csv"), readingInRealCSV("BF_Comparison/kalnajsPerturbationG.csv")
	
	xpos, xlabels = [0, 5, 10, 15, 20], ["0", r"$5r_{0}$", r"$10r_{0}$", r"$15r_{0}$", r"$20r_{0}$"]

	axs[0].plot(radii, pertK, color = 'royalblue', linestyle = '--', label = 'Kalnajs Perturbation')
	axs[0].plot(radii, pertG, color = 'firebrick', label = 'Gaussian Fit')
	axs[0].set_ylabel(r"Potential")
	axs[0].set_title(r"Kalnajs Perturbation")
	axs[0].legend(loc = 'upper right')
	axs[0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	axs[0].set_xticks(xpos)
	axs[0].set_xticklabels(xlabels)


	pertK, pertG = readingInRealCSV("BF_Comparison/gaussianPerturbationK.csv"), readingInRealCSV("BF_Comparison/gaussianPerturbationG.csv")
	axs[1].plot(radii, pertK, color = 'royalblue', label = 'Kalnajs Fit')
	axs[1].plot(radii, pertG, color = 'firebrick', linestyle = '--', label = "Gaussian Perturbation")
	axs[1].set_title(r"Gaussian Perturbation")
	axs[1].legend(loc = 'lower right')
	axs[1].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	axs[1].set_xticks(xpos)
	axs[1].set_xticklabels(xlabels)


	## Now do the density plots
	spacing = 20/(kalnajsK.nCols-1)
	centre = (kalnajsK.nCols-1)*0.5

	x = np.arange(-10,10+spacing, spacing)
	y = np.arange(-10,10+spacing, spacing)
	XX, YY = np.meshgrid(x, y)

	for i in range(len(times)):
		pcm1 = axs[2 + i*4].contourf(XX, YY, kalnajsK.densityAtTime(times[i]), levels = 100)
		pcm2 = axs[3 + i*4].contourf(XX, YY, gaussianK.densityAtTime(times[i]), levels = 100)
		pcm3 = axs[4 + i*4].contourf(XX, YY, kalnajsG.densityAtTime(times[i]), levels = 100)
		pcm4 = axs[5 + i*4].contourf(XX, YY, gaussianG.densityAtTime(times[i]), levels = 100)
		axs[3 + i*4].tick_params(labelleft=False) 
		axs[4 + i*4].tick_params(labelleft=False) 
		axs[5 + i*4].tick_params(labelleft=False)
		
		if i != (len(times)-1):  
			axs[2 + i*4].tick_params(labelbottom=False)
			axs[3 + i*4].tick_params(labelbottom=False, labelleft=False) 
			axs[4 + i*4].tick_params(labelbottom=False, labelleft=False) 
			axs[5 + i*4].tick_params(labelbottom=False, labelleft=False)

	
	
	cb_ax = fig.add_axes([0.13,.04,.36,.02])
	cbar = fig.colorbar(pcm1, cax=cb_ax, orientation='horizontal', format=ticker.FuncFormatter(fmt))# fig.colorbar(pcm1, ax = [axs[-4], axs[-3]], orientation="horizontal", format=ticker.FuncFormatter(fmt))
	cbar.set_ticks([-2*(10**-4), 0, 2*(10**-4)])
	
	cb_ax = fig.add_axes([0.53,.04,.36,.02])
	cbar = fig.colorbar(pcm3, cax = cb_ax, orientation="horizontal", format=ticker.FuncFormatter(fmt))
	cbar.set_ticks([-1*(10**-2), 0, 1*(10**-2)])

	
	axs[2].set_ylabel(r"$5.0T_{dyn}$")
	axs[6].set_ylabel(r"$7.5T_{dyn}$")
	axs[10].set_ylabel(r"$10.0T_{dyn}$")



	axs[2].set_title("Kalnajs")
	axs[3].set_title("Gaussian")
	axs[4].set_title("Kalnajs")
	axs[5].set_title("Gaussian")
	
	plt.show()



def densityComparison1D(filenames):
	kalnajsK, gaussianK, kalnajsG, gaussianG = twoDdensity(filenames[0]), twoDdensity(filenames[1]), twoDdensity(filenames[2]), twoDdensity(filenames[3])
	times = [100, 125,150]

	fig, axs = plt.subplots(nrows = 3, ncols = 2, sharex = True)
	radius = np.linspace(0, 10, len(kalnajsK.densityCutThrough(times[0])))

	for i in range(len(times)):
		axs[i,0].plot(radius, kalnajsK.densityCutThrough(times[i]), label = "Kalnajs BF")
		axs[i,0].plot(radius, gaussianK.densityCutThrough(times[i]), label = "Gaussian BF")

		axs[i, 1].plot(radius, kalnajsG.densityCutThrough(times[i]), label = "Kalnajs BF")
		axs[i, 1].plot(radius, gaussianG.densityCutThrough(times[i]), label = "Gaussian BF")

	axs[-1, 0].legend()
	axs[1, 0].set_ylabel("Density")

	axs[0,0].set_title("Global Perturbation")
	axs[0,1].set_title("Local Perturbation")


	plt.show()


def fittingStraightline(x, y):
	logY = np.log10(np.absolute(y))
	m, c = linregress(x, logY)[0:2]
	'''plt.plot(x, x*m + c)
	plt.show()'''
	print(m, c)
	return 10**(x*m + c)
	
	


def varyingN(filename):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	data = readingInRealCSV(filename)
	timeDyn = np.linspace(0, 80, np.shape(data)[1])

	fig, axs = plt.subplots(nrows = 1, ncols = 1)

	#for i in range(np.shape(data)[0]):
	colors, lst, startFit  = ["midnightblue", "#1f77b4", "cornflowerblue"], [0,2, 3], [820, 1072, 1144]
	for i in range(len(lst)):
		axs.plot(timeDyn[600:1330], 0.1*np.absolute(data[lst[i],600:1330]), label = r"$n_{upper} = $ " + str(int(data[lst[i],0])), color = colors[i])
		axs.plot(timeDyn[startFit[i]:1330], 0.1*fittingStraightline(timeDyn[startFit[i]:1330], np.absolute(data[lst[i],startFit[i]:1330])), color = colors[i], linestyle = '--')

	



	axs.legend()

	


	axs.set_ylabel(r"Energy")
	axs.set_yscale('log')
	
	plt.show()

def comparingDifferentN(filename):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	data = np.absolute(readingInRealCSV(filename))
	fig, axs = plt.subplots(nrows=1, ncols=1)
	timeDyn = np.linspace(0, 80, np.shape(data)[1]-1)
	

	axs.plot(timeDyn, data[0,1:], label = r"$N_{max}=1$ ")
	axs.plot(timeDyn, data[1,1:], label = r"$N_{max}=2$ ")

	index, label = [5,3], ["3","5"]

	for i in range(len(index)):
		axs.plot(timeDyn, data[index[i],1:], label = r"$N_{max}=$ "+str(label[i]))


	'''colors, lst, startFit, label  = ["midnightblue","cornflowerblue", "#1f77b4", ], [0,2, 3], [1500, 500, 800], [9, 7,8]
	for i in [8,9,7]:
		axs.plot(timeDyn, data[i,1:], label = r"$N_{max}=$ "+str(label[i-7]), color = colors[i-7])

		start = startFit[i-7]
		#print(fittingStraightline(timeDyn[start:-1], data[lst[i-7],start:1999]))
		axs.plot(timeDyn[start:-1], fittingStraightline(timeDyn[start:-1], data[i,start:1999]), color = colors[i-7], linestyle = '--')
	'''
	
	axs.set_yscale('log')
	axs.legend(fontsize = 14)

	axs.set_ylabel(r"Potential Energy", fontsize =14)
	axs.set_xticks(range(0, 81, 20))
	axs.set_xticklabels([0]+[str(time) + r"$T_{dyn}$" for time in range(20, 81, 20)])
	axs.set_xlabel("Time", fontsize = 14)
	plt.show()

def comparingDifferentNG(filename):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	data = np.absolute(readingInRealCSV(filename))
	fig, axs = plt.subplots(nrows=1, ncols=1)

	labels = [r"$N_{max}=$ "+str(i) for i in range(30, 51, 10)]
	timeDyn = np.linspace(0, 80, np.shape(data)[1]-1)
	startFit = [400,750, 1200]
	color, colorFit = ["cornflowerblue", "#1f77b4", "midnightblue"], ["firebrick", "indianred","lightcoral"]

	for i in range(3):
		axs.plot(timeDyn[125:1749], data[i,126:1750], label = labels[i], color = color[i])
		axs.plot(timeDyn[startFit[i]:1749], fittingStraightline(timeDyn[startFit[i]:1749], data[i,startFit[i]:1749]), color = 'firebrick', linestyle='--')

	axs.set_yscale("log")
	axs.legend(fontsize=14)
	axs.set_ylabel(r"Potential Energy", fontsize =14)

	axs.set_xlabel("Time", fontsize = 14)
	axs.set_xticks([5, 20, 35, 50, 65])
	axs.set_xticklabels([str(time) + r"$T_{dyn}$" for time in range(5, 71, 15)])
	plt.show()

def varyingPokeN(filename):
	#plt.rc('text', usetex=True)
	#plt.rc('font', family='serif')
	data = readingInRealCSV(filename)
	timeDyn = np.linspace(0, 80, np.shape(data)[1])

	fig, axs = plt.subplots(nrows = 1, ncols = 1)

	#for i in range(np.shape(data)[0]):
	for i in range(6, np.shape(data)[0]):
		axs.plot(timeDyn[20:1330], np.absolute(data[i,20:1330]), label = str(i))
		#axs.plot(timeDyn[startFit[i]:1330], fittingStraightline(timeDyn[startFit[i]:1330], np.absolute(data[lst[i],startFit[i]:1330])), linestyle = '--')
	
	
	plt.show()


def powerspectrum(nvalues):
	fig, axs = plt.subplots(nrows = 1, ncols = len(nvalues))

	for i in range(len(nvalues)):
		data = readingInComplexCSV("BF_Comparison/Coefficent_" + str(nvalues[i])+".csv")
		data = np.abs(data)


		for j in range(0, np.shape(data)[1]):
			axs[i].plot(data[50:, j], label = str(j))

		axs[i].set_yscale("log")
		
		axs[i].legend()

	plt.show()

def varyingCoupling(filename):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(nrows = 1, ncols = 1)
	

	data = np.absolute(readingInRealCSV(filename))
	time = np.linspace(0, 40, len(data[0,1:1000]))
	colors, lst, startFit  = ["cornflowerblue", "#1f77b4", "midnightblue"], [0,2, 3], [820, 1072, 1144]

	for i in range(6,6+3):
		axs.plot(time, data[i,1:1000], label = r"$N_{cut}=$ "+str(int(data[i,0])), color = colors[i-6])
		if i ==6:
			axs.plot(time[300:1000], fittingStraightline(time[300:999], data[i,300:999]), color = colors[i-6], linestyle = '--')


		#fittingStraightline(time[500:1000], data[lst[i-6],500:999])
		#axs.plot(time[500:1000], fittingStraightline(time[500:1000], data[lst[i-6],500:999]), color = colors[i-6], linestyle = '--')

	axs.legend(fontsize = 14)
	axs.set_yscale("log")	

	axs.set_xlabel(r"Time", fontsize = 14)
	axs.set_ylabel(r"Potential Energy", fontsize = 14)

	axs.set_xticks(range(0, 41, 10))
	axs.set_xticklabels([0]+[str(time) + r"$T_{dyn}$" for time in range(10, 41, 10)])

	plt.show()

def differentTimeStep():
	data = readingInRealCSV("KalnajsTesting/VaryingTimestep.csv")
	fig, axs = plt.subplots(nrows=1, ncols=1)
	for i in range(np.shape(data)[0]):
		
		time = np.linspace(0, 50, len(data[i])-1)
		axs.plot(time, -np.asarray(data[i][1:]), label = str(data[i][0])+r"$t_{0}$")


	axs.legend()
	axs.set_yscale('log')

	axs.set_xlabel("Time")
	axs.set_ylabel("Potential Energy")
	plt.show()




#comparingDifferentN("KalnajsTesting/kalnajsVaryingNbasis.csv")
#varyingCoupling("KalnajsTesting/kalnajsVaryingCoupling.csv")
#powerspectrum([5,10])
#varyingPokeN("BF_Comparison/kalnajsVaryingNpoke.csv")
#bfComparison(["BF_Comparison/kalnajsPertKalnajsRepDensity.csv", "BF_Comparison/kalnajsPertGaussianRepDensity.csv"])
#bfComparison(["BF_Comparison/gaussianPertKalnajsRepDensity.csv", "BF_Comparison/gaussianPertGaussianRepDensity.csv"]) 
#plottingPerturbation()
#densityComparison1D(["BF_Comparison/kalnajsPertKalnajsRepDensity.csv", "BF_Comparison/kalnajsPertGaussianRepDensity.csv", "BF_Comparison/gaussianPertKalnajsRepDensity.csv", "BF_Comparison/gaussianPertGaussianRepDensity.csv"])
#densityComparison(["BF_Comparison/kalnajsPertKalnajsRepDensity.csv", "BF_Comparison/kalnajsPertGaussianRepDensity.csv", "BF_Comparison/gaussianPertKalnajsRepDensity.csv", "BF_Comparison/gaussianPertGaussianRepDensity.csv"])

#powerspectrum([9,10])


#plottingPerturbation()
comparingDifferentNG("GaussianTesting/gaussianVaryingNbasis.csv")
#plottingFittedGaussian()

