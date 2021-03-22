from generalPlottingFunctions import *


def coefficentFilename(littleSigma, radius, angHarmonic, testParticle = False):
	filename = "../Disk_Kicking/littleSigma_" + str(round(littleSigma*100)) + "/Coeff" + str(round(radius*10)) +"_" + str(angHarmonic) 
	if testParticle:
		filename += "_Test"
	return (filename + ".csv")


def findRoots(n, someTimeSeries):
	times = []
	for i in range(np.shape(someTimeSeries)[0]-1):
		if (someTimeSeries[i, n] * someTimeSeries[i+1, n] < 0):
			times.append(i * 50 /np.shape(someTimeSeries)[0])

	return times

def coeffEvolutionPlot(littleSigma, radius):
	modes, angHarmonics = range(6), range(3)
	data = [readingInComplexCSV(coefficentFilename(littleSigma, radius, l)) for l in angHarmonics]

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(len(modes), len(angHarmonics), sharex = True, sharey = False)
	time = np.linspace(0,50, np.shape(data[0])[0])
	for i in modes:
		for j in angHarmonics:
			axs[i,j].ticklabel_format(axis='y', style = 'sci', scilimits = (0,0))
			axs[i,j].plot(time, np.real(data[j][:, i]), label = 'Real')

			if j != 0:
				axs[i,j].plot(time, np.imag(data[j][:, i]), label = 'Imag')

	
	fig.suptitle("Evolution of Coefficents, $r_{n}=$" + str(radius))
	axs[-1, 1].set_xlabel(r"Time $t_{0}$")
	axs[-1,-1].legend()
	for n in modes:
		axs[n, 0].set_ylabel(r"$n_{p}=$ " + str(n))
	for l in angHarmonics:
		axs[0,l].set_title(r"$\ell_{p} = $ " + str(l))


	plt.show()



def rootPlotting(littleSigma):
	modes, angHarmonics, radii = range(7), range(3), range(1,6)
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(len(modes), len(angHarmonics), sharex = True, sharey = True)
	data = []
	for n in modes:
		holdingModes = []
		for l in angHarmonics:
			holdingHarmonics = []
			for r in radii:
				holdingData = readingInComplexCSV(coefficentFilename(littleSigma, r, l))
				holdingHarmonics.append(findRoots(n, np.real(holdingData)))
				
			
			holdingModes.append(holdingHarmonics)
		data.append(holdingModes)


	for i in modes:
		for j in angHarmonics:
			for r in range(len(data[0][0])):

				plotX = [r+1]*len(data[i][j][r])
				axs[i,j].scatter(plotX, data[i][j][r], s = 10)

	
	axs[-1,1].set_xlabel(r"$r_{n}$")
	axs[3, 0].set_ylabel(r"Roots of Coefficent Oscillation $t_{0}$")
	for l in angHarmonics:
		axs[0,l].set_title(r"$\ell_{p}=$ " + str(l))

	fig.suptitle(r"Time Period Of Oscillation, $\sigma_{r}=$ " +str(littleSigma))

	plt.show()






coeffEvolutionPlot(.25, 5)


