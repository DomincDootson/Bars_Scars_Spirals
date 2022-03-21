from generalPlottingFunctions import *
import matplotlib.animation as animation

def updateResonant(data, i):
	phi, Lz = [], []
	for j in range(np.shape(data)[0]):
		if data[j,2*i]/pi > 0:
			phi.append(data[j,2*i]/pi-1)
			Lz.append(data[j,2*i+1])

		elif data[j,2*i]/pi < 0:
			phi.append(data[j,2*i]/pi+1)
			Lz.append(data[j,2*i+1])

	return phi, Lz


def isTrapped(rapo): 
	rapo = rapo / pi 

	if (abs(rapo[0]) < 0.5):
		whichRegion = 0
	else:
		whichRegion = 1

	for i in range(1, np.shape(rapo)[0]):
		if (abs(rapo[i]) < 0.5):
			checkregion = 0
		else:
			checkregion = 1

		if (checkregion != whichRegion):
			return False
		whichRegion = checkregion

	return True
	 

	

def sectionPlot(filename, nPlot = 6):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(1, 1, sharey = False, sharex = False)	

	nplot = [9]
	#nplot = [0]
	#nplot = [0,2]
	#nplot = [0,2, 3, 4, 5, 7, 9, 10, 11, 12, 13]

	n = 25
	data = readingInRealCSV(filename)


	n, count = 25, 0
	
	for i in range(0,20):
		
		if isTrapped(data[n:,2*i]):
			axs.scatter(data[n:,2*i]/pi, data[n:,2*i+1], s = 1, color = 'navy')
			count += 1
		else:
			axs.scatter(data[n:,2*i]/pi, data[n:,2*i+1], s = 1, color = 'firebrick')

	
	#axs.set_ylim([1.8, 2.1])	


	axs.set_xlabel(r"$\phi_{apo}/\pi$")
	axs.set_ylabel(r"$L_{z}$ $[r_{0}v_{c}]$")
	#axs.legend()
	print(count)
	plt.show()



filenames = ["Orbit_Sections/SectionFastSmall.csv", "Orbit_Sections/SectionSlowSmall.csv", "Orbit_Sections/SectionFastLarge.csv", "Orbit_Sections/SectionSlowLarge.csv"]


sectionPlot("Orbit_Sections/ResonantSections.csv")
