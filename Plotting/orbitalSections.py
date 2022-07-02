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
	nplot = [0,2,4,6] + [i for i in range(8,28)] + [28, 30, 31]

	n = 24
	data = readingInRealCSV(filename)


	n, count = 25, 0
	
	for i in nplot:
		
		if isTrapped(data[n:,2*i]):
			axs.scatter(data[n:,2*i]/pi, data[n:,2*i+1], s = .1, color = 'navy')
			count += 1
		else:
			axs.scatter(data[n:,2*i]/pi, data[n:,2*i+1], s = .1, color = 'firebrick')


	axs.text(-.1,2.12,"Area 1", fontsize = 15)
	axs.text(-.85,2.12,"Area 2", fontsize = 15)
	axs.text(.65,2.12,"Area 2", fontsize = 15)

	
	axs.set_ylim([1.79, 2.11])	
	axs.set_xlim([-1, 1])	

	axs.axvline(-0.5, color = 'k', linestyle = '--')
	axs.axvline(0.5, color = 'k', linestyle = '--')


	axs.set_xlabel(r"$\phi_{apo}/\pi$", fontsize = 15)
	axs.set_ylabel(r"$L_{z}$", fontsize = 15)
	#axs.legend()
	print(count)
	plt.show()

def ptleChangeInAngMom():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(1, 1, sharey = False, sharex = False)	


	data = readingInRealCSV("Orbit_Sections/angmom.csv")
	xlabels, xlabelsPositions = [0] + [str(i) + r"$T_{dyn}$" for i in range(1,6)], [0] +[i for i in range(20, 101, 20)]

	time = np.linspace(0, 100, np.shape(data[:2000])[0])


	axs.plot(time, (data[:2000,6] - data[0,6])/data[0,6], color = 'navy', label = "Outer Particle")
	axs.plot(time, (data[:2000,28] - data[0,28])/data[0,28], color = 'firebrick', label = "Inner Particle")
	
	axs.set_xticks(xlabelsPositions)
	axs.set_xticklabels(xlabels)

	#axs.set_xlabel(r"Time")
	axs.set_ylabel(r"$\Delta L_{z}/L_{z}(t=0)$ $[r_{0}v_{c}]$")
	axs.legend()




	plt.show()



filenames = ["Orbit_Sections/SectionFastSmall.csv", "Orbit_Sections/SectionSlowSmall.csv", "Orbit_Sections/SectionFastLarge.csv", "Orbit_Sections/SectionSlowLarge.csv"]

#ptleChangeInAngMom()
sectionPlot("Orbit_Sections/ResonantSections.csv")
