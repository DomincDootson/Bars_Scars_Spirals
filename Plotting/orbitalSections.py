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

	if (rapo[0] < 0):
		whichRegion = 0
	else:
		whichRegion = 1

	for r_a in rapo:
		if (r_a < 0.0):
			checkregion = 0
		else:
			checkregion = 1

		if (checkregion != whichRegion):
			return False


	return True
	 

	

def sectionPlot(filename = "Orbit_Sections/ResonantSections.csv", nPlot = 6):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(1, 1, sharey = False, sharex = False)	

	data = readingInRealCSV(filename)


	n, count = 10, 0
	
	for i in range(0,48,2):
		
		if isTrapped(data[n:,2*i]):
			
			axs.scatter(data[n:,2*i]/pi, data[n:,2*i+1], s = .1, color = 'navy')
			axs.scatter(-data[n:,2*i]/pi, data[n:,2*i+1], s = .1, color = 'navy')
			count += 1
		else:
			axs.scatter(data[n:,2*i]/pi, data[n:,2*i+1], s = .1, color = 'firebrick')




	
	

	axs.set_xlim([-1, 1])	

	axs.axvline(0, color = 'k', linestyle = '--')
	


	axs.set_xlabel(r"$\phi_{apo}/\pi$", fontsize = 15)
	axs.set_ylabel(r"$L_{z}$", fontsize = 15)
	#axs.legend()
	
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





#ptleChangeInAngMom()
#sectionPlot("Orbit_Sections/ResonantSections.csv")

sectionPlot("../nBody/sormaniILRSections.csv") # Bar_Data/Sections/
