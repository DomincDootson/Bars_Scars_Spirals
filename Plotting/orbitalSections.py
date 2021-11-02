from generalPlottingFunctions import *
import matplotlib.animation as animation

def updateResontant(data, i):
	phi, Lz = [], []
	for j in range(np.shape(data)[0]):
		if data[j,2*i]/pi > 0:
			phi.append(data[j,2*i]/pi-1)
			Lz.append(data[j,2*i+1])

		elif data[j,2*i]/pi < 0:
			phi.append(data[j,2*i]/pi+1)
			Lz.append(data[j,2*i+1])

	return phi, Lz


def updateLibrating(data, i):
	pass

def sectionPlot(filename, nPlot = 6):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(1, 1, sharey = False, sharex = False)	

	nplot = [0,1,2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14]



	data = readingInRealCSV(filename)
	for i in nplot:
		if i <= 2:
			axs.scatter(data[:,2*i]/pi, data[:,2*i+1], s = .1, color = 'navy')

		elif i <=13:
			axs.scatter(data[:,2*i]/pi, data[:,2*i+1], s = .1, color = 'firebrick')
			phi, Lz = updateResontant(data, i)
			axs.scatter(phi, Lz, s=.1, color='firebrick')
		else:
			axs.scatter(data[:,2*i]/pi, data[:,2*i+1], s = .1, color = 'navy')

	axs.set_ylim([1.8, 2.1])	


	axs.set_xlabel(r"$\phi_{apo}/\pi$")
	axs.set_ylabel(r"$L_{z}$ $[r_{0}v_{c}]$")

	plt.show()




sectionPlot("Orbit_Sections/ResonantSections.csv")
