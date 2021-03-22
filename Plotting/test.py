from generalPlottingFunctions import *


def plottingImpulsePerturbations():	

	data = readingInRealCSV("../Disk_Kicking/PlottingPerturbation.csv")
	radii = np.linspace(0,10, np.shape(data)[1]-1)
	rNudge = [6, 46]
	
	
	for r in rNudge:
		plt.plot(radii, data[r, 1:])
		print(np.amax(data[r, 1:]))
		plt.plot([data[r, 0], data[r, 0]], [0, (np.amax(data[r, 1:]))], color = 'grey')
	
	plt.xlabel("Radius")
	plt.ylabel("Potential")
	plt.show()

def differentTimeSteps():
	n = 5
	data2000 = readingInComplexCSV("evolution2000.csv")
	plt.plot(np.linspace(0,50, np.shape(data2000)[0]), np.absolute(data2000[:, n]))

	data100 = readingInComplexCSV("evolution100.csv")
	plt.plot(np.linspace(0,50, np.shape(data100)[0]), np.absolute(data100[:, n]))
	
	data200 = readingInComplexCSV("evolution200.csv")
	plt.plot(np.linspace(0,50, np.shape(data200)[0]), np.absolute(data200[:, n]))
	
	plt.show()


#plottingImpulsePerturbations()
differentTimeSteps()