import numpy as np
import matplotlib.pyplot as plt
import csv
#from cmath import *
from math import *
import matplotlib.tri as tri # this allows for contour plot with non-uniform spacing 
from matplotlib.colors import * 

def readingInRealCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [float(i) for i in row]
			data.append(lst)

		return np.asarray(data)

def readingInComplexCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [complex(i) for i in row]
			#lst = [print(row) for i in row] 
			
			data.append(lst)

		return np.asarray(data)

def readingInRealOUT(filename):
	data = []
	with open(filename, 'r') as file: 
		#data.append((file.readline()).split())
		for line in file:
			lst = [float(number) for number in line.split()]
			data.append(lst)

		file.close()
	del data[0] 
	return np.asarray(data)



ENERGY_DIR = "../Disk_Kicking/Energy_Evolution/"

class EnergyEvolutionData(): # This holds the data output by the C++ code
	
	def __init__(self, filename, timeStep = 0.25):
		self.m_data = readingInRealCSV(ENERGY_DIR + filename)
		self.m_timeStep = timeStep

	def __init__(self, filename, rInner, rOuter, timeStep = 0.25):
		self.m_data = readingInRealCSV(ENERGY_DIR + filename)
		self.m_rInner = rInner
		self.m_rOuter = rOuter
		self.m_timeStep = timeStep


	def check_Agreement(self, row, littleSigma, angHarmonic, radius):
		if (self.m_data[row, 0] == littleSigma) & (self.m_data[row, 1] == angHarmonic) & (self.m_data[row,2] == radius):
			return True
		else:
			return False

	def energy_evolution(self, littleSigma, angHarmonic, radius):
		for i in range(np.shape(self.m_data)[0]):
			if self.check_Agreement(i, littleSigma, angHarmonic, radius):
				return self.m_data[i, 3:]
		
		print("Incorrect values.")
		print("Values entered: ", littleSigma, angHarmonic, radius)
		exit(0)


	def little_sigma(self):
		return np.unique(self.m_data[:,0])

	def ang_harmonic(self):
		return np.unique(self.m_data[:,1])

	def radii(self):
		return np.unique(self.m_data[:,2])


	def max_Energy(self, littleSigma, angHarmonic):
		radii = self.radii()
		lst = []
		for radius in radii:
			lst.append(np.amin(self.energy_evolution(littleSigma, angHarmonic, radius)))
		return np.asarray(lst)

	def numb_Time_Steps(self):
		return np.shape(self.m_data[0,3:])[0]

	def index_Max_Energy(self, littleSigma, angHarmonic, radius):
		maxValue = np.amin(self.energy_evolution(littleSigma, angHarmonic, radius))
		energy = self.energy_evolution(littleSigma, angHarmonic, radius)
		for i in range(self.numb_Time_Steps()):
			if energy[i] == maxValue:
				return i

	def time_Max_Energy(self, littleSigma, angHarmonic, radius):
		maxIndex = self.index_Max_Energy(littleSigma, angHarmonic, radius)
		energy = self.energy_evolution(littleSigma, angHarmonic, radius)
		if energy[maxIndex] == energy[-1]:
			return self.m_timeStep * maxIndex

		gradBefore, gradAfter = abs((energy[maxIndex]- energy[maxIndex-1])), abs((energy[maxIndex+1]- energy[maxIndex]))
		return self.m_timeStep *( (maxIndex-.5) +  gradBefore/(gradAfter+gradBefore)) 

	def time_Max_Energy_All_Radii(self, littleSigma, angHarmonic):
		return [self.time_Max_Energy(littleSigma, angHarmonic, radius) for radius in self.radii()]

class DFClass(object):
	"""docstring for DFClass"""
	def __init__(self, filename):
		data = readingInRealCSV(filename)

		self.coords = [data[0,:], data[1,:]]
		self.dfFunction = data[2:,:]

	def dfAtTime(self, timeIndex):
		print (self.dfFunction[timeIndex, -1], np.amax(self.dfFunction[timeIndex, :]))
		return self.dfFunction[timeIndex, :]

	def apoCentrePlots(self, timeIndex):
		fig, axs = plt.subplots(1, 1)

		axs.tricontour(self.coords[0], self.coords[1], self.dfAtTime(timeIndex), colors = 'black')
		cntr = axs.tricontourf(self.coords[0], self.coords[1], self.dfAtTime(timeIndex), levels=50)
		fig.colorbar(cntr, ax=axs)

		axs.set(xlim=(0.2, 7), ylim=(0.4, 7))

		axs.set_xlabel(r"$r_{-}$")
		axs.set_ylabel(r"$r_{+}$")

		plt.show()

	def apoDFEvolution(self):
		fig, axs = plt.subplots(1, np.shape(self.dfFunction)[0])
		maxValue, minValue = np.amax(self.dfFunction), np.min(self.dfFunction)
		print(minValue, maxValue)

		for i in range(np.shape(self.dfFunction)[0]):
			axs[i].tricontour(self.coords[0], self.coords[1], self.dfAtTime(i), colors = 'black')
			axs[i].set(xlim=(0.2, 7), ylim=(0.4, 7))
			
			cntr = axs[i].tricontourf(self.coords[0], self.coords[1], self.dfAtTime(i), levels=50, vmin = minValue, vmax = maxValue)
			cntr.set_clim(minValue, maxValue)
			axs[i].set_xlabel(r"$r_{-}$")
			axs[i].set_ylabel(r"$r_{+}$")

		
		norm = Normalize(vmin=minValue, vmax=maxValue)

		cbar = fig.colorbar(cntr, boundaries = np.linspace(minValue, maxValue, 10))
		cbar.formatter.set_powerlimits((0, 4))
		cbar.update_ticks()

		plt.show()

	def elPlots(self, timeIndex): # Why is this so weird?? 
		fig, axs = plt.subplots(1, 1)

		axs.tricontour(self.coords[0], self.coords[1], self.dfAtTime(timeIndex), colors = 'black')
		cntr = axs.tricontourf(self.coords[0], self.coords[1], self.dfAtTime(timeIndex), levels=50)
		fig.colorbar(cntr, ax=axs)

		axs.set(xlim=(2, 3.5))
		axs.set_xlabel(r"$E$")
		axs.set_ylabel(r"$L$")

		plt.show()



class twoDdensity(object):
	
	def __init__(self, filename):	
		flatternedDensity = list(readingInRealCSV(filename))
		self.nRows, self.nCols = int(sqrt(np.size(flatternedDensity[0]))), int(sqrt(np.size(flatternedDensity[0]))) # Assume square

		self.maxValue, self.minValue = np.amax(flatternedDensity), np.amin(flatternedDensity)
		self.density2D = [np.reshape(array[:self.nRows*self.nCols], (self.nRows, self.nCols,)) for array in flatternedDensity] 	

	def densityAtTime(self, timeIndex):
		return self.density2D[timeIndex]

	def densityEvolution(self, rMax):

		spacing = 20/(self.nCols-1)
		centre = (self.nCols-1)*0.5

		x = np.arange(-10,10+spacing, spacing)
		y = np.arange(-10,10+spacing, spacing)
		XX, YY = np.meshgrid(x, y)

		fig, axs = plt.subplots(1, len(self.density2D))

		for i in range(len(axs)):
			#axs[i].contour(XX, YY, self.densityAtTime(i), 4, colors = 'k')
			cntr = axs[i].contourf(XX, YY, self.densityAtTime(i), levels = 50)
			cntr.set_clim(self.minValue, self.maxValue)

		plt.show()





def maxDensity():
	data = readingInRealCSV("../Disk_Kicking/maxDensity.csv")
	rInner = data[1:,0]
	print(rInner)
	for s in range(1, np.shape(data)[1]):
		plt.plot(rInner, data[1:, s]-rInner, label = r"$\sigma_{r}=$ " + str(data[0,s]))

	plt.ylabel(r"$R_{Max Density}- R_{i}$")
	plt.xlabel(r"$R_{i}$")

	plt.legend()
	plt.show()