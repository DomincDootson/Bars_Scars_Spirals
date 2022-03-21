import numpy as np
import matplotlib.pyplot as plt
import csv
#from cmath import *
from math import *
import matplotlib.tri as tri # this allows for contour plot with non-uniform spacing 
from matplotlib.colors import * 
import matplotlib.animation as animation

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



class twoDdensity(object):
	
	def __init__(self, filename):	
		flatternedDensity = list(readingInRealCSV(filename))
		self.nRows, self.nCols = int(sqrt(np.size(flatternedDensity[0]))), int(sqrt(np.size(flatternedDensity[0]))) # Assume square

		self.maxValue, self.minValue = np.amax(flatternedDensity), np.amin(flatternedDensity)
		self.density2D = [np.reshape(array[:self.nRows*self.nCols], (self.nRows, self.nCols,)) for array in flatternedDensity] 	
		self.nSteps = len(self.density2D)

	def densityAtTime(self, timeIndex):
		return self.density2D[timeIndex]

	def densityCutThrough(self, timeIndex):
		return self.density2D[timeIndex][round(self.nRows*0.5),round(self.nRows*0.5):]

	def maxDensityAtTime(self, timeIndex, rMax = 10):
		spacing = 2*rMax/(self.nCols-1)
		centre = (self.nCols-1)*0.5

		ind = np.unravel_index(np.argmax(self.densityAtTime(timeIndex), axis=None), self.densityAtTime(timeIndex).shape)

		return spacing * sqrt((ind[0]-centre)**2+(ind[1]-centre)**2)

	def maxDensityRadius(self, rMax):
		lst = []
		for i in range(self.nSteps):
			lst.append(self.maxDensityAtTime(i, rMax))

		return lst

	def maxDensityEvolution(self):
		return [np.amax(each) for each in self.density2D]


	def densityEvolution(self, rMax):

		spacing = 2*rMax/(self.nCols-1)
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

	def maxDensityTime(self, time):
		return np.amax(self.densityAtTime(time))
	def maxDensity(self): 
		return max([np.amax(each) for each in self.density2D])

	def minDensity(self):
		return min([np.amin(each) for each in self.density2D])

	def densityPower(self, timeIndex):
		return np.sum(np.square(self.densityAtTime(timeIndex)))

	def densityPowerEvolution(self):
		lst = []
		for i in range(self.nSteps):
			lst.append(self.densityPower(i))
		
		return 	[i/self.densityPower(0) for i in lst]



	def densityAnimation(self, rMax=10, filename = None):
		#plt.rc('text', usetex=True)
		#plt.rc('font', family='serif')
		
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs = plt.subplots(1,1)
		ims = []

		spacing = 20/(self.nCols-1)
		centre = (self.nCols-1)*0.5

		x = np.arange(-10,10+spacing, spacing)
		y = np.arange(-10,10+spacing, spacing)
		XX, YY = np.meshgrid(x, y)

		maxRho, minRho, R = self.maxDensity(), self.minDensity(), self.maxDensityAtTime(0,rMax)
		for time in range(len(self.density2D)):
			contourFilled = axs.imshow(self.densityAtTime(time), vmax = maxRho, vmin = minRho, extent = (-rMax,rMax,-rMax,rMax,))
			#point = axs.scatter(R * cos(time*0.25/R), R * sin(time*0.25/R), color = 'firebrick')

			title = fig.text(.4,.9,(r"Time: " +str(0.01*time) + r"$T_{dyn}$"))
			#ims.append([contourFilled, title, point])
			ims.append([contourFilled, title])


		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		cbar = plt.colorbar(contourFilled, cax=cbar_ax)


		ani = animation.ArtistAnimation(fig, ims, interval=30)
		if (filename):
			ani.save(filename, writer = writer)
			print("Animation save to: " + filename)
		else:
			plt.show()

	def densityAnimation1D(self, rMax=10, filename = None):
		#plt.rc('text', usetex=True)
		#plt.rc('font', family='serif')
		
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs = plt.subplots(1,1)
		ims = []

		spacing = 20/(self.nCols-1)
		centre = (self.nCols-1)*0.5

		x = np.arange(-10,10+spacing, spacing)
		y = np.arange(-10,10+spacing, spacing)
		XX, YY = np.meshgrid(x, y)

		maxRho, minRho, R = self.maxDensity(), self.minDensity(), self.maxDensityAtTime(0,rMax)
		for time in range(len(self.density2D)):
			contourFilled, = axs.plot(self.densityCutThrough(time), color = 'royalblue')
			#point = axs.scatter(R * cos(time*0.25/R), R * sin(time*0.25/R), color = 'firebrick')

			title = fig.text(.4,.9,(r"Time: " +str(0.01*time) + r"$T_{dyn}$"))
			#ims.append([contourFilled, title, point])
			ims.append([contourFilled, title])


		ani = animation.ArtistAnimation(fig, ims, interval=30)
		if (filename):
			ani.save(filename, writer = writer)
			print("Animation save to:" + filename)
		else:
			plt.show()

	def densityPlots(self, lst, rMax = 10): 
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')

		fig, axs = plt.subplots(1,len(lst))
		ims = []

		spacing = 20/(self.nCols-1)
		centre = (self.nCols-1)*0.5

		x = np.arange(-10,10+spacing, spacing)
		y = np.arange(-10,10+spacing, spacing)
		XX, YY = np.meshgrid(x, y)

		maxRho, minRho = self.maxDensity(), self.minDensity() 
		for i in range(len(lst)):
			#contours  = axs.contour(XX, YY, self.densityAtTime(time), 10, colors = 'k', animated = True)
			contourFilled = axs[i].imshow(self.densityAtTime(lst[i]), vmax = maxRho, vmin = minRho, extent = (-rMax,rMax,-rMax,rMax,))
			title = axs[i].set_title(r"Time: " +str(lst[i]) + r"$T_{dyn}$")
			

		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		cbar = plt.colorbar(contourFilled, cax=cbar_ax)

		plt.show()


	def densityPolar(self, time, R, phi, rMax):
		spacing, centre = 2*rMax/(self.nCols-1), (self.nCols-1)*0.5

		iIndex = lambda R, phi : (R * sin(phi))/spacing + centre
		jIndex = lambda R, phi : (R * cos(phi))/spacing + centre

		i0, i, i1 = floor(iIndex(R, phi)), iIndex(R, phi), ceil(iIndex(R, phi))
		j0, j, j1 = floor(jIndex(R, phi)), jIndex(R, phi), ceil(jIndex(R, phi))

		if (i0==i1 & j0==j1):
			return self.densityAtTime(time)[i0,j0]

		elif (i0==i1):
			return self.densityAtTime(time)[i0,j0] * (j1-j) + self.densityAtTime(time)[i0,j1] * (j-j0)

		elif (j0==j1):
			return self.densityAtTime(time)[i0,j0] * (i1-i) + self.densityAtTime(time)[i1,j0] * (i-i0)

		else:
			return self.densityAtTime(time)[i0,j0] * (i1-i)*(j1-j) + self.densityAtTime(time)[i1,j0] * (i-i0)*(j1-j) + self.densityAtTime(time)[i0,j1] * (i1-i)*(j-j0) + self.densityAtTime(time)[i1,j1] * (i-i0)*(j-j0)			


	def densityPolar2D(self, time, R, phi, rMax): # note that R and phi are lists 
		density = np.zeros((np.shape(R)[0], np.shape(phi)[0]))

		for i in range(np.shape(R)[0]):
			for j in range(np.shape(phi)[0]): 
				density[i,j] = self.densityPolar(time, R[i], phi[j], rMax)

		return density

	def densityAnimationPolar(self, rMax=10, filename = None):
		#plt.rc('text', usetex=True)
		#plt.rc('font', family='serif')
		
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs = plt.subplots(1,1)
		ims = []

		
		R, phi = np.linspace(0,10, 100), np.linspace(0, 2*3.14,100)
		

		#maxRho, minRho, R = self.maxDensity(), self.minDensity(), self.maxDensityAtTime(0,rMax)
		for time in range(len(self.density2D)):
			#contourFilled = axs.imshow(self.densityAtTime(time), vmax = maxRho, vmin = minRho, extent = (-rMax,rMax,-rMax,rMax,))
			contourFilled = axs.imshow(self.densityPolar2D(time, R, phi, rMax), extent = (0,2 * 3.14,0,rMax,))

			title = fig.text(.4,.9,(r"Time: " +str(time) + r"$T_{dyn}$"))
			ims.append([contourFilled, title])


		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		cbar = plt.colorbar(contourFilled, cax=cbar_ax)


		ani = animation.ArtistAnimation(fig, ims, interval=30)
		if (filename):
			ani.save(filename, writer = writer)
			print("Animation save to:" + filename)
		else:
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