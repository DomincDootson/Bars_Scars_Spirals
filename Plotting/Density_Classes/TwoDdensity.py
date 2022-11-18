import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from math import *

import sys
sys.path.append('/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/Plotting')
from generalPlottingFunctions import readingInRealCSV

def findMaxima(array):
	positions = []
	for i in range(1, np.shape(array)[0]-1):
		for j in range(1, np.shape(array)[1]-1):
			cc, uc, cl, cr, lc = (array[i,j]), (array[i-1, j]), (array[i, j-1]),(array[i, j+1]), (array[i+1, j])
			d1, d2, d3, d4 =  (array[i+1, j+1]), (array[i+1, j-1]),(array[i-1, j+1]), (array[i-1, j-1])
			if ((cc > uc) and (cc > lc) and (cc > cr) and (cc > cl) and (cc > d1) and (cc > d2) and (cc > d3) and (cc > d4)):
				positions.append([i,j]) 

	return positions



class TwoDdensity(object):
	
	def __init__(self, filename):	
		flatternedDensity = list(readingInRealCSV(filename))

		self.nRows, self.nCols = int(sqrt(np.size(flatternedDensity[0]))), int(sqrt(np.size(flatternedDensity[0]))) # Assume square
		self.maxValue, self.minValue = np.amax(flatternedDensity), np.amin(flatternedDensity)
		size = self.nRows*self.nCols
		self.density2D = [np.reshape(array[:size], (self.nRows, self.nCols,)) for array in flatternedDensity] 	
		self.nSteps = len(self.density2D) 
		self.centre = (self.nCols-1)*0.5

	def __getitem__(self, index):
		return self.density2D[index]

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



	def densityAnimation(self, rMax=10, filename = None, lines = []):
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
			#contourFilled = axs.imshow(self.densityAtTime(time), vmax = maxRho, vmin = minRho, extent = (-rMax,rMax,-rMax,rMax,))
			contourFilled = axs.imshow(self.densityAtTime(time), extent = (-rMax,rMax,-rMax,rMax,))
			#point = axs.scatter(R * cos(time*0.25/R), R * sin(time*0.25/R), color = 'firebrick')

			title = fig.text(.4,.9,(r"Time: " +str(0.01*time) + r"$T_{dyn}$"))
			#ims.append([contourFilled, title, point])
			ims.append([contourFilled, title])


		for line in lines:
			axs.plot([line * cos(th) for th in np.linspace(0, 2*pi)], [line * sin(th) for th in np.linspace(0, 2*pi)], linestyle = '--', color = "firebrick")

		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		cbar = plt.colorbar(contourFilled, cax=cbar_ax)


		ani = animation.ArtistAnimation(fig, ims, interval=30)
		if (filename):
			print("Animation saved to: " + filename)
			ani.save(filename, writer = writer)
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
			#contourFilled = axs[i].imshow(self.densityAtTime(lst[i]), vmax = maxRho, vmin = minRho, extent = (-rMax,rMax,-rMax,rMax,))
			contourFilled = axs[i].imshow(self.densityAtTime(lst[i]), extent = (-rMax,rMax,-rMax,rMax,))
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

	def densityPolarPlots(self, lst, rMax = 10): 
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')

		fig, axs = plt.subplots(1,len(lst))
		ims = []

		spacing = 20/(self.nCols-1)
		centre = (self.nCols-1)*0.5

		R, phi = np.linspace(0,10, 100), np.linspace(0, 2*3.14,100)

		maxRho, minRho = self.maxDensity(), self.minDensity() 
		for i in range(len(lst)):
			#contours  = axs.contour(XX, YY, self.densityAtTime(time), 10, colors = 'k', animated = True)
			#contourFilled = axs[i].imshow(self.densityAtTime(lst[i]), vmax = maxRho, vmin = minRho, extent = (-rMax,rMax,-rMax,rMax,))
			contourFilled = axs[i].imshow(self.densityPolar2D(lst[i], R, phi, rMax), extent = (-rMax,rMax,-rMax,rMax,))
			title = axs[i].set_title(r"Time: " +str(lst[i]) + r"$T_{dyn}$")
			

		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		cbar = plt.colorbar(contourFilled, cax=cbar_ax)

		plt.show()

	def densityAnimationPolar(self, rMax=10, filename = None):
		#plt.rc('text', usetex=True)
		#plt.rc('font', family='serif')
		
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs = plt.subplots(1,1)
		ims = []

		
		R, phi = np.linspace(0,rMax, 100), np.linspace(0, 2*3.14,100)
		

		#maxRho, minRho, R = self.maxDensity(), self.minDensity(), self.maxDensityAtTime(0,rMax)
		for time in range(len(self.density2D)):
			#contourFilled = axs.imshow(self.densityAtTime(time), vmax = maxRho, vmin = minRho, extent = (-rMax,rMax,-rMax,rMax,))
			contourFilled = axs.imshow(self.densityPolar2D(time, R, phi, rMax), extent = (0,2 * 3.14,rMax,0))

			title = fig.text(.4,.9,(r"Time: " +str(0.5* time)))
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



	def fourierCoeffT(self, angHarmonic, time, rMax = 15): # We fit rho = A(R)Cos(l phi + psi(R)) 
		
		phase, magnitude = [], []
		radii, angles  = np.linspace(0.01, rMax, 100), np.linspace(0, 2 * 3.14, 100)

		dtheta = angles[1] - angles[0]

		for radius in radii:

			sinI, cosI = 0, 0
			for theta in angles:
				sinI += sin(angHarmonic * theta) * (dtheta / pi) * self.densityPolar(time, radius, theta, radii[-1]) # Factor of pi 
				cosI += cos(angHarmonic * theta) * (dtheta / pi) * self.densityPolar(time, radius, theta, radii[-1]) # to normalise the integrals 
		
			phase.append(atan(-sinI/cosI))
			if isnan(phase[-1]): # For the case when the amplitude is zero everywhere at the start
				phase[-1] = 0
			magnitude.append(sqrt(sinI**2  + cosI**2))

		return radii, magnitude, phase

	def fourierAnimations(self, angHarmonic, filename = None, rMax = 15):
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs1 = plt.subplots(1,1)
		axs2 = axs1.twinx()
		ims = []

		for time in range(self.nSteps):
			radii, mag, ph = self.fourierCoeffT(angHarmonic, time, rMax)

			magnitude, = axs1.plot(radii, mag, color = 'royalblue', label = "Magnitude") 

			#phase, = axs2.plot(radii, ph, color = 'firebrick', label = "Magnitude") 
			title = fig.text(.4,.9,(r"Time: " +str(0.5*time)))
			
			ims.append([magnitude, title])


		ani = animation.ArtistAnimation(fig, ims, interval=30)
		if (filename):
			ani.save(filename, writer = writer)
			print("Animation save to:" + filename)
		else:
			plt.show()

		
	def fouierEvolution(self, angHarmonic, rMax = 10):
		lst = [[*self.fourierCoeffT(angHarmonic, timeIndex, rMax)] for timeIndex in range(self.nSteps)]

		return np.asarray(lst[0][0]), np.asarray([each[1] for each in lst]), np.asarray([each[2] for each in lst])

	def maximaEvolution(self):
		positions = []

		for time in range(self.nSteps):
			positions.append(findMaxima(self.densityAtTime(time)- self.densityAtTime(0)))

		radii, time = [], []
		for i, timeSet in enumerate(positions):
			for each in timeSet:
				radii.append(sqrt((each[0]-self.centre)**2 + (each[1]-self.centre)**2))
				time.append(i)
				

		return time, radii



