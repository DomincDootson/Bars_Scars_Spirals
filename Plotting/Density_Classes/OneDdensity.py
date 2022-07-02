import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import sys
sys.path.append('/Users/dominicdootson/Documents/PhD/phd/Linear_Stability_Clean/Plotting')
from generalPlottingFunctions import readingInRealCSV

class OneDdensity(object):
	def __init__(self, filename):
		holding = readingInRealCSV(filename)
		self.densityArray = holding[1:, :]
		self.radii = holding[0,:]
		
		self.nStep = np.shape(self.densityArray)[0]
		self.nRows = np.shape(self.densityArray)[1]
		

	def density(self, timeIndex):
		return self.densityArray[timeIndex, :]	

	def densityAnimation(self, write2file = None):
		
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs = plt.subplots(1,1)
		ims = []

		for time in range(self.nStep):		
			title = fig.text(.4,.9,(r"Time: " +str(0.01*time) + r"$T_{dyn}$"))
			line, = plt.plot(self.radii, self.density(time), animated = True, color = 'navy')
			ims.append([line, title]) 


		ani = animation.ArtistAnimation(fig, ims, interval=30)
		if (write2file):
			print("Animation saved to: " + write2file)
			ani.save(write2file, writer = writer)
		else:
			plt.show()