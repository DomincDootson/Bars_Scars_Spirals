from generalPlottingFunctions import *
import matplotlib.animation as animation
from scipy.stats import linregress


class CoefficientClass(object):

	def __init__(self, filename, endTime = 1):
		self.coeff = readingInComplexCSV(filename)
		self.nTime, self.nMax = np.shape(self.coeff)
		self.time = np.linspace(0, endTime, self.nTime)

	def __getitem__(self, index):
		return self.coeff[index, time]


	def modePower(self):
		return np.square(np.absolute(self.coeff))

	def arg(self):
		return np.angle(self.coeff)

	def totalModePower(self):
		return self.modePower().sum(axis = 1)

	def powerPlot(self, timeEnd):
		time = np.linspace(0, timeEnd, self.nTime)

		plt.plot(time, np.sum(self.modePower(), axis = 1))
		plt.yscale('log')
		plt.show()

	## Mode Fitting ##
	## ------------ ##

	def expoFit(self, time):
		timeFit, logPower = time[self.nTime//2:], np.log(self.totalModePower()[self.nTime//2:])
		fit = linregress(timeFit, logPower)
		return timeFit, np.exp(fit.intercept + fit.slope * timeFit), fit.slope/2

	def periodFitter4Index(self, index, time):
		timeFit, arg = time[self.nTime//2:], self.arg()[self.nTime//2:, index]
		roots = [index for index in range(1,np.shape(arg)[0]) if (arg[index]>0 and arg[index-1]<0)]
		
		timePeriod = (timeFit[roots[-1]] - timeFit[roots[0]]) / (len(roots) - 1)

		return (2 * pi) / timePeriod

	def periodFitter(self, time):
		omega4index = [self.periodFitter4Index(i, time) for i in range(self.nMax)]
		return sum(omega4index)/len(omega4index)
		

	def powerEvolutions(self, filename = None, index2Normalise = -1):
		modePower = self.modePower()

		if index2Normalise != -1:
			modePower[:, index2Normalise] -= modePower[0,index2Normalise]
		bins, minValue, maxValue = range(self.nMax), np.amin(modePower), np.amax(modePower)
		colors = ["firebrick" if i == index2Normalise else "royalblue" for i in range(0, self.nMax)]

		fig, axs = plt.subplots()
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		def animate(time):
			axs.clear()
			axs.set_title("Time = " + str(time*0.5))
			axs.set_ylim([minValue, maxValue])
			axs.bar(bins, modePower[time, :], color = colors)
			return axs

		ani = animation.FuncAnimation(fig, animate, frames = 250)
		if (filename):
			ani.save(filename, writer = writer)
			print("Animation saved to: " + filename) 
		else:
			plt.show()

	## Fitting Modes ##
	## ------------- ##

	def fittingGrowthIndex(self, index):
		timeFit, logMagnitude = self.time[self.nTime//2:], np.log(np.absolute(self.coeff[self.nTime//2:, index]))
		fit = linregress(timeFit, logMagnitude)
		return fit.slope

	def fitGrowth(self):
		etas = [self.fittingGrowthIndex(i) for i in range(self.nMax)]
		return sum(etas)/len(etas)

	def fitFreqIndex(self, index):
		timeFit, phase = self.time[self.nTime//2:], np.angle(self.coeff[self.nTime//2:, index])
		_ , gradient = getGradient(timeFit, phase, multiple = 8, includeMinus = 1)

		return sum(gradient)/len(gradient)

	def fitFreq(self):
		omega_0s = [self.fitFreqIndex(i) for i in range(self.nMax)]
		return -sum(omega_0s)/len(omega_0s)

	def fitUnstableMode(self):
		return self.fitFreq() + 1j * self.fitGrowth()

'''
coeff = CoefficientClass("Waves_Data/Coefficent_SelfConsistent_10.out")
coeff.powerEvolutions(index2Normalise = 10)

coeff = CoefficientClass("Waves_Data/Coefficent_SelfConsistent_20.out")
coeff.powerEvolutions("Waves_Plots/Waves_Videos/SelfConsistent_20.mp4", 20)

coeff = CoefficientClass("Waves_Data/Coefficent_SelfConsistent_30.out")
coeff.powerEvolutions("Waves_Plots/Waves_Videos/SelfConsistent_30.mp4", 30)'''

