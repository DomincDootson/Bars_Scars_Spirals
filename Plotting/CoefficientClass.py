from generalPlottingFunctions import *
import matplotlib.animation as animation


class CoefficientClass(object):

	def __init__(self, filename):
	
		self.coeff = readingInComplexCSV(filename)
		self.nTime, self.nMax = np.shape(self.coeff)

	def modePower(self):
		return np.absolute(self.coeff)

	def totalModePower(self):
		return self.modePower().sum(axis = 1)

	def powerPlot(self, timeEnd):
		time = np.linspace(0, timeEnd, self.nTime)

		plt.plot(time, np.sum(self.modePower(), axis = 1))
		plt.yscale('log')
		plt.show()

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



'''
coeff = CoefficientClass("Waves_Data/Coefficent_SelfConsistent_10.out")
coeff.powerEvolutions(index2Normalise = 10)

coeff = CoefficientClass("Waves_Data/Coefficent_SelfConsistent_20.out")
coeff.powerEvolutions("Waves_Plots/Waves_Videos/SelfConsistent_20.mp4", 20)

coeff = CoefficientClass("Waves_Data/Coefficent_SelfConsistent_30.out")
coeff.powerEvolutions("Waves_Plots/Waves_Videos/SelfConsistent_30.mp4", 30)'''