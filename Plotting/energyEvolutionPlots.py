from generalPlottingFunctions import *

ENERGY_DIR = "../Disk_Kicking/Energy_Evolution/"

class EnergyEvolutionData(object): # This hols the data output by the C++ code
	"""docstring for EnergyEvolutionData"""
	def __init__(self, filename):
		self.m_data = readingInRealCSV(ENERGY_DIR + filename)

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

	def littleSigma(self):
		return np.unique(self.m_data[:,0])

	def angHamonics(self):
		return np.unique(self.m_data[:,1])

	def radii(self):
		return np.unique(self.m_data[:,2])

test = EnergyEvolutionData("KalnajsEnergyEvolutionTest.csv")
selfConsistent = EnergyEvolutionData("KalnajsEnergyEvolution.csv")

plt.plot(test.energy_evolution(0.25, 1, 5))
plt.plot(selfConsistent.energy_evolution(0.25, 1, 5))

plt.show()