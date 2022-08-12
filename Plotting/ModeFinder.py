from generalPlottingFunctions import *
from matplotlib import ticker, cm

class ModeFinder():

	def __init__(self, filename):
		with open(filename) as csv_file:
			csv_reader = csv.reader(csv_file, delimiter = ',')
			det, omega = [], []
			readingInOmega = True

			for row in csv_reader:
				lst = [complex(i) for i in row]
				
				if (len(row) == 0):
					readingInOmega = False
					
				if (readingInOmega):
					omega.append(lst)
				elif (len(row)!=0):
					det.append(lst)

				#data.append(lst)
		

		self.omegas, self.det = np.asarray(omega, dtype=np.cdouble), np.asarray(det, dtype=np.cdouble)
		self.nO, self.nE = np.shape(self.omegas)[1], np.shape(self.omegas)[0]
		print(self.nO)

	def abs(self):
		return np.absolute(self.det)

	def contourPlot(self):
		XX, YY = np.real(self.omegas), np.imag(self.omegas)
		plt.contourf(XX, YY, np.absolute(self.det), levels = 100, locator=ticker.LogLocator())
		plt.colorbar()

		minVal = self.minValue()
		plt.scatter(minVal.real, minVal.imag)

		plt.show()

	def omegaPlot(self):
		fig, axs = plt.subplots()
		for i in range(self.nE):
			axs.plot(np.real(self.omegas[i,:]) , self.abs()[i,:], label = str(self.omegas[i,0].imag))
			

		
		axs.set_yscale('log')
		plt.legend(title = r"$\eta$")
		axs.set_xlabel(r"$\omega_{0}$")
		axs.set_ylabel(r"$\det\left[\hat{\mathcal{M}}(\omega) - I \right]$")
		plt.show()

	def minValue(self):
		coords = np.where(self.abs() == np.amin(self.abs()))
		print(coords, np.amin(self.abs()))
		return self.omegas[coords[0][0], coords[1][0]]




# det = ModeFinder("Modes_Data/Unstable_Mode_Search_nu_4.csv")
# det = ModeFinder("../test.csv")
# print(np.real(det.det)[0,:])
# det.omegaPlot()

det = ModeFinder("Modes_Data/Unstable_Mode_Search.csv")
det.omegaPlot()