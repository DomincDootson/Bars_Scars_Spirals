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
		
		return np.real(np.absolute(self.det))

	def contourPlot(self, axs = None):
		if axs == None:
			fig, axs = plt.subplots()

		XX, YY, logDet = np.real(self.omegas), np.imag(self.omegas), np.log(np.absolute(self.det))
		axs.contourf(XX, YY, logDet, levels = 100)
		#cbar = axs.colorbar()
		#cbar.axs.set_title(r"$\log\left[\det\left[\hat{\mathcal{M}}(\omega) - I \right]\right] $")

		minVal = self.mostUnstableMode()
		print(minVal)

		if minVal != None:
			axs.scatter(minVal.real, minVal.imag, color = 'firebrick')

		if (axs != None):
			return axs
		else:
			plt.show()

	def contourPlotShow(self):
		
		fig, axs = plt.subplots()

		XX, YY, logDet = np.real(self.omegas), np.imag(self.omegas), np.log(np.absolute(self.det))
		contour = axs.contourf(XX, YY, logDet, levels = 100)
		cbar = fig.colorbar(contour)
		cbar.ax.set_title(r"$\log\left[\det\left[\hat{\mathcal{M}}(\omega) - I \right]\right] $")

		minVal =self.modes()
		

		for each in minVal:
			print(each)
			axs.scatter(each.real, each.imag, color = 'firebrick')

		
		plt.show()

	def omegaPlot(self, title = None):
		#plt.rc('text', usetex=True)
		#plt.rc('font', family='serif')
		fig, axs = plt.subplots()
		for i in range(self.nE):
			axs.plot(np.real(self.omegas[i,:]) , self.abs()[i,:], label = str(self.omegas[i,0].imag))
			
		axs.set_yscale('log')
		plt.legend(title = r"$\eta$")
		axs.set_xlabel(r"$\omega_{0}$", fontsize = 13)
		axs.set_ylabel(r"$\det\left[\hat{\mathcal{M}}(\omega) - I \right]$", fontsize = 13)
		if title != None:
			plt.title(title, fontsize = 15)
		plt.show()

	def cutThrough(self, eta):
		plt.plot(self.omegas[:, eta], self.abs()[:, eta])
		plt.show()

	def minValue(self):
		coords = np.where(self.abs() == np.amin(self.abs()))
		return self.omegas[coords[0][0], coords[1][0]]

	def mostUnstableMode(self):
		det = self.abs()

		for i in range(self.nE-2, 0, -1):
			for j in range(1, self.nO-1):
				cc, cl, cr, uc, lc = det[i,j], det[i,j-1], det[i,j+1], det[i-1,j], det[i+1,j]
				if ((cc < cl) & (cc < cr) & (cc < uc) & (cc < lc)):

					return self.omegas[i,j]


	def modes(self):
		det = self.abs()
		modes = []
		for i in range(self.nE-2, 0, -1):
			for j in range(1, self.nO-1):
				cc, cl, cr, uc, lc = det[i,j], det[i,j-1], det[i,j+1], det[i-1,j], det[i+1,j]
				if ((cc < cl) & (cc < cr) & (cc < uc) & (cc < lc)):
					modes.append(self.omegas[i,j])

		return modes
				






#det = ModeFinder("Modes_Data/Unstable_Mode_Search.csv")
#det = ModeFinder("../test.csv") Unstable_Mode_Search_kernel_nu_4
#det = ModeFinder("Modes_Data/Unstable_Mode_Search_kernel_nu_4.csv")
#det = ModeFinder("../scarredModeSearch.csv")
#det = ModeFinder("../SM_5_Upper.csv")
#det = ModeFinder("../scarredMode_Depth_0_Lower_Eta.csv")
#det = ModeFinder("../scarredMode_Depth_10.csv") #
#det = ModeFinder("Modes_Data/Kernel_Search/Single_Scar/AM_Scarred_Kernel_R_20_W_25_D_-95_G.csv") #(0.813924+0.102703j)

#det.contourPlotShow()
# det.cutThrough(20)


