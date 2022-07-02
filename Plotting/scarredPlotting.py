from generalPlottingFunctions import *
from Density_Classes.OneDdensity import *

from matplotlib import ticker, cm
from scipy.stats import linregress

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')


def circularOrbitCutThrough():
	scarred, unscarred = readingInRealCSV("Scar_Data/withScar.csv"), readingInRealCSV("Scar_Data/withoutScar.csv")
	fig, axs = plt.subplots(ncols = 1)

	l = np.linspace(0, 3, np.shape(scarred)[0])

	axs.plot(l,scarred, label = "Ungrooved Mestel", color = 'firebrick')
	axs.plot(l,unscarred, label = "Grooved Mestel", color = 'royalblue', linestyle = '--')

	axs.set_xlabel(r"$J_{\phi}$ $[R_{0}v_{c}]$")
	axs.set_ylabel(r"$F_{s}(J_{\phi}, J_{r} = 0)$")
	plt.legend()

	plt.show()

def densityDifferentTemp(temp = [10, 15, 20, 25, 30, 35]):
	fig, axs = plt.subplots()

	for t in temp:
		filename = "Scar_Data/Angular_Momentum_Densities/scarredDensity" + str(t) + ".csv"
		data = readingInRealCSV(filename)

		axs.plot(data[:350,0],data[:350,1], label = r"$\sigma_{R}=$ " + str(t/100))

	plt.legend()
	axs.set_xlabel(r"Radius $[r_{0}]$")
	axs.set_ylabel(r"Density")
	axs.set_title("Angular Momentum Ridges")
	plt.show()


def densityPlot():
	scarred, unscarred   = readingInRealCSV("Scar_Data/scarredDensity.csv"), readingInRealCSV("Scar_Data/unscarredDensity.csv")
	scarredP, unscarredP = readingInRealCSV("Scar_Data/scarredPotential.csv"), readingInRealCSV("Scar_Data/unscarredPotential.csv")
	fig, axs = plt.subplots(ncols = 2)

	axs[0].plot(unscarred[:,0],  unscarred[:,1], label = "Ungrooved Mestel", color = 'firebrick')
	axs[0].plot(scarred[:,0],   scarred[:,1], label = "Grooved Mestel", color = 'royalblue', linestyle = '--')

	axs[0].set_xlabel(r"$R$ $[R_{0}]$")
	axs[0].set_ylabel(r"$\rho(R)$ $[R_{0}^{-2}]$")
	axs[0].set_title(r"Density")
	axs[0].legend()

	axs[1].plot(unscarredP[:,0],  unscarredP[:,1], label = "Ungrooved Mestel", color = 'firebrick')
	axs[1].plot(scarredP[:,0],    scarredP[:,1], label = "Grooved Mestel", color = 'royalblue', linestyle = '--')

	axs[1].set_xlabel(r"$R$ $[R_{0}]$")
	axs[1].set_ylabel(r"$\Phi(R)$ $[v_{c}^{2}]$")
	axs[1].set_title(r"Potential")
	plt.show()

def scarredDFPlot():
	unscarred, scarred = readingInRealCSV("Scar_Data/unscarredDF.csv"), readingInRealCSV("Scar_Data/scarredDF.csv")
	data = np.log(scarred)
	energy = np.linspace(0.5 +log(0.1), 0.5 +log(30), np.shape(data)[0])
	angMom = np.linspace(0, 20, np.shape(data)[1])

	L, E = np.meshgrid(angMom, energy)
	
	plt.contourf(L,E, data, levels = 100)
	#plt.imshow(np.log(data))
	plt.colorbar()
	plt.show()


def fittingStraightline(x, y):
	logY = np.log(np.absolute(y))
	m, c = linregress(x, logY)[0:2]
	'''plt.plot(x, x*m + c)
	plt.show()'''
	print(m, c)
	return x, np.exp(x*m + c)



def scarredModes(filename):
	fig, axs = plt.subplots()
	data = readingInComplexCSV(filename)
	energy = [np.sum(np.absolute(data[i,:])) for i in range(np.shape(data)[0])]
	
	t = np.linspace(0,np.shape(data)[0]*0.25, np.shape(data)[0])

	tFit, fit = fittingStraightline(t[1:], np.absolute(data[1:,5]))

	axs.plot(t, np.absolute(data[:,-1]))
	axs.plot(t, np.absolute(data[:,5]))
	axs.plot(tFit, fit)
	print(np.absolute(data[1,-1]), np.absolute(data[-1,-1]))



	axs.set_yscale("log")
	plt.show()

def angularMomentumScarNudges(filenames = ["Scar_Data/scarredEvolutionModesS.csv", "Scar_Data/scarredEvolutionModesT.csv", "Scar_Data/evolutionModesS.csv", "Scar_Data/evolutionModesT.csv"]): 
	fig, axs = plt.subplots(nrows = 2, ncols = 2, sharex = True)

	for i in range(len(filenames)):
		data = readingInComplexCSV(filenames[i])
		time = np.linspace(0,20, np.shape(data)[0])
		[axs[i//2, i%2].plot(time, np.absolute(data[:, n]), label = str(n)) for n in range(0,10, 2)]
		axs[i//2, i%2].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))

	axs[0,0].set_ylabel(r"Scarred")
	axs[1,0].set_ylabel(r"Unscarred")

	axs[0,0].set_title("Self Consistent")
	axs[0,1].set_title("Test Particle")

	axs[1,0].set_xlabel(r"Time $[t_{0}]$")
	axs[1,1].set_xlabel(r"Time $[t_{0}]$")

	axs[0,0].legend()


	plt.show()

def infallingSatelliteEvolution(filenames = ["Scar_Data/infallingRingSS.csv", "Scar_Data/infallingRingUS.csv", "Scar_Data/infallingRingST.csv", "Scar_Data/infallingRingUT.csv"]):
	infallers = [OneDdensity(filename) for filename in filenames] 
	
	
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=20, metadata=dict(artist='Me'))

	fig, axs = plt.subplots(2,2)
	ims = []

	
	
	for time in range(infallers[0].nStep):
		ss, = axs[0,0].plot(infallers[0].radii, infallers[0].density(time), color = 'royalblue')
		comp, = axs[0,0].plot(infallers[1].radii, infallers[1].density(time), color = 'firebrick')
		us, = axs[0,1].plot(infallers[1].radii, infallers[1].density(time), color = 'royalblue')
		st, = axs[1,0].plot(infallers[2].radii, infallers[2].density(time), color = 'royalblue')
		ut, = axs[1,1].plot(infallers[3].radii, infallers[3].density(time), color = 'royalblue')
		
		ims.append([ss,comp, us, st, ut])


	axs[0,0].set_title("Scarred")
	axs[0,1].set_title("Unscarred")

	axs[0,0].set_ylabel("Self Consistent")
	axs[1,0].set_ylabel("Test Particle")

	ani = animation.ArtistAnimation(fig, ims, interval=30)

	#ani.save(filename, writer = writer)
	plt.show()




#scarredDFPlot()
#scarredModes("scarredEvolution.csv")
#densityPlot()

infallingSatelliteEvolution()
#angularMomentumScarNudges()
#densityDifferentTemp()

