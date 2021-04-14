from generalPlottingFunctions import *

def barDensity():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	density2D = readingInRealCSV("density2d.csv")
	nCols = np.shape(density2D)[0]
	spacing = 20/(nCols-1)
	centre = (nCols-1)*0.5

	x = np.linspace(-5,5, nCols)
	y = np.linspace(-5,5, nCols)#np.arange(-10,10, spacing)
	XX, YY = np.meshgrid(x, y)
	
	fig, axs = plt.subplots(1,1)

	contours = [-.45, -.35, -.15, .15, .35, .45]
	axs.contour(XX, YY, density2D, contours, colors = 'k')
	density = axs.contourf(XX, YY, density2D, 100)
	cbar = plt.colorbar(density, ticks = [0.1*each for each in range(-5, 6)])
	cbar.set_label(r"$\rho_{b}/\left(\epsilon m_{b}\right)$", fontsize = 15)

	xCir, yCir = [2*cos(theta) for theta in np.linspace(0, 2*3.142)], [2*sin(theta) for theta in np.linspace(0, 2*3.142)]
	axs.plot(xCir, yCir, color = 'firebrick')
	axs.set_title("Toy Bar Model", fontsize = 15)
	
	plt.show()


def TorqueEvolution():
	data = readingInRealCSV("barEvolution.csv")
	plt.plot(data[:,-1])
	plt.show()

def deltaTorquePlot():
	data = readingInRealCSV("TorqueDelta.csv")
	plt.plot(data)
	plt.show()

def smearedTorque():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	data = readingInRealCSV("SmearedTorque.csv")/3.14
	t = np.linspace(0, 20, len(data))
	plt.plot(t, data, label = "Smeared Torque")

	data = readingInRealCSV("TorqueDelta.csv")
	t = np.linspace(0, 20, len(data))
	plt.plot(t, data, label = "Delta Density", alpha = .8)

	data = readingInRealCSV("barEvolution.csv")
	t = np.linspace(0, 20, len(data[:,0]))
	plt.plot(t, data[:,-1], label = "Full Torque")

	data = readingInRealCSV("FullTorque.csv")
	t = np.linspace(0, 20, len(data))
	plt.plot(t, data, label = "Integrated Torque")

	plt.legend()
	plt.xlabel(r"Time")
	plt.ylabel(r"Torque / $\epsilon m_{b}$")
	plt.show()

#TorqueEvolution()
#deltaTorquePlot()
smearedTorque()

