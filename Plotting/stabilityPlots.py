from generalPlottingFunctions import *

def stabilityPlotsKalnajs(filename):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(ncols = 2, nrows = 2)

	axs[0,0].set_title(r"Varying $m_{1}$")
	data = readingInRealCSV(filename[0])
	timeL = np.linspace(0, 10, np.shape(data)[1]-1)
	for row in range(np.shape(data)[0]):
		axs[0,0].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
		axs[0,0].plot(timeL, data[row, 1:], label = str(data[row, 0]))
		axs[0,0].legend()

	
	axs[0,1].set_title("Varying Time Step")
	data = readingInRealCSV(filename[1])
	for row in range(np.shape(data)[0]):
		time = np.linspace(0, 10, len(data[row])-1)
		axs[0,1].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
		axs[0,1].plot(time, data[row][1:], label = str(data[row][0]))	
		axs[0,1].legend()


	axs[1,0].set_title(r"Varying $\Delta X_{grid}$")
	data = readingInRealCSV(filename[2])
	for row in range(np.shape(data)[0]):
		axs[1,0].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
		axs[1,0].plot(timeL, data[row, 1:], label = str(round(20/data[row, 0],2) ))
		axs[1,0].legend()

	axs[1,1].set_title("Varying Number of Basis Functions")
	data = readingInRealCSV(filename[3])
	for row in range(np.shape(data)[0]):
		axs[1,1].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
		axs[1,1].plot(timeL, data[row, 1:], label = str(data[row, 0]))
		axs[1,1].legend()

	axs[0,0].set_ylabel("Energy")
	axs[1,0].set_ylabel("Energy")

	axs[1,0].set_xlabel(r"$T_{dyn}$")
	axs[1,1].set_xlabel(r"$T_{dyn}$")


	plt.show()


stabilityPlotsKalnajs(["KalnajsTesting/kalnajsVaryingm1.csv", "KalnajsTesting/kalnajsVaryingTimestep.csv", "KalnajsTesting/kalnajsVaryingGridSize.csv", "KalnajsTesting/kalnajsVaryingNbasis.csv"])