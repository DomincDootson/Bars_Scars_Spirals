from Density_Classes.OneDdensity import *
from Density_Classes.TwoDdensity import *

from generalPlottingFunctions import *
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker

def spiralInitialConditions(fileLeading = "Spiral_data/spiralLeading.csv", fileTrailing = "Spiral_data/spiralTrailing.csv"):
	spiralLeading, spiralTrailing = TwoDdensity(fileLeading), TwoDdensity(fileTrailing)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 4)

	R, phi = np.linspace(0,10, 100), np.linspace(0, 2*3.14,100)
	x, y = np.linspace(-10,10, spiralLeading.nRows), np.linspace(-10,10, spiralLeading.nRows)

	xC, yC= [5 *cos(ang) for ang in phi], [5 *sin(ang) for ang in phi]

	maxValue = spiralLeading.maxValue
		
	levels = np.linspace(0.01*maxValue, 0.99*maxValue, 4)
	axs[0].contour(phi, R, spiralTrailing.densityPolar2D(-1, R, phi, R[-1]), levels = levels, colors = 'k', linewidths = 0.5)
	axs[1].contour(x, y, spiralTrailing.densityAtTime(-1), levels = levels, colors = 'k', linewidths = 0.5)
	axs[0].contourf(phi, R, spiralTrailing.densityPolar2D(-1, R, phi, R[-1]), levels = 50)
	axs[1].contourf(x, y, spiralTrailing.densityAtTime(-1), levels = 50)

	axs[2].contour(phi, R, spiralLeading.densityPolar2D(-1, R, phi, R[-1]), levels = levels, colors = 'k', linewidths = 0.5)
	axs[3].contour(x, y, spiralLeading.densityAtTime(-1), levels = levels, colors = 'k', linewidths = 0.5)
	axs[2].contourf(phi, R, spiralLeading.densityPolar2D(-1, R, phi, R[-1]), levels = 50)
	axs[3].contourf(x, y, spiralLeading.densityAtTime(-1), levels = 50)


	for i in range(2):
		axs[i*2].set_xticks([0, 3.14, 2*3.14])
		axs[2*i].set_xticklabels([0, r"$\pi$",r"$2\pi$"])

		axs[i*2].set_xlabel(r"$\phi$", fontsize = 15)
		axs[i*2].set_ylabel(r"$R$", fontsize = 15)
		axs[i*2+1].set_xlabel(r"$x$", fontsize = 15)
		axs[i*2+1].set_ylabel(r"$y$", fontsize = 15)

	axs[0].set_title("Trailing Spiral")
	axs[1].set_title("Trailing Spiral")
	axs[2].set_title("Leading Spiral")
	axs[3].set_title("Leading Spiral")

	axs[0].axhline(5, color = 'firebrick', linestyle = '--')
	axs[1].plot(xC, yC, color = 'firebrick', linestyle = '--')
	axs[2].axhline(5, color = 'firebrick', linestyle = '--')
	axs[3].plot(xC, yC, color = 'firebrick', linestyle = '--')


	plt.show()	

	


def spiralAnimation(filename = "Spiral_Data/spiral.csv", write2file = None):
	spiral = twoDdensity(filename)
	spiral.densityAnimation(10, write2file)		


def maxRadius(filename = "spiral.csv"):
	density = twoDdensity("spiral.csv")
	fig, ax = plt.subplots()
	#ax.plot(density.maxDensityRadius(10))
	ax.plot(density.maxDensityEvolution())
	# ax2 = ax.twinx()
	# ax2.plot(density.maxDensityEvolution())
	
	plt.show()

def spiralEvolutionWithRhoSquared(filename = "VaryingSpiralTightness.csv", filenameDensity = "spiral.csv"):
	data = readingInRealCSV(filename)
	spiral = twoDdensity(filenameDensity)

	timeIndex = [0,28, 74, 148,799]

	fig = plt.figure(constrained_layout=True)
	gs = GridSpec(2, len(timeIndex), figure=fig)

	axsTop = fig.add_subplot(gs[0, :])
	axsBottom = [fig.add_subplot(gs[-1, i]) for i in range(len(timeIndex))]

	time = np.linspace(0, 50, spiral.nSteps)
	for i in range(np.shape(data)[0]):
		#plt.plot(data[i, 1:], label = str(i))
		axsTop.plot(time, data[i, 1:], label = str(data[i, 0]))		
	axsTop.axhline(1, linestyle ='--')
	axsTop.legend()

	axsTop.set_xlabel(r"$t_{0}$")
	axsTop.set_ylabel(r"$|\rho/\rho(t=0)|^{2}$")
	

	for i in range(len(timeIndex)):
		im = axsBottom[i].imshow(spiral.densityAtTime(timeIndex[i]), extent = [-10,10,-10,10])
		
		fig.colorbar(im, orientation='horizontal', ax = axsBottom[i])
		axsBottom[i].set_title(str(timeIndex[i]*0.25) + r"$T_{dyn}$")

	plt.show()


def varyingK(filenames = "Spiral_Data/VaryingK.csv"):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 2, sharey = True)

	data = readingInRealCSV(filenames)
	time = np.linspace(0,50, np.shape(data)[1]-1)
	colors = ['navy', 'forestgreen', 'firebrick']

	for i, row in enumerate(data):
		column = int(0.5*(row[0]/abs(row[0]) + 1))
		axs[column].plot(time, row[1:], label = r"$|k|=$"+ str(abs(row[0])), color = colors[i%3])
		axs[column].axhline(row[-1], linestyle = '--', color = colors[i%3])
		print(row[0], row[-1])
	

	

	axs[0].set_title("Initally Leading")
	axs[1].set_title("Initally Trailing")

	axs[0].set_ylabel(r"$M_{\mathrm{RMS}}$", fontsize = 15)

	labelPositions, labelsValues = range(0,51, 10), ["0"] + [str(i) + r"$\times 2\pi$" for i in range(10, 51, 10)]
	axs[1].legend()

	for ax in axs:
		ax.set_xticks(labelPositions)
		ax.set_xticklabels(labelsValues)
		
		
		ax.set_xlabel("Time")
	





	plt.show()

def varyingSigma(filename = "Spiral_Data/VaryingSpiralTemperature.csv"):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	data = readingInRealCSV(filename)
	time = np.linspace(0, 25, np.shape(data)[1]-1)

	colors = ['navy', 'cornflowerblue', 'lightcoral', 'firebrick']

	fig, axs = plt.subplots(ncols = 1)
	for i in range(np.shape(data)[0]):
		#plt.plot(data[i, 1:], label = str(i))
		axs.plot(time, np.sqrt(data[i, 1:]), label = str(data[i, 0]) + r"$v_{c}$", color = colors[i])		
	


	axs.axhline(1, linestyle ='--')
	#axs.set_yscale("log")
	plt.xlabel(r"Time", fontsize = 15)
	plt.ylabel(r"$M_{RMS}$", fontsize = 15)

	labelPositions, labelsValues = range(0,26, 5), ["0"] + [str(i) + r"$\times 2\pi$" for i in range(5, 26, 5)]

	axs.set_xticks(labelPositions)
	axs.set_xticklabels(labelsValues)

	axs.legend()
	plt.show()


def polarDensity(rMax = 10, lstK = [.5, .75, 1, -.5, -.75, -1]):
		fig, axs = plt.subplots(ncols = 3, nrows = 2)

		filename = lambda k : "Spiral_Data/Spiral_" + str(round(k*100)) + ".csv"



		spiral = []
		for i in range(np.shape(axs)[0]):
			sp = []
			for j in range(np.shape(axs)[1]):
				sp.append(twoDdensity(filename(lstK[j + i * 3])))
			spiral.append(sp)


		ims = []		
		R, phi = np.linspace(0,10, 100), np.linspace(0, 2*3.14,100)
		

		for time in range(spiral[0][0].nSteps):
			im = []

			for i in range(np.shape(axs)[0]):
				for j in range(np.shape(axs)[1]):
					im.append(axs[i,j].imshow(spiral[i][j].densityPolar2D(time, R, phi, rMax), extent = (0,2 * 3.14,0,rMax,)))
					if (i == 0 & j ==1):
						title = fig.text(.4,.9,(r"Time: " +str(time * 5 * 0.25) + r"$T_{0}$"))
						im.append(title)

			
			ims.append(im)


		
		for i in range(np.shape(axs)[0]):
				for j in range(np.shape(axs)[1]):
					axs[i,j].set_title(r"k = " + str(-lstK[j +i *3]))
		

		fig.suptitle(r"Time step = 0.1 $t_{0}$")			
		ani = animation.ArtistAnimation(fig, ims, interval=30)

		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))
		ani.save("Spiral_Plots/Spiral_Videos/DiffKShortTimeSmallPeak.mp4", writer = writer)
		
		plt.show()
		

def leadingAndTrailing(filenames = ["Spiral_Data/Spiral_100.csv", "Spiral_Data/Spiral_-100.csv"]):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs =plt.subplots(nrows = 5, ncols = 2 * len(filenames), sharex = 'col')

	spirals = [TwoDdensity(file) for file in filenames]
	R, phi = np.linspace(0,10, 100), np.linspace(0, 2*3.14,100)
	x, y = np.linspace(-10,10, spirals[0].nRows), np.linspace(-10,10, spirals[0].nRows)

	nskip = 2
	for i in range(np.shape(axs)[0]):
		if i == 0:
			axs[i,0].set_ylabel("Time = 0" +'\n' + r"$R$")
		else:
			axs[i,0].set_ylabel("Time = " +str(i*nskip) + r"$\times 2\pi$" +'\n' + r"$R$")
		for j in range(len(filenames)):
			maxValue = spirals[j].maxDensityTime(nskip*i)
			absMaxValue = max([spiral.maxDensity() for spiral in spirals])
			levels = np.linspace(0.01*maxValue, 0.99*maxValue, 4)
			axs[i, 2*j].contour(phi, R, spirals[j].densityPolar2D(i*nskip, R, phi, R[-1]), levels = levels, colors = 'k', linewidths = 0.5)
			axs[i,2*j+1].contour(x, y, spirals[j].densityAtTime(i*nskip), levels = levels, colors = 'k', linewidths = 0.5)

			axs[i, 2*j].contourf(phi, R, spirals[j].densityPolar2D(i*nskip, R, phi, R[-1]), levels = 50, vmin = -absMaxValue, vmax = absMaxValue)
			axs[i,2*j+1].contourf(x, y, spirals[j].densityAtTime(i*nskip), levels = 50, vmin = -absMaxValue, vmax = absMaxValue)


	for i in range(2):
		axs[-1, i*2].set_xticks([0, 3.14, 2*3.14])
		axs[-1, 2*i].set_xticklabels([0, r"$\pi$",r"$2\pi$"])

		axs[-1, i*2].set_xlabel(r"$\phi$")
		axs[-1, i*2+1].set_xlabel(r"$x$")

	axs[0,0].set_title("Trailing Spiral" + '\n' + r"R vs $\phi$")
	axs[0,1].set_title("Trailing Spiral" + '\n' + r"Cartesian")
	axs[0,2].set_title("Leading Spiral" + '\n' + r"R vs $\phi$")
	axs[0,3].set_title("Leading Spiral" + '\n' + r"Cartesian")


	plt.show()


def final_coeff(filename = "Spiral_Data/FinalCoeff.csv"):
	data = readingInComplexCSV(filename)

	fig, axs = plt.subplots(ncols = 2)

	for row in data:
		mag = np.sum(np.abs(row))
		axs[0].plot(np.real(row)/mag)
		axs[1].plot(np.imag(row)/mag)

	plt.show()

def fitAmplificationFactor(densityArray, densityArrayPrime):
	return np.sum(densityArray*densityArrayPrime)/np.sum(densityArrayPrime*densityArrayPrime)

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if x ==0:
    	return "0"
    return r'${} \times 10^{{{}}}$'.format(a, b)

def finalResiduals(kValues = [-.5, .5, -.75, .75, -1, 1], fileStem = "Spiral_Data/Spiral_"):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	
	filenames = [fileStem + str(int(100*k)) +".csv" for k in kValues]
	spirals = [TwoDdensity(file) for file in filenames]
	
	fig, axs = plt.subplots(nrows = int(len(kValues)/2), ncols = 4, sharex = 'col')

	R, phi = np.linspace(0,10, 100), np.linspace(0, 2*3.14,100)
	x, y = np.linspace(-10,10, spirals[0].nRows), np.linspace(-10,10, spirals[0].nRows)

	for i, ax in enumerate(axs):
		ax[0].set_ylabel(f"$|k| = {abs(kValues[i*2])}$")
		for j, a in enumerate(ax):
			maxValue = spirals[j].maxDensityTime(-1)
			absMaxValue = max([spiral.maxDensity() for spiral in spirals])
			levels = np.linspace(0.01*maxValue, 0.99*maxValue, 4)
			if j < 2: 
				axs[i, j].contourf(phi, R, spirals[2*i+j].densityPolar2D(-1, R, phi, R[-1]), levels = 50)
				axs[i, j].contour(phi, R, spirals[2*i+j].densityPolar2D(-1, R, phi, R[-1]), levels = levels, colors = 'k', linewidths = 0.5)
			else:
				axs[i,j].contourf(x, y, spirals[2*i+j-2].densityAtTime(-1), levels = 50)
				axs[i, j].contour(x, y, spirals[2*i+j-2].densityAtTime(-1), levels = levels, colors = 'k', linewidths = 0.5)

	axs[0, 0].set_title("Polar Leading")
	axs[0, 1].set_title("Polar Trailing")
	axs[0, 2].set_title("Cartesian Leading")
	axs[0, 3].set_title("Cartesian Trailing")

	for i in range(2):
		axs[-1, i].set_xticks([0, 3.14, 2*3.14])
		axs[-1, i].set_xticklabels([0, r"$\pi$",r"$2\pi$"])

		axs[-1, i].set_xlabel(r"$\phi$")
		axs[-1, i+2].set_xlabel(r"$x$")


	plt.show()


#maxRadius()
#final_coeff()
finalResiduals()
#varyingK()
#spiralInitialConditions()

#finalResiduals()

#density = TwoDdensity("Spiral_Data/Spiral_50.csv")
#density.densityAnimationPolar(10, "Spiral_Plots/Spiral_Videos/LongTrailingWave.mp4")

#finalResiduals()