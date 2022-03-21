from generalPlottingFunctions import *
from matplotlib.gridspec import GridSpec

def spiralInitialConditions(filename = "spiral.csv"):
	spiral = twoDdensity(filename)

	x, y = np.linspace(-10,10, spiral.nRows), np.linspace(-10,10, spiral.nCols)
	XX, YY = np.meshgrid(x,y)

	xCir, yCir = [5*cos(theta) for theta in np.linspace(0, 2*3.14)], [5*sin(theta) for theta in np.linspace(0, 2*3.14)]

	plt.contour(XX, YY, spiral.densityAtTime(0),colors = 'k')
	plt.contourf(XX, YY, spiral.densityAtTime(0), levels = 100)

	plt.plot(xCir, yCir, color = 'firebrick', linestyle = '--')

	plt.colorbar()
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


def varyingK(filenames = "VaryingSpiralTightness.csv"):
	fig, axs = plt.subplots(ncols = len(filenames), sharey = True)

	for i in range(len(filenames)):
		filename = filenames[i]
		data = readingInRealCSV(filename)
		time = np.linspace(0, 50, np.shape(data)[1]-1)
		for j in range(np.shape(data)[0]):
			#plt.plot(data[i, 1:], label = str(i))
			axs[i].plot(time, data[j, 1:], label = str(data[j, 0]))

			axs[i].axhline(1, linestyle ='--')
			axs[i].set_xlabel(r"$t_{0}$")
			axs[i].set_ylabel(r"$|\rho/\rho(t=0)|^{2}$")

		plt.legend()

	axs[0].set_title("Negative K")
	axs[1].set_title("Positive K")
	plt.show()

def varyingSigma(filename = "VaryingSpiralTemperature.csv"):
	data = readingInRealCSV(filename)
	time = np.linspace(0, 50, np.shape(data)[1]-1)

	fig, axs = plt.subplots(ncols = 1)
	for i in range(np.shape(data)[0]):
		#plt.plot(data[i, 1:], label = str(i))
		axs.plot(time, data[i, 1:], label = str(data[i, 0]) + r"$v_{c}$")		
	axs.axhline(1, linestyle ='--')
	axs.set_yscale("log")
	plt.xlabel(r"$t_{0}$")
	plt.ylabel(r"$|\rho/\rho(t=0)|^{2}$")

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
		

def leadingAndTrailing(filenames = ["Spiral_Data/Spiral_-100.csv", "Spiral_Data/Spiral_100.csv"]):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs =plt.subplots(nrows = 5, ncols = 2 * len(filenames), sharex = 'col')

	spirals = [twoDdensity(file) for file in filenames]
	R, phi = np.linspace(0,10, 100), np.linspace(0, 2*3.14,100)
	x, y = np.linspace(-10,10, spirals[0].nRows), np.linspace(-10,10, spirals[0].nRows)

	nskip = 40
	for i in range(np.shape(axs)[0]):
		axs[i,0].set_ylabel(str(nskip*i*0.25) + r"$t_{0}$" +'\n' + r"$R$ $[R_{0}]$")
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
		axs[-1, i*2+1].set_xlabel(r"$x$ $[R_{0}]$")

	axs[0,0].set_title("Trailing Spiral" + '\n' + r"R vs $\phi$")
	axs[0,1].set_title("Trailing Spiral" + '\n' + r"Cartesian")
	axs[0,2].set_title("Leading Spiral" + '\n' + r"R vs $\phi$")
	axs[0,3].set_title("Leading Spiral" + '\n' + r"Cartesian")


	plt.show()



spiralAnimation("Spiral_Data/TightlyWoundNoIC.csv")
#leadingAndTrailing()
#varyingK(["Spiral_Data/VaryingSpiralTightnessTightlyWoundNegativeK.csv", "Spiral_Data/VaryingSpiralTightnessTightlyWoundNegativeK24BF.csv"])
#spiralEvolutionWithRhoSquared()
#maxRadius()

#varyingSigma()

#density = twoDdensity("spiral.csv")
#density.densityAnimation() #It seems 
#density.densityPlots([0,140,165,-1])