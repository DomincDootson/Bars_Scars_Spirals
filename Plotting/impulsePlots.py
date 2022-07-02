import matplotlib.animation as animation
from Density_Classes.OneDdensity import *

from generalPlottingFunctions import *

def impulseDensityName(littleSigma, radius, angHarmonic, testParticle = False):
	filename = "../Disk_Kicking/littleSigma_" + str(round(littleSigma*100)) + "/Density" + str(round(radius*10)) +"_" + str(angHarmonic) 
	if testParticle:
		filename += "_Test"
	return (filename + ".csv")

def twoDdensityName(littleSigma, radius, angHarmonic, testParticle = False):
	filename = "../Disk_Kicking/littleSigma_" + str(round(littleSigma*100)) + "/Density2D" + str(round(radius*10)) +"_" + str(angHarmonic) 
	if testParticle:
		filename += "_Test"
	return (filename + ".csv")

def densityEvolutionPlot(littleSigma, timestep = 0.25):
	
	radii = [1, 5]
	colors = ["#1F75FE", "#191970"]
	angHarmonic = [0, 1, 2]
	times = [25, 50, 75, 100]


	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	data = []

	for time in times:
		holdingAng = []
		for l in angHarmonic:
			holdingRadius = []
			for r in radii:
				holdingRadius.append(readingInRealCSV(impulseDensityName(littleSigma, r, l)))
			holdingAng.append(holdingRadius)
		data.append(holdingAng)

	fig, axs = plt.subplots(len(times),len(angHarmonic), sharex = True)

	for i in range(len(times)):
		for j in range(len(angHarmonic)):
			for r in range(len(radii)):
				axs[i,j].plot(data[i][j][r][0,:], data[i][j][r][times[i],:], label = str(radii[r])+r"$r_{0}$", color = colors[r])

	

	for l in angHarmonic:
		axs[0,l].set_title(r"$\ell_{p} = $"+str(l))

	for t in range(len(times)):
		axs[t, 0].set_ylabel(str(times[t]*timestep)+r"$t_{0}$")
	

	axs[-1,1].set_xlabel(r"Radius")
	
	axs[0,0].legend()
	fig.suptitle(r"$\sigma_{r} =$ "+str(littleSigma))
	
	plt.show()

	


def densityAtMaxEnergy(energySelf, energyTest, radii):
	littleSigma = energySelf.little_sigma()
	angHarmonic = energySelf.ang_harmonic()

	maxTimeIndexSelf, maxTimeIndexTest = [], []
	for s in range(np.shape(littleSigma)[0]):
		holdingSelf, holdingTest =[], []
		for l in range(np.shape(angHarmonic)[0]):
			holdingSelf.append(energySelf.index_Max_Energy(littleSigma[s], angHarmonic[l], radii[l]))
			holdingTest.append(energyTest.index_Max_Energy(littleSigma[s], angHarmonic[l], radii[l]))

		maxTimeIndexSelf.append(holdingSelf)
		maxTimeIndexTest.append(holdingTest)

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(np.shape(littleSigma)[0], np.shape(angHarmonic)[0], sharex = True)

	data = []
	for s in littleSigma:
		holdingAng = []
		for l in angHarmonic:
			holdingRadius = []
			holdingRadius.append(readingInRealCSV(impulseDensityName(s, radii[int(l)], int(l))))
			holdingRadius.append(readingInRealCSV(impulseDensityName(s, radii[int(l)], int(l), True)))
			
			holdingAng.append(holdingRadius)
		data.append(holdingAng)

	for i in range(np.shape(littleSigma)[0]):
		for j in range(np.shape(angHarmonic)[0]):
			print(maxTimeIndexSelf[i][j])
			axs[i,j].plot(data[i][j][0][0,:], data[i][j][0][maxTimeIndexSelf[i][j],:], color = 'firebrick')
			axs[i,j].plot(data[i][j][1][0,:], data[i][j][1][maxTimeIndexTest[i][j],:], color = 'cornflowerblue')

	for l in range(np.shape(angHarmonic)[0]):
		axs[0,l].set_title(r"$\ell_{p}=$ " + str(int(angHarmonic[l])) + r" $r_{n}=$ "+str(radii[l])+r"$r_{0}$")

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r}=$ " + str(littleSigma[s]) + r"$v_{c}$")

	axs[-1, 1].set_xlabel(r"Radius $[r_{0}]]$")
	fig.suptitle("Density at Peak Energy")

	plt.show()



def densityAnimation(littleSigma, radius, timestep = 0.25):

	angHarmonic = [0,1,2]
	data = []
	furtherData = []

	for l in angHarmonic:
		data.append(readingInRealCSV(impulseDensityName(littleSigma, 1, l)))
		furtherData.append(readingInRealCSV(impulseDensityName(littleSigma, 5, l)))

	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=300, metadata=dict(artist='Me'))

	fig, axs = plt.subplots(1,3)
	ims = []
	axs[0].plot([0],[0], color = 'royalblue', label = str(1)+r"$r_{0}$")
	axs[0].plot([0],[0], color = 'firebrick', label = str(5)+r"$r_{0}$")

	for t in range(1,np.shape(data[0])[0]-1):
		
		line0, = axs[0].plot(data[0][0,:], data[0][t,:], color = "royalblue")
		line3, = axs[0].plot(furtherData[0][0,:], furtherData[0][t,:], color = "firebrick")
		
		line1, = axs[1].plot(data[0][0,:], data[1][t,:], color = "royalblue")
		line4, = axs[1].plot(furtherData[0][0,:], furtherData[1][t,:], color = "firebrick")

		
		line2, = axs[2].plot(data[0][0,:], data[2][t,:], color = "royalblue")
		line5, = axs[2].plot(furtherData[0][0,:], furtherData[2][t,:], color = "firebrick")

		ims.append([line0, line3, line1, line4, line2, line5])

	
	axs[0].legend(loc = 'lower right')

	axs[1].set_xlabel("Radius")
	axs[0].set_ylabel("Density")

	axs[0].set_title(r"$\ell_{p} = 0$")
	axs[1].set_title(r"$\ell_{p} = 1$")
	axs[2].set_title(r"$\ell_{p} = 2$")

	fig.suptitle(r"$\sigma_r =$ " +str(littleSigma))

	ani = animation.ArtistAnimation(fig, ims, interval=30, blit=True)
	plt.show()

def densityAnimation2D(littleSigma, radius, angHarmonic):
	flatternedDensity = list(readingInRealCSV(twoDdensityName(littleSigma, radius, angHarmonic, False)))
	nRows, nCols = int(sqrt(np.size(flatternedDensity[0]))), int(sqrt(np.size(flatternedDensity[0]))) # Assume square


	density2D = [np.reshape(array[:nRows*nCols], (nRows, nCols,)) for array in flatternedDensity] 	
	maxValues, minValues = [np.amax(each) for each in density2D],  [np.amin(each) for each in density2D]
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=20, metadata=dict(artist='Me'))

	fig, axs = plt.subplots(1,1)
	ims = []



	spacing = 20/(nCols-1)
	centre = (nCols-1)*0.5

	x = np.arange(-5,5+spacing, spacing)
	y = np.arange(-5,5+spacing, spacing)
	XX, YY = np.meshgrid(x, y)
	
	#axs.xaxis.set_animated(True)
	xCir, yCir = [radius*cos(theta) for theta in np.linspace(0,2*pi)], [radius*sin(theta) for theta in np.linspace(0,2*pi)]
	#plt.plot(xCir, yCir, color = 'firebrick')
	for time in range(len(density2D)):
		#axis,  = axs.contour(XX, YY, density2D[time], 6, colors = 'k')
		contourFilled = axs.imshow(density2D[time], vmin = min(minValues) , vmax = max(maxValues), extent = (-5,5,-5,5,))
		title = fig.text(.4,.9,(r"Time: " +str(round(8*time/len(density2D), 1))) + r"$T_{dyn}$")
		ims.append([contourFilled, title])


	#fig.subplots_adjust(right=0.8)
	#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	#cbar = plt.colorbar(contourFilled, cax=cbar_ax)

	ani = animation.ArtistAnimation(fig, ims, interval=30)
	#plt.show()
	ani.save("GreensFunctionUnstable.mp4", writer = writer)



densityAnimation2D(.25, 2, 2)


#densityAnimation(0.25, 1)
