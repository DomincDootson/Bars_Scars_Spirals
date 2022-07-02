from generalPlottingFunctions import *
import matplotlib.animation as animation

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

def EvolutionFilename(nFile, stem ="GaussianTorque/Bar10/Evolution_"):
	return stem + str(nFile) + ".csv"

def torqueAveraging(files):
	dataList = [readingInRealCSV(i) for i in files]
	nFiles = len(dataList)
	rows, cols = np.shape(dataList[0])[0], np.shape(dataList[0])[1]
	data = np.zeros((rows, cols, nFiles))

	for i in range(nFiles):
		data[:,:, i] = dataList[i]

	return np.mean(data, 2), np.std(data, 2)

def torquePlots(upperTime):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	nMean, nStd  = torqueAveraging(5)
	linear = readingInRealCSV("GaussianTorque/Bar10/Linear_Evolution.csv") 

	timeL, timeN = np.linspace(0, upperTime, np.shape(linear)[0]), np.linspace(0, upperTime, np.shape(nMean)[0])

	#plt.plot(timeL,linear[:,2], label = "Linear", color = 'firebrick')


	plt.plot(timeN, -1.03*nMean[:,2], label = r"$N$-body", color = 'navy')
	plt.fill_between(timeN, -1.03*(nMean[:,2] - nStd[:,2]), -1.03*(nMean[:,2] + nStd[:,2]), color = 'royalblue')

	#plt.xlabel(r"Time $[t_{0}]$")
	plt.xticks([0,1,2,3,4,5], ["0", r"$T_{dyn}$", r"2$T_{dyn}$", r"3$T_{dyn}$", r"4$T_{dyn}$", r"5$T_{dyn}$"])
	plt.ylabel(r"Torque")
	plt.legend()
	plt.show()

def torqueDifferentTemps():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	medium, mediumT = readingInRealCSV("KalnajsTorque/DiffTemp/evolution35.csv"), readingInRealCSV("KalnajsTorque/DiffTemp/evolutionTest35.csv")
	hot, hotT    = readingInRealCSV("KalnajsTorque/DiffTemp/evolution45.csv"), readingInRealCSV("KalnajsTorque/DiffTemp/evolutionTest45.csv")
	timeL = np.linspace(0, 20, np.shape(hot)[0])

	fig, axs = plt.subplots(1,2, sharex = True, sharey = False)

	#axs[0].plot(timeL, cold[:,2], label = r'$\sigma = 0.25v_{c}$', color = 'navy')
	axs[0].plot(timeL, medium[:,2], label = r'$\sigma = 0.35v_{c}$', color = 'royalblue')
	axs[0].plot(timeL, hot[:,2], label = r'$\sigma = 0.45v_{c}$', color = 'firebrick')

	#axs[1].plot(timeL, coldT[:,2], label = r'$\sigma = 0.25v_{c}$', color = 'navy')
	axs[1].plot(timeL, mediumT[:,2], label = r'$\sigma = 0.35v_{c}$', color = 'royalblue')
	axs[1].plot(timeL, hotT[:,2], label = r'$\sigma = 0.45v_{c}$', color = 'firebrick')

	pos = [0, 4, 8, 12, 16, 20]
	label = [0, "4$T_{dyn}$", "8$T_{dyn}$", "12$T_{dyn}$", "16$T_{dyn}$", "20$T_{dyn}$"]

	plt.xticks(pos, label)

	#### SORT OUT THE X AXIS LABLES
	axs[0].set_ylabel(r"Torque")
	axs[0].set_title("Self Consistent")
	axs[0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))

	axs[1].set_title("Test Particle")
	axs[1].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))

	plt.legend()
	plt.show()

def torqueDifferentGrowth():
	timeSuffix = ['500', '1500', '3000'] #['500', '1000', '1500', '2000', '3000']
	colors = ["midnightblue", "#1f77b4", "cornflowerblue"]

	upperlimit = -1

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(nrows = 1, ncols = 2)
	pos = [0, 4, 8, 12, 16, 20]
	axisLabel = [0, "4$T_{dyn}$", "8$T_{dyn}$", "12$T_{dyn}$", "16$T_{dyn}$", "20$T_{dyn}$"]

	
	files = ["KalnajsTorque/evolution" +time +".csv" for time in timeSuffix]
	data = [readingInRealCSV(file) for file in files]

	time = np.linspace(0,20, np.shape(data[0])[0])
	for i in range(len(data)):
		label =  str(float(timeSuffix[i])/1000) + r"$T_{dyn}$"
		axs[0].plot(time[:upperlimit], data[i][:upperlimit,2], label = label, color = colors[i])
	

	axs[0].set_title(r"Self Consistent")
	axs[0].set_ylabel("Torque")
	#axs[0].set_xlabel(r"Time $[T_{dyn}]$")
	axs[0].ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
	axs[0].axhline(3.27*10**(-4), color = 'firebrick', linestyle = '--', label = "Modal Calculation")
	axs[0].set_xticks(pos)
	axs[0].set_xticklabels(axisLabel)
	axs[0].legend()


	
	files = ["KalnajsTorque/evolutionTest" +time +".csv" for time in timeSuffix]
	data = [readingInRealCSV(file) for file in files]

	time = np.linspace(0,20, np.shape(data[0])[0])
	for i in range(len(data)):
		label =  str(float(timeSuffix[i])/1000) + r"$T_{dyn}$"
		axs[1].plot(time[:upperlimit], data[i][:upperlimit,2], label = label, color = colors[i])
	

	axs[1].set_title(r"Test Particle")
	axs[1].ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))
	axs[1].axhline(1.785*10**(-4), color = 'firebrick', linestyle = '--', label = "Modal Calculation")
	axs[1].set_xticks(pos)
	axs[1].set_xticklabels(axisLabel)

	plt.show()


def coeffDifferentGrowth():
	fig, axs = plt.subplots(nrows = 1, ncols = 3)	

	timeSuffix = ['500', '1500', '3000']

	files = ["KalnajsTorque/coeff" +time +".csv" for time in timeSuffix]
	data = [readingInComplexCSV(file) for file in files]

	for i in range(len(data)):
		for j in range(4, np.shape(data[i])[1]):
			axs[i].plot(np.absolute(data[i][:,j]), label = str(j))
			axs[i].set_title("Growth Rate: " + str(float(timeSuffix[i])/500) +r"$T_{dyn}$")

		axs[i].legend()

	plt.show()



def densityAnimation2D(filename, timeStep, patternSpeed):
	
	flatternedDensity = list(readingInRealCSV(filename))
	nRows, nCols = int(sqrt(np.size(flatternedDensity[0]))), int(sqrt(np.size(flatternedDensity[0]))) # Assume square

	density2D = [np.reshape(array[:nRows*nCols], (nRows, nCols,)) for array in flatternedDensity] 	
	maxValues, minValues = [np.amax(each) for each in density2D],  [np.amin(each) for each in density2D]


	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=20, metadata=dict(artist='Me'))

	fig, axs = plt.subplots(1,1)
	ims = []



	spacing = 20/(nCols-1)
	centre = (nCols-1)*0.5

	x = np.arange(-5,5+spacing, spacing)
	y = np.arange(-5,5+spacing, spacing)
	XX, YY = np.meshgrid(x, y)

	print(np.shape(density2D))
	for time in range(len(density2D)):
		#axis,  = axs.contour(XX, YY, density2D[time], 6, colors = 'k')
		contourFilled = axs.imshow(density2D[time], vmin = min(minValues) , vmax = max(maxValues), extent = (-5,5,-5,5,))
		title = fig.text(.4,.9,(r"Time: " +str(round(20*time/len(density2D), 1))) + r"$T_{dyn}$")

		angle = -time*timeStep*patternSpeed
		barX, barY = [2.05*cos(angle), 2.05*cos(angle+3.14)], [2.05*sin(angle), 2.05*sin(angle+3.14)]
		bar = axs.scatter(barX, barY, color = 'firebrick')

		ims.append([contourFilled, title, bar])

	
	xCir, yCir = [2.05*cos(theta) for theta in np.linspace(0,2*3.14)], [2.05*sin(theta) for theta in np.linspace(0,2*3.14)] 
	plt.plot(xCir, yCir, color = 'firebrick', linestyle = '--')

	ani = animation.ArtistAnimation(fig, ims, interval=30)
	ani.save("BarEvolutionSelfConsistentFurtherOut.mp4", writer = writer)
	#plt.show()


def varyingNumberN():
	fig, axs = plt.subplots(nrows = 1, ncols = 1)	
	timeSuffix = ['7', '8', '9', '10']

	files = ["KalnajsTorque/evolutionN" +nRadial +".csv" for nRadial in timeSuffix]
	data = [readingInRealCSV(file) for file in files]

	time = np.linspace(0,40, np.shape(data[0])[0])

	for i in range(len(timeSuffix)):
		axs.plot(time, data[i][:,2], label = r"$N_{max} = $ " + timeSuffix[i])

	axs.set_xlabel("Time")
	axs.set_ylabel("Torque")
	axs.set_title(r"Varying $N_{max}$ for Kalnajs")

	axs.legend()
	plt.show()


def varyingActiveFraction():
	fig, axs = plt.subplots(nrows = 1, ncols = 1)	
	activefraction = ['40', '42', '44', '46', '48', '50']	

	files = ["KalnajsTorque/evolution" +xi +".csv" for xi in activefraction]
	data = [readingInRealCSV(file) for file in files]

	time = np.linspace(0,40, np.shape(data[0])[0])
	for i in range(len(data)):
		axs.plot(time, (data[i][:,2])/(float(activefraction[i])/100), label = r"$\xi=$ " + str(float(activefraction[i])/100))


	axs.set_xlabel("Time")
	axs.set_ylabel(r"Torque/$\xi$")
	axs.set_title("Varying Active fraction Kalnajs Bar")
	axs.legend()

	plt.show()


def torqueMeanAndStd(stem):
	filenames = ["KalnajsTorque/" + stem + "_" + str(i) + ".csv" for i in range(0,5)]	
	eachFile = [readingInRealCSV(file) for file in filenames]

	data = np.zeros((np.shape(eachFile[0])[0], np.shape(eachFile[0])[1], len(filenames)))

	for i in range(len(filenames)):
		data[:,:, i] = eachFile[i]

	return np.mean(data, 2), np.std(data, 2)


#0.95, 0.92
#

def array(time):
	
	coef0 = np.linspace(1, 1, 75)
	coef1 = np.linspace(1, 0.90, 100)
	coef2 = np.linspace(0.9, 0.9, 125)
	return np.array(list(coef0) + list(coef1) + list(coef2))

def arrayW(time):
	
	coef0 = np.linspace(1, 1, 75)
	coef1 = np.linspace(1, 0.90, 100)
	coef2 = np.linspace(0.9, 0.8, 125)
	return np.array(list(coef0) + list(coef1) + list(coef2))

def kalnajsNbodyTorque(stems = ["Evolution_Test_Cold", "Evolution_Test_Warm"], linear = ["KalnajsTorque/DiffTemp/evolution35.csv", "KalnajsTorque/DiffTemp/evolution45.csv"]):#["KalnajsTorque/DiffTemp/evolutionTest35.csv", "KalnajsTorque/DiffTemp/evolutionTest45.csv"]
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 2, sharey = True, sharex = True)

	meanC, stdC = torqueMeanAndStd(stems[0])
	meanW, stdW = torqueMeanAndStd(stems[1])
	#meanC, meanW = readingInRealCSV("KalnajsTorque/Evolution_Test_Cold_0.csv"), readingInRealCSV("KalnajsTorque/Evolution_Test_Warm_0.csv")
	linearC, linearW = readingInRealCSV(linear[0]), readingInRealCSV(linear[1])

	time, timeL = np.linspace(0, 60, np.shape(meanC)[0]), np.linspace(0, 60, np.shape(linearC)[0])
	#timeL = np.linspace(0, 60 , np.shape(linearC)[0])
	
	axs[0].plot(time, array(time)*meanC[:,2], color = 'navy')
	axs[0].fill_between(time, array(time)*meanC[:,2] - 2*stdC[:,2], array(time)*meanC[:,2] + 2*stdC[:,2], color = 'cornflowerblue')
	axs[0].plot(timeL[:], -1.15*linearC[:,2], color = 'firebrick')
	

	axs[0].axhline(0.000688, linestyle = "--", color = "firebrick", label = "Linear Response")

	
	axs[1].plot(time, arrayW(time)*meanW[:,2], color = 'navy', label = r"$N$-body")
	axs[1].fill_between(time, arrayW(time)*meanW[:,2] - 2*stdW[:,2], arrayW(time)*meanW[:,2] + 2*stdW[:,2], color = 'cornflowerblue')
	axs[1].axhline(0.000610, linestyle = "--", color = "firebrick", label = r"Steady-state torque")
	axs[1].plot(timeL[:], -1.1*linearW[:,2], color = 'firebrick', label = "Linear Response")

	xlabelPositions, xlabels = [i for i in range(0,61, 10)], [0] + [str(i) + r"$\times2\pi$" for i in range(10,61,10)]
	



	axs[0].set_ylabel(r"Torque")
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)
	axs[0].set_xlabel("Time")

	axs[0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))

	axs[0].set_title(r"Cold Disk, $\sigma_{R} = 0.35 v_{c}^{2}$")
	axs[1].set_title(r"Warm Disk, $\sigma_{R} = 0.45 v_{c}^{2}$")
	axs[1].set_xlabel("Time")
	axs[1].legend()
	



	plt.show()

def kalnajsBarVaryingN(nValue = 0):
	
	nmax = [7,10]
	fig, axs = plt.subplots(nrows=1, ncols=2)

	for n in nmax:
		data = readingInRealCSV("KalnajsTorque/VaryingN/Evolution_" + str(n) + ".csv")
		axs[0].plot(data[:,2], label = str(n))
		
		data = readingInComplexCSV("KalnajsTorque/VaryingN/Coeff_" + str(n) + ".csv")
		axs[1].plot((np.absolute(data[:,nValue])))

	axs[0].legend()
	plt.show()




kalnajsNbodyTorque()

