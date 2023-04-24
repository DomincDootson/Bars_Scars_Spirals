from generalPlottingFunctions import *
import matplotlib.animation as animation
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from Density_Classes.TwoDdensity import *

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




def kalnajsNbodyTorque(stems = ["Monari_Bar/Evolution_Consistent_Cold", "Monari_Bar/Evolution_Test_Warm"], linear = ["Nbody_Sormani_Data/Varying_Ep/Linear_Consistent_Particle.csv", "Nbody_Sormani_Data/Varying_Ep/Linear_Consistent_Particle_Warm.csv", "Nbody_Sormani_Data/Varying_Ep/Linear_Test_Particle.csv", "Nbody_Sormani_Data/Varying_Ep/Linear_Test_Particle_Warm.csv"]):#["KalnajsTorque/DiffTemp/evolutionTest35.csv", "KalnajsTorque/DiffTemp/evolutionTest45.csv"]
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 2, nrows = 2, sharey = 'row', sharex = True)#,  sharey = True, sharex = True

	# meanC, stdC = torqueMeanAndStd(stems[0])
	# meanW, stdW = torqueMeanAndStd(stems[1])
	
	linear  = [readingInRealCSV(file) for file in linear]
	#timeL = np.linspace(0,100, np.shape(linear[2])[0])

	
	# axs[0,1].plot(timeL[:], linear[1][:,2], color = 'firebrick')
	
	# axs[1,1].plot(timeL[:], linear[3][:,2], color = 'firebrick')
	ep2  = 0.005**2
	axs[0,0].axhline(-0.000970/ep2, linestyle = "--", color = "firebrick")
	axs[0,1].axhline(-0.000910/ep2, linestyle = "--", color = "firebrick")
	axs[1,0].axhline(-0.000865/ep2, linestyle = "--", color = "firebrick")
	axs[1,1].axhline(-0.000732/ep2, linestyle = "--", color = "firebrick")
	# axs[1,0].axhline(-0.00200, linestyle = "--", color = "firebrick", label = r"Steady-state Torque")
	# axs[1,1].axhline(-0.00195, linestyle = "--", color = "firebrick")

	## Cold Consistent ##
	stem = "evolutionC_001_18"
	mean, std = torqueMeanAndStd(stem, upperIndex = 3, file = "Nbody_Sormani_Data/Varying_Ep/")
	time = np.linspace(0,100, np.shape(mean)[0])
	axs[0,0].plot(time, mean[:,2]/(0.001**2), color = 'navy')
	axs[0,0].fill_between(time, (mean[:,2] - std[:,2])/(0.001**2), (mean[:,2] + std[:,2])/(0.001**2), color = 'cornflowerblue')
	timeL = np.linspace(0,125, np.shape(linear[0])[0])
	axs[0,0].plot(timeL[:-2], 1.25*linear[0][2:,2]/ep2, color = 'firebrick')

	# mean, std = torqueMeanAndStd("Sormani_Bar/Evolution_Consistent_Warm")
	# time = np.linspace(0,100, np.shape(mean)[0])
	# axs[0,1].plot(time, mean[:,2], color = 'navy')
	# axs[0,1].fill_between(time, mean[:,2] - std[:,2], mean[:,2] + std[:,2], color = 'cornflowerblue')

	## Cold Test ## 
	stem = "evolution_005_18"
	mean, std = torqueMeanAndStd(stem, upperIndex = 4, file = "Nbody_Sormani_Data/Varying_Ep/")
	time = np.linspace(0,100, np.shape(mean)[0])
	axs[1,0].plot(time, mean[:,2]/ep2, color = 'navy', label = r"$N$-body")
	axs[1,0].fill_between(time, (mean[:,2] - std[:,2])/ep2, (mean[:,2] + std[:,2])/ep2, color = 'cornflowerblue')
	timeL = np.linspace(0,100, np.shape(linear[2])[0])
	axs[1,0].plot(timeL[:-1], linear[2][1:,2]/ep2, color = 'firebrick', label = "Linear Response")
	axs[1,0].legend(fontsize = 12)

	## Warm Consistent ## 
	stem = "evolutionCW_001_18"
	mean, std = torqueMeanAndStd(stem, upperIndex = 3, file = "Nbody_Sormani_Data/Varying_Ep/")
	time = np.linspace(0,100, np.shape(mean)[0])
	axs[0,1].plot(time, 0.95*0.8*mean[:,2]/(0.001**2), color = 'navy')
	axs[0,1].fill_between(time, 0.95*0.8*(mean[:,2] - std[:,2])/(0.001**2), 0.95*0.8*(mean[:,2] + std[:,2])/(0.001**2), color = 'cornflowerblue')
	timeL = np.linspace(0,370, np.shape(linear[1])[0])
	axs[0,1].plot(timeL[:-2], linear[1][2:,2]/ep2, color = 'firebrick')


	## Warm Test ##
	stem = "evolutionW_005_18"
	mean, std = torqueMeanAndStd(stem, upperIndex = 4, file = "Nbody_Sormani_Data/Varying_Ep/")
	time = np.linspace(0,100, np.shape(mean)[0])
	axs[1,1].plot(time, mean[:,2]/ep2, color = 'navy')
	axs[1,1].fill_between(time, (mean[:,2] - std[:,2])/ep2,(mean[:,2] + std[:,2])/ep2, color = 'cornflowerblue')
	timeL = np.linspace(0,300, np.shape(linear[3])[0])
	axs[1,1].plot(timeL[:-2], linear[3][2:,2]/ep2, color = 'firebrick')

	



	axs[0,0].set_xlim([0,100])
	axs[0,1].set_xlim([0,100])
	axs[1,0].set_xlim([0,100])
	axs[1,1].set_xlim([0,100])

	axs[1,1].ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,4))
	
	xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	axs[0,0].set_ylabel("Self Consistent\n" + r"$\tau/\epsilon^{2}$", fontsize = 15)
	axs[1,0].set_ylabel("Test Particle\n" +r"$\tau/\epsilon^{2}$", fontsize = 15)
	axs[1,0].set_xticks(xlabelPositions)
	axs[1,0].set_xticklabels(xlabels)
	axs[1,1].set_xticks(xlabelPositions)
	axs[1,1].set_xticklabels(xlabels)
	axs[1,0].set_xlabel("Time", fontsize = 15)
	axs[1,1].set_xlabel("Time", fontsize = 15)

	axs[0,0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	axs[1,0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))

	axs[0,0].set_title(r"Cold", fontsize = 15)
	axs[0,1].set_title(r"Warm", fontsize = 15)
	
	#axs[1,0].legend(loc = 'upper right', fontsize = 12)
	

	plt.show()


def densityWakeDiffTemp(rMax = 5):
	twoDden = [TwoDdensity("Bar_Data/Sormani_Diff_Temp/evolution_Sormani_35.csv"), 
	TwoDdensity("Bar_Data/Sormani_Diff_Temp/evolutionTest_Sormani_35.csv"), 
	TwoDdensity("Bar_Data/Sormani_Diff_Temp/evolution_Sormani_45.csv"),
	TwoDdensity("Bar_Data/Sormani_Diff_Temp/evolutionTest_Sormani_45.csv")]
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(nrows = 4, ncols = 4, sharex = True, sharey = True )
	XX, YY = twoDden[0].meshgrid(rMax)
	time = [1, 2, 3, 4]
	print([np.amax((den[time[0]])) for den in twoDden])
	print([np.amax((den[time[1]])) for den in twoDden])
	print([np.amax((den[time[2]])) for den in twoDden])
	print([np.amax((den[time[3]])) for den in twoDden])
	max_values = max([np.amax((den[time[i]])) for i in range(0,4) for den in twoDden])
	for i in range(4):
		xBar, yBar = np.linspace(-sqrt(2)*rMax, sqrt(2)*rMax)*cos(2*(1/5.5) * time[i] * 10), np.linspace(-sqrt(2)*rMax, sqrt(2)*rMax)*sin((1/5.5) * time[i] * 10)
		axs[i, 0].set_ylabel(f"{(i+1)*4}" +r"$\times 2\pi$", fontsize = 15)
		for j in range(4):
			m = np.amax((twoDden[j][time[i]]))/max_values
			levels = [m*0.3, m*0.90]
			pcm1 = axs[i,j].contourf(XX, YY, (twoDden[j][time[i]])/max_values, levels  =50)
			axs[i,j].contour(XX, YY, (twoDden[j][time[i]])/max_values, levels  =levels, colors = 'k')
			axs[i,j].set(aspect = 1),
			axs[i,j].plot(xBar, yBar, color = 'firebrick', linestyle =':')
			axs[i,j].set_xlim([-rMax, rMax])
			axs[i,j].set_ylim([-rMax, rMax])

	cb_ax = fig.add_axes([0.92,.11,0.02,.76])
	cbar = fig.colorbar(pcm1, cax=cb_ax, orientation='vertical')
	cbar.set_ticks([-0.6,-0.3, 0, 0.3, 0.6])


	axs[0,0].set_title("Self Consistent", fontsize = 13)
	axs[0,1].set_title("Test Particle", fontsize = 13)
	axs[0,2].set_title("Self Consistent", fontsize = 13)
	axs[0,3].set_title("Test Particle", fontsize = 13)

	plt.figtext(0.30,0.93, "Cold", va="center", ha="center", fontsize=17)
	plt.figtext(0.71,0.93,"Warm", va="center", ha="center", fontsize=17)

	plt.show()



## Sormani Bar ##
## ----------- ##

def sormaniFits():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	dar = readingInRealCSV("Bar_Data/Sormani_Fitting/Sormani_Small.csv")
	plt.plot(dar[:,0], dar[:,1], label = r"$1.0$ [Kpc]", color = 'royalblue')

	dar = readingInRealCSV("Bar_Data/Sormani_Fitting/Sormani_Medium.csv")
	plt.plot(dar[:,0], dar[:,1], label = r"$1.5$ [Kpc]", color = 'cornflowerblue')

	dar = readingInRealCSV("Bar_Data/Sormani_Fitting/Sormani_Large.csv")
	plt.plot(dar[:,0], dar[:,1], label = r"$2.0$ [Kpc]", color = 'firebrick')

	dar = readingInRealCSV("Bar_Data/Sormani_Fitting/Sormani_Small_Fit.csv")
	plt.plot(dar[:,0], dar[:,1], linestyle = '--', color = 'royalblue')

	dar = readingInRealCSV("Bar_Data/Sormani_Fitting/Sormani_Medium_Fit.csv")
	plt.plot(dar[:,0], dar[:,1], linestyle = '--' , color = 'cornflowerblue')

	dar = readingInRealCSV("Bar_Data/Sormani_Fitting/Sormani_Large_Fit.csv")
	plt.plot(dar[:,0], dar[:,1], linestyle = '--', color = 'firebrick')

	plt.xlim([0,8])

	plt.ylabel(r"$\Phi_{2}(R)$", fontsize = 15)
	plt.xlabel(r"$R$ [Kpc]", fontsize = 15)
	plt.legend(title = r"$R_{q}$", fontsize = 15, title_fontsize = 15)
	plt.show()

def sormaniTorque():
	data = readingInRealCSV("KalnajsTorque/Sormani_Bar/evolution_Sormani_35.csv")
	plt.plot(data[:, 2], color = 'royalblue')

	data = readingInRealCSV("KalnajsTorque/Sormani_Bar/evolution_Sormani_45.csv")
	plt.plot(data[:, 2], color = 'firebrick')

	data = readingInRealCSV("KalnajsTorque/Sormani_Bar/evolutionTest_Sormani_35.csv")
	plt.plot(data[:, 2], color = 'royalblue', linestyle = '--')

	data = readingInRealCSV("KalnajsTorque/Sormani_Bar/evolutionTest_Sormani_45.csv")
	plt.plot(data[:, 2], color = 'firebrick', linestyle = '--')
	plt.xlim([0,100])
	plt.show()

def sormaniShape():
	# plt.rc('text', usetex=True)
	# plt.rc('font', family='serif')

	# fig, axs = plt.subplots(ncols=2, nrows=2, sharey = True)
	# gs = axs[0, 0].get_gridspec()
	# # remove the underlying axes
	# for ax in axs[0, :]:
	#     ax.remove()
	# axbig = fig.add_subplot(gs[0, :])
	

	# data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_Small.csv")
	# axbig.plot(data[:, 2], color = 'royalblue', label = "1.0 Kpc")
	# data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_test_Small.csv")
	# axbig.plot(data[:, 2], color = 'royalblue', linestyle = '--')
	
	# data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_density_Small.csv")
	# r, a, _ = data.fourierCoeffT(2, 5, 10)
	# axs[1,0].plot(r,a, color = 'royalblue')
	# data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_test_density_Small.csv")
	# r, a, _ = data.fourierCoeffT(2, 5, 10)
	# axs[1,1].plot(r,a, linestyle = '--', color = 'royalblue')
	

	# data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_Medium.csv")
	# axbig.plot(data[:, 2], color = 'cornflowerblue', label = "1.5 Kpc")
	# data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_test_Medium.csv")
	# axbig.plot(data[:, 2], color = 'cornflowerblue', linestyle = '--')
	
	# data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_density_Medium.csv")
	# r, a, _ = data.fourierCoeffT(2, 5, 10)
	# axs[1,0].plot(r/1.5,a, color = 'cornflowerblue')
	# data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_test_density_Medium.csv")
	# r, a, _ = data.fourierCoeffT(2, 5, 10)
	# axs[1,1].plot(r/1.5,a, linestyle = '--', color = 'cornflowerblue')


	# data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_Large.csv")
	# axbig.plot(data[:, 2], color = 'firebrick', label = "2.0 Kpc")
	# data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_test_Large.csv")
	# axbig.plot(data[:, 2], color = 'firebrick', linestyle = '--')
	
	# data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_density_Large.csv")
	# r, a, _ = data.fourierCoeffT(2, 5, 10)
	# axs[1,0].plot(r/2.0,a,color = 'firebrick')
	# data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_test_density_Large.csv")
	# r, a, _ = data.fourierCoeffT(2, 5, 10)
	# axs[1,1].plot(r/2.0,a, linestyle = '--', color = 'firebrick')




	# axbig.set_xlim([0,100])
	# axbig.set_ylabel(r"$\tau/\epsilon^{2}$", fontsize = 15)
	# axbig.set_xlabel(r"Time", fontsize = 15)

	# xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	# axbig.set_xticks(xlabelPositions)
	# axbig.set_xticklabels(xlabels)
	# axbig.legend(title = r"$R_{q}$", title_fontsize = 15, fontsize = 15, loc = 'lower right') # This might want to go on the other set of axis

	# axs[1,0].set_xlim([0,5])
	# axs[1,1].set_xlim([0,5])

	# axs[1,1].set_xlabel(r"$R/R_{q}$", fontsize = 15)
	# axs[1,1].set_title("Test Particle")
	# axs[1,0].set_title("Self Consistent")

	# axs[1,0].set_xlabel(r"$R/R_{q}$", fontsize = 15)
	# axs[1,0].set_ylabel(r"$\rho(R)/\epsilon$", fontsize = 15)

	# plt.show()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(ncols=2)

	data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_Small.csv")
	axs[0].plot(data[:, 2], color = 'royalblue', label = "1.0 Kpc")
	
	
	data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_density_Small.csv")
	r, a, _ = data.fourierCoeffT(2, 5, 10)
	axs[1].plot(r*0.18,a, color = 'royalblue', label = "1.0 Kpc")
	
	

	data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_Medium.csv")
	axs[0].plot(data[:, 2], color = 'cornflowerblue', label = "1.5 Kpc")
	
	
	data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_density_Medium.csv")
	r, a, _ = data.fourierCoeffT(2, 5, 10)
	axs[1].plot(r*0.18,a, color = 'cornflowerblue', label = "1.5 Kpc")
	


	data = readingInRealCSV("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_Large.csv")
	axs[0].plot(data[:, 2], color = 'firebrick', label = "2.0 Kpc")
	
	
	data = TwoDdensity("Bar_Data/Kalnajs_Shape/Kalnajs_consistent_density_Large.csv")
	r, a, _ = data.fourierCoeffT(2, 5, 10)
	axs[1].plot(r*0.18,a,color = 'firebrick', label = "2.0 Kpc")
	

	axs[0].set_xlim([0,100])
	axs[0].set_ylabel(r"$\tau/\epsilon^{2}$", fontsize = 15)
	axs[0].set_xlabel(r"Time", fontsize = 15)

	xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)
	axs[1].legend(title = r"$R_{q}$", title_fontsize = 15, fontsize = 15, loc = 'upper right') # This might want to go on the other set of axis

	axs[1].set_xlim([0,1])
	axs[1].set_ylabel(r"$\rho(R)/\epsilon$", fontsize = 15)
	axs[1].set_xlabel(r"$R/R_{CR}$", fontsize = 15)

	plt.show()

def sormaniSpeed():
	# plt.rc('text', usetex=True)
	# plt.rc('font', family='serif')
	# filename_C, filename_T = [f"Bar_Data/Pattern_Speed/Kalanajs_evolution_consistent_{i}.csv" for i in range(0, 25, 6)], [f"Bar_Data/Pattern_Speed/Kalanajs_evolution_test_{i}.csv" for i in range(0, 25, 6)]
	# evolution_C, evolution_T = [readingInRealCSV(file) for file in filename_C], [readingInRealCSV(file) for file in filename_T]

	# filename_C, filename_T = [f"Bar_Data/Pattern_Speed/Kalanajs_density_consistent_{i}.csv" for i in range(0, 25, 6)], [f"Bar_Data/Pattern_Speed/Kalanajs_density_test_{i}.csv" for i in range(0, 25, 6)]
	# density_C, density_T = [TwoDdensity(file) for file in filename_C], [TwoDdensity(file) for file in filename_T]

	# fig, axs = plt.subplots(ncols = 2, nrows = 2, sharey = 'row')
	# cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-0.0, vmax=0.30))

	# for c, t, dc, dt, op, in zip(evolution_C, evolution_T, density_C, density_T, [0, 0.06, 0.12, 0.18, 0.24]):
	# 	axs[0,0].plot(c[:,2]/(1), color =  cmap.to_rgba(op))
	# 	r, a, _ = dc.fourierCoeffT(2, 5, 5)
	# 	axs[1,0].plot(r, a, color =  cmap.to_rgba(op))
		
	# 	axs[0,1].plot(t[:,2]/(1), label = f"{op:.2f}", color =  cmap.to_rgba(op))
	# 	r, a, _ = dt.fourierCoeffT(2, 5, 5)
	# 	axs[1,1].plot(r,a, label = f"{op:.2f}",color =  cmap.to_rgba(op))

	# axs[0,0].set_xlim([0,100])
	# axs[0,1].set_xlim([0,100])

	# axs[0,0].set_ylabel(r"$\tau/\epsilon^{2}$", fontsize = 15)
	# axs[0,0].set_xlabel(r"Time", fontsize = 15)
	# axs[0,1].set_xlabel(r"Time", fontsize = 15)
	# axs[0,0].set_title("Self Consistent", fontsize =15)
	# axs[0,1].set_title("Test Particle", fontsize =15)


	# xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	# axs[0,0].set_xticks(xlabelPositions)
	# axs[0,0].set_xticklabels(xlabels)
	# axs[0,1].set_xticks(xlabelPositions)
	# axs[0,1].set_xticklabels(xlabels)

	# axs[1,1].legend(title = r"$\Omega_{p}$", title_fontsize = 12, fontsize = 12, loc = 'lower right')


	# plt.show()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	filename_C, filename_T = [f"Bar_Data/Pattern_Speed/Kalanajs_evolution_consistent_{i}.csv" for i in range(0, 25, 6)], [f"Bar_Data/Pattern_Speed/Kalanajs_evolution_test_{i}.csv" for i in range(0, 25, 6)]
	evolution_C, evolution_T = [readingInRealCSV(file) for file in filename_C], [readingInRealCSV(file) for file in filename_T]

	filename_C, filename_T = [f"Bar_Data/Pattern_Speed/Kalanajs_density_consistent_{i}.csv" for i in range(0, 25, 6)], [f"Bar_Data/Pattern_Speed/Kalanajs_density_test_{i}.csv" for i in range(0, 25, 6)]
	density_C, density_T = [TwoDdensity(file) for file in filename_C], [TwoDdensity(file) for file in filename_T]

	fig, axs = plt.subplots(ncols = 2)
	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-0.0, vmax=0.30))

	for c, t, dc, dt, op, in zip(evolution_C, evolution_T, density_C, density_T, [0, 0.06, 0.12, 0.18, 0.24]):
		axs[0].plot(c[:,2]/(1), color =  cmap.to_rgba(op))
		
		r, a, _ = dc.fourierCoeffT(2, 5, 10)
		axs[1].plot(op*r, a, color =  cmap.to_rgba(op), label = f"{op:.2f}")
		

	axs[0].set_xlim([0,100])


	axs[0].set_ylabel(r"$\tau/\epsilon^{2}$", fontsize = 15)
	axs[0].set_xlabel(r"Time", fontsize = 15)
	axs[0].set_title("Torque", fontsize =15)
	


	xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)

	axs[1].legend(title = r"$\Omega_{p}$", title_fontsize = 12, fontsize = 12, loc = 'upper right')
	axs[1].set_xlim([0,1])
	axs[1].set_xlabel(r"$R/R_{CR} = R \Omega_{p}/v_{c}$", fontsize = 15)
	axs[1].set_ylabel(r"$\rho(R)/\epsilon$", fontsize = 15)
	axs[1].set_title("Density", fontsize = 15)


	plt.show()

def sormaniRatio():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(ncols=2)

	data = readingInRealCSV("Bar_Data/Sormani_Ratio/Kalnajs_consistent_Small.csv")
	axs[0].plot(data[:, 2], color = 'royalblue', label = "1.0 Kpc")
	
	
	data = TwoDdensity("Bar_Data/Sormani_Ratio/Kalnajs_consistent_density_Small.csv")
	r, a, _ = data.fourierCoeffT(2, 5, 10)
	axs[1].plot(r*0.27,a, color = 'royalblue', label = "1.0 Kpc")
	
	

	data = readingInRealCSV("Bar_Data/Sormani_Ratio/Kalnajs_consistent_Medium.csv")
	axs[0].plot(data[:, 2], color = 'cornflowerblue', label = "1.5 Kpc")
	
	
	data = TwoDdensity("Bar_Data/Sormani_Ratio/Kalnajs_consistent_density_Medium.csv")
	r, a, _ = data.fourierCoeffT(2, 5, 10)
	axs[1].plot(r*0.18,a, color = 'cornflowerblue', label = "1.5 Kpc")
	


	data = readingInRealCSV("Bar_Data/Sormani_Ratio/Kalnajs_consistent_Large.csv")
	axs[0].plot(data[:, 2], color = 'firebrick', label = "2.0 Kpc")
	
	
	data = TwoDdensity("Bar_Data/Sormani_Ratio/Kalnajs_consistent_density_Large.csv")
	r, a, _ = data.fourierCoeffT(2, 5, 10)
	axs[1].plot(r*0.135,a,color = 'firebrick', label = "2.0 Kpc")
	

	axs[0].set_xlim([0,100])
	axs[0].set_ylabel(r"$\tau/\epsilon^{2}$", fontsize = 15)
	axs[0].set_xlabel(r"Time", fontsize = 15)

	xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)
	axs[1].legend(title = r"$R_{q}$", title_fontsize = 15, fontsize = 15, loc = 'upper right') # This might want to go on the other set of axis

	axs[1].set_xlim([0,1])
	axs[1].set_ylabel(r"$\rho(R)/\epsilon$", fontsize = 15)
	axs[1].set_xlabel(r"$R/R_{CR}$", fontsize = 15)

	plt.show()

def sormaniRatioDensity():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, (axs, cbax) = plt.subplots(nrows = 2, ncols = 3, gridspec_kw={"height_ratios":[1, 0.03]})
	ILR = 1 - 0.5 *sqrt(2)
	

	data = TwoDdensity("Bar_Data/Sormani_Ratio/Kalnajs_consistent_density_Small.csv")
	XX, YY = data.meshgrid(5)
	m = np.amax(data[5])
	print(m)
	axs[0].contour(XX, YY, data[5], levels = [m*0.3, m *0.9], colors = 'black')
	im = axs[0].contourf(XX, YY, data[5], levels = 100)
	axs[0].plot([3.7 * ILR * cos(t) for t in np.linspace(0,2*pi)], [3.7 * ILR * sin(t) for t in np.linspace(0,2*pi)], color = 'firebrick', linestyle = '--')
	axs[0].set(aspect = 1)
	axs[0].set_xlim([-5,5])
	axs[0].set_ylim([-5,5])
	axs[0].plot(sqrt(2)*5 *np.linspace(-1,1) * cos(2 * 10 * 100 * 0.27), sqrt(2)*5 *np.linspace(-1,1) * sin(2 * 10 * 100 * 0.27), color = 'firebrick', linestyle = ':')
	axs[0].set_title(r"$R_{q} = 1.0$", fontsize = 15)
	cb = fig.colorbar(im, cax=cbax[0], orientation="horizontal")
	cb.set_ticks([-6.0, -3.0, 0, 3, 6])
	cbax[0].set_xlabel(r"$\psi(R,\phi)/\epsilon$", fontsize = 15)
	
	data = TwoDdensity("Bar_Data/Sormani_Ratio/Kalnajs_consistent_density_Medium.csv")
	XX, YY = data.meshgrid(5)
	m = np.amax(data[5])
	print(m)
	axs[1].contour(XX, YY, data[5], levels = [m*0.3, m *0.9], colors = 'black')
	im = axs[1].contourf(XX, YY, data[5], levels = 100)
	axs[1].plot([5.5 * ILR * cos(t) for t in np.linspace(0,2*pi)], [5.5 * ILR * sin(t) for t in np.linspace(0,2*pi)], color = 'firebrick', linestyle = '--')
	axs[1].set(aspect = 1)
	axs[1].set_xlim([-5,5])
	axs[1].set_ylim([-5,5])
	axs[1].plot(sqrt(2)*5 *np.linspace(-1,1) * cos(2 * 10 * 100 * 0.18), sqrt(2)*5 *np.linspace(-1,1) * sin(2 * 10 * 100 * 0.18), color = 'firebrick', linestyle = ':')
	axs[1].set_title(r"$R_{q} = 1.5$", fontsize = 15)
	cb = fig.colorbar(im, cax=cbax[1], orientation="horizontal")
	cb.set_ticks([-6.0, -3.0, 0, 3, 6])
	cbax[1].set_xlabel(r"$\psi(R,\phi)/\epsilon$", fontsize = 15)


	data = TwoDdensity("Bar_Data/Sormani_Ratio/Kalnajs_consistent_density_Large.csv")
	XX, YY = data.meshgrid(5)
	m = np.amax(data[5])
	print(m)
	axs[2].contour(XX, YY, data[5], levels = [m*0.3, m *0.9], colors = 'black')
	im = axs[2].contourf(XX, YY, data[5], levels = 100)
	axs[2].plot([7.4 * ILR * cos(t) for t in np.linspace(0,2*pi)], [7.4 * ILR * sin(t) for t in np.linspace(0,2*pi)], color = 'firebrick', linestyle = '--')
	axs[2].plot(5*sqrt(2))
	axs[2].set(aspect = 1)
	axs[2].set_xlim([-5,5])
	axs[2].set_ylim([-5,5])
	axs[2].plot(sqrt(2)*5 *np.linspace(-1,1) * cos(2 * 10 * 100 * 0.135), sqrt(2)*5 *np.linspace(-1,1) * sin(2 * 10 * 100 * 0.135), color = 'firebrick', linestyle = ':')
	axs[2].set_title(r"$R_{q} = 2.0$", fontsize = 15)
	cb = fig.colorbar(im, cax=cbax[2], orientation="horizontal")
	cb.set_ticks([-6.0, -3.0, 0, 3, 6])
	cbax[2].set_xlabel(r"$\psi(R,\phi)/\epsilon$", fontsize = 15)

	plt.show()



def coefficentEvolution(): # Maybe Could include a shaded region for the amplitude? 
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	linear = readingInComplexCSV("Nbody_Sormani_Data/Varying_Ep/Linear_Test_Particle_Coeff.csv")

	fig,axs = plt.subplots(ncols = 2)
	lst_n = range(0,15,3)
	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-1, vmax=17))

	for n in lst_n:
		axs[0].plot(np.linspace(0, 100, np.shape(linear[1:100,n])[0]), np.absolute(linear[1:100,n]), color = cmap.to_rgba(n))
		deriv = getGradient(np.linspace(0, 100, np.shape(linear[1:100,n])[0]), np.angle(linear[1:100,n]), multiple = 3, includeMinus = 1)

		axs[1].plot(*deriv, label = f"{n}", color = cmap.to_rgba(n))

	axs[0].set_ylabel(r"$|B_{p}(t)|$", fontsize = 15)
	axs[1].set_ylabel(r"$\partial\arg \left(B_{p}\right)/\partial t$", fontsize = 15)

	axs[1].axhline(-2*0.177, color = 'k', linestyle = '--')
	axs[1].set_ylim([-.6,-.2])
	axs[1].legend(title = r"Index, $p$", title_fontsize = 15, fontsize =12)

	axs[0].set_title("Amplitude", fontsize = 15)
	axs[1].set_title("Phase", fontsize = 15)

	xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)
	axs[1].set_xticks(xlabelPositions)
	axs[1].set_xticklabels(xlabels)
	axs[0].set_xlabel("Time", fontsize = 15)
	axs[1].set_xlabel("Time", fontsize = 15)
	
	axs[0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	plt.show()







#densityWakeDiffTemp()
#kalnajsNbodyTorque()
#torquePlots()
#sormaniFits()
#sormaniTorque()
#sormaniShape()
#sormaniRatioDensity()

#sormaniSpeed()

#coefficentEvolution()
sormaniRatioDensity() 