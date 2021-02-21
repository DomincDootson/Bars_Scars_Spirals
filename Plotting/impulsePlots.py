import matplotlib.animation as animation

from generalPlottingFunctions import *


def impulseDensityName(littleSigma, radius, angHarmonic):
	return "../Disk_Kicking/littleSigma_" + str(round(littleSigma*100)) + "/Density" + str(round(radius*10)) +"_" + str(angHarmonic) + ".csv" 

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


densityAnimation(0.25, 1)
