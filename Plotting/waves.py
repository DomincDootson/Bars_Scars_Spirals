from generalPlottingFunctions import *
from Density_Classes.OneDdensity import *
from Density_Classes.TwoDdensity import *
from Density_Classes.WaveFitter import * 
from CoefficientClass import *

import matplotlib.animation as animation
from statistics import median

def modesComparison(oneDdensity, filename = None, index2Normalise = -1):
	modePower = [each.modePower() for each in oneDdensity]

	if index2Normalise != -1:
		for each in modePower:
			each[:, index2Normalise] -= each[0,index2Normalise]

	bins, minValue, maxValue = [range(each.nMax) for each in oneDdensity], min([np.amin(each) for each in modePower]), max([np.max(each) for each in modePower])

	fig, axs = plt.subplots()
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=20, metadata=dict(artist='Me'))

	color, alpha = ['firebrick', "royalblue"], [1, 0.8]


	def animate(time):
		axs.clear()
		axs.set_title("Time = " + str(time*0.2))
		axs.set_ylim([minValue, maxValue])
		
		for i in range(len(modePower)):
			axs.bar(bins[i], modePower[i][time, :], color = color[i], alpha = alpha[i], label =len(bins[i])-1)

		axs.legend()
		return axs

	ani = animation.FuncAnimation(fig, animate, frames = 250)
	if (filename):
		ani.save(filename, writer = writer)
		print("Animation saved to: " + filename) 
	else:
		plt.show()


def savingSplitWaves():
	fitter = WaveFitter("Waves_Data/Stirring/CR_5_OLR_W_N.csv", 10)
	fitter.split_waves()
	fitter.density_animation("Waves_Plots/Waves_Videos/Stirring/CR_5_OLR_W_N.mp4")

	fitter = WaveFitter("Waves_Data/Stirring/CR_5_OLR_W_P.csv", 10)
	fitter.split_waves()
	fitter.density_animation("Waves_Plots/Waves_Videos/Stirring/CR_5_OLR_W_P.mp4")

	fitter = WaveFitter("Waves_Data/Stirring/CR_5_OLR_R_N.csv", 10)
	fitter.split_waves()
	fitter.density_animation("Waves_Plots/Waves_Videos/Stirring/CR_5_OLR_R_N.mp4")

	fitter = WaveFitter("Waves_Data/Stirring/CR_5_OLR_R_P.csv", 10)
	fitter.split_waves()
	fitter.density_animation("Waves_Plots/Waves_Videos/Stirring/CR_5_OLR_R_P.mp4")


def decomposed1Danimations(readStem, writeStem, file2read):
	for f2r in file2read:
		wavefitter = WaveFitter(readStem + f2r + '.csv', 10)
		wavefitter.split_waves()
		
		write2 = writeStem + f2r +".mp4"
		wavefitter.density_animation(filename = write2)

def getGradient(radius, phase, multiple = 8):
	rad, grad = [r1 for r1 in radius[1:]], [(ph1-ph0)/(r1-r0) for ph1, ph0, r1, r0 in zip(phase[1:], phase, radius[1:], radius)]
	r_o, g_o =[], []
	med = median(grad)
	for r, g in zip(rad, grad):
		if abs(g) < multiple * med:
			r_o.append(r)
			g_o.append(-g)


	return r_o, g_o

def wavesFittingTest():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(nrows = 4, ncols = 3, sharex = 'col')

	wave = WaveFitter("Waves_Data/test_spiral.csv", timeEnd = 130*0.3, l = 2, rmax = 15)


	wave.split_waves()
	times, rad = [0,42,84,126], [0, 3.53, 5.01, 6.14]
	for i in range(4):
		time = times[i]
		axs[i,0].plot(wave.radii, wave.amp[time, :], color = 'firebrick', linestyle = '--')
		axs[i,0].plot(wave.radii, wave.amp[time, :] * np.cos(wave.phase[time, :]), color = 'firebrick')
		axs[i,1].plot(*getGradient(wave.radii, wave.phase[time, :]), color = 'firebrick', label = 'Original')
		axs[i,0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	

		amp, phase = wave.total_T(time)
		axs[i,0].plot(wave.radii, amp, color = 'royalblue', linestyle = '--')
		axs[i,0].plot(wave.radii, amp * np.cos(phase), color = 'royalblue')
		axs[i,1].plot(*getGradient(wave.radii, phase), color = 'royalblue', label = "Reconstructed")
		axs[i,1].axhline(0, linestyle = '--')

		axs[i,1].axvline(rad[i], linestyle = '--', color = 'black')


		XX, YY = wave.density.meshgrid(rMax = 15)
		axs[i,2].contourf(XX, YY, wave.density[time], levels = 100)
		m = np.amax(wave.density[time])
		axs[i,2].contour(XX, YY, wave.density[time], levels = [0.5 *m], colors = 'firebrick')
		axs[i,2].set(aspect = 1)
		axs[i,2].plot([rad[i]*cos(th) for th in np.linspace(0,2*pi)], [rad[i]*sin(th) for th in np.linspace(0,2*pi)], color = 'black', linestyle = '--')

		axs[i,0].set_xlim([1, 10])
		axs[i,1].set_xlim([1, 10])
		

	
	axs[-1,0].set_xlabel(r"$R$", fontsize  = 15)
	axs[-1,1].set_xlabel(r"$R$", fontsize  = 15)
	axs[-1,2].set_xlabel(r"$x$", fontsize  = 15)

	axs[0,0].set_title(r"$\rho(R,t)/M$", fontsize = 15)
	axs[0,1].set_title(r"$k(R,t)$", fontsize = 15)
	axs[0,2].set_title(r"Density", fontsize = 15)

	axs[0,0].set_ylabel(r"$t = 0$", fontsize = 15)
	axs[1,0].set_ylabel(r"$t = 2 \times 2\pi$", fontsize = 15)
	axs[2,0].set_ylabel(r"$t = 4 \times 2\pi$", fontsize = 15)
	axs[3,0].set_ylabel(r"$t = 6 \times 2\pi$", fontsize = 15)

	axs[0,1].legend(fontsize = 10)


	plt.show()

def wavesDirectionTest():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(nrows = 4, ncols = 3, sharex = True, sharey = 'col')

	waves = WaveFitter("Waves_Data/test_spiral.csv", timeEnd = 130*0.3, l = 2, rmax = 15)


	waves.split_waves()
	times = [45, 60, 75, 90]
	for i in range(4):
		time = times[i]
		
		ingoing, outgoing = waves.ingoing_T(time), waves.outgoing_T(time)

		axs[i,0].plot(waves.radii, np.real(ingoing), color = 'royalblue', label = r"$\rho(R,\phi=0,t)$")
		axs[i,0].plot(waves.radii, np.absolute(ingoing), color = 'firebrick', linestyle = "--", label =  r"$\rho(R,t)$")
		axs[i,0].set_xlim([waves.radii[5],8])

		axs[i,1].plot(waves.radii, np.real(outgoing), color = 'royalblue')
		axs[i,1].plot(waves.radii, np.absolute(outgoing), color = 'firebrick', linestyle = "--")
		axs[i,1].set_xlim([waves.radii[5],8])

		axs[i,2].plot(waves.radii, np.real(ingoing+outgoing), color = 'royalblue')
		axs[i,2].plot(waves.radii, np.absolute(ingoing+outgoing), color = 'firebrick', linestyle = "--")
		axs[i,2].set_xlim([waves.radii[5],8])


	
	axs[-1,0].set_xlabel(r"$R$", fontsize  = 15)
	axs[-1,1].set_xlabel(r"$R$", fontsize  = 15)
	axs[-1,2].set_xlabel(r"$R$", fontsize  = 15)

	axs[0,0].set_title(r"In-going", fontsize = 15)
	axs[0,1].set_title(r"Out-going", fontsize = 15)
	axs[0,2].set_title(r"Total", fontsize = 15)

	axs[0,0].set_ylabel(r"$t = 2.1 \times 2\pi$", fontsize = 15)
	axs[1,0].set_ylabel(r"$t = 2.9 \times 2\pi$", fontsize = 15)
	axs[2,0].set_ylabel(r"$t = 3.6 \times 2\pi$", fontsize = 15)
	axs[3,0].set_ylabel(r"$t = 4.3 \times 2\pi$", fontsize = 15)

	axs[0,0].legend(fontsize = 10)


	plt.show()
#saveAnimations()

#decomposed1Danimations("Waves_Data/Stirring/", "Waves_Plots/Waves_Videos/Stirring/", [f"CR_5_OLR_{w}_{d}" for w, d  in zip(['W', 'W', 'R', 'R'], ['P', 'N', 'P', 'N'])])
#twoD.fourierAnimations(2, rMax = 10, filename ="Another_vid_4_ali.mp4")

#den = TwoDdensity("Waves_Data/Stirring/CR_5_OLR_W_N.csv")


#savingSplitWaves()
# density = TwoDdensity("test_spiral.csv")
# plt.plot(np.linspace(-15,15,201),density[15][100,:])

#wavesFittingTest()

wavesDirectionTest()
# wave = WaveFitter("Waves_Data/test_spiral.csv", timeEnd = 40, l = 2, rmax = 15)
# wave.split_waves()
# wave.density_animation("spiral_video.mp4")




