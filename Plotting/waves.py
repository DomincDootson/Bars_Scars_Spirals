from generalPlottingFunctions import *
from Density_Classes.OneDdensity import *
from Density_Classes.TwoDdensity import *
from Density_Classes.WaveFitter import * 
from CoefficientClass import *

import matplotlib.animation as animation


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



def wavesFittingTest():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(nrows = 4, ncols = 3, sharex = 'col', sharey = 'col')

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
		axs[i,1].axhline(0, linestyle = '--', color = 'k')

		#axs[i,1].axvline(rad[i], linestyle = '--', color = 'black')


		XX, YY = wave.density.meshgrid(rMax = 15)
		s = slice(40,160)
		axs[i,2].contourf(XX[s, s], YY[s, s], wave.density[time][s, s], levels = 100)
		m = np.amax(wave.density[time][s, s])
		axs[i,2].contour(XX[s, s], YY[s, s], wave.density[time][s, s], levels = [0.5 *m], colors = 'black')
		axs[i,2].set(aspect = 1)
		#axs[i,2].plot([rad[i]*cos(th) for th in np.linspace(0,2*pi)], [rad[i]*sin(th) for th in np.linspace(0,2*pi)], color = 'firebrick', linestyle = '--')

		axs[i,0].set_xlim([1, 10])
		axs[i,1].set_xlim([1, 10])
		

	
	axs[-1,0].set_xlabel(r"$R$", fontsize  = 15)
	axs[-1,1].set_xlabel(r"$R$", fontsize  = 15)
	axs[-1,2].set_xlabel(r"$x$", fontsize  = 15)

	axs[0,0].set_title(r"$\rho(R,t)/M$", fontsize = 15)
	axs[0,1].set_title(r"$k(R,t)$", fontsize = 15)
	axs[0,2].set_title(r"$\rho(R,\phi)$", fontsize = 15)

	axs[0,0].set_ylabel(r"$t = 0$", fontsize = 15)
	axs[1,0].set_ylabel(r"$t = 2 \times 2\pi$", fontsize = 15)
	axs[2,0].set_ylabel(r"$t = 4 \times 2\pi$", fontsize = 15)
	axs[3,0].set_ylabel(r"$t = 6 \times 2\pi$", fontsize = 15)

	axs[0,1].legend(fontsize = 10)


	plt.show()

def wavesDirectionTest():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(nrows = 4, ncols = 2, sharex = True)

	waves = WaveFitter("Waves_Data/test_spiral.csv", timeEnd = 130*0.3, l = 2, rmax = 15)


	#waves.split_waves()
	waves.split_waves_k(1.297)
	times = [45, 60, 75, 90]
	for i in range(4):
		time = times[i]
		
		i_s, i_l, o_s, o_l = waves.ingoing_S(time), waves.ingoing_L(time), waves.outgoing_S(time), waves.outgoing_L(time)

		
		## Outgoing ##
		axs[i,0].plot(waves.radii, np.real(o_s), color = 'firebrick', label = r"Short")
		axs[i,0].plot(waves.radii, np.absolute(o_s), color = 'firebrick', linestyle = "--")

		axs[i,0].plot(waves.radii, np.real(o_l), color = 'royalblue', label = r"Long")
		axs[i,0].plot(waves.radii, np.absolute(o_l), color = 'royalblue', linestyle = "--")
		axs[i,0].set_xlim([waves.radii[8],8])
		
		axs[i,1].plot(*getGradient(waves.radii, np.angle(o_s), multiple = 5), color = 'firebrick', label = r"Short")
		axs[i,1].plot(*getGradient(waves.radii, np.angle(o_l), multiple = 5), color = 'royalblue', label = r"Long")
		axs[i,1].axhline(0, color = 'k', linestyle = '--')
		axs[i,1].set_xlim([waves.radii[8],8])

		## Ingoing ##
		# axs[i,2].plot(waves.radii, np.real(i_s), color = 'firebrick', label = r"Short")
		# axs[i,2].plot(waves.radii, np.absolute(i_s), color = 'firebrick', linestyle = "--")

		# axs[i,2].plot(waves.radii, np.real(i_l), color = 'royalblue', label = r"Long")
		# axs[i,2].plot(waves.radii, np.absolute(i_l), color = 'royalblue', linestyle = "--")
		# axs[i,2].set_xlim([waves.radii[5],8])
		
		# axs[i,3].plot(*getGradient(waves.radii, np.angle(i_s), multiple = 5), color = 'firebrick', label = r"Short")
		# axs[i,3].plot(*getGradient(waves.radii, np.angle(i_l), multiple = 5), color = 'royalblue', label = r"Long")
		# axs[i,3].axhline(0, color = 'k', linestyle = '--')
		# axs[i,3].set_xlim([waves.radii[5],8])
		


	
	axs[-1,0].set_xlabel(r"$R$", fontsize  = 15)
	axs[-1,1].set_xlabel(r"$R$", fontsize  = 15)
	# axs[-1,2].set_xlabel(r"$R$", fontsize  = 15)
	# axs[-1,3].set_xlabel(r"$R$", fontsize  = 15)

	axs[0,0].set_title("Out-going\nDensity", fontsize = 15)
	axs[0,1].set_title("Out-going\nWavenumber", fontsize = 15)
	# axs[0,2].set_title("In-going\nAmplitude", fontsize = 15)
	# axs[0,3].set_title("In-going\nPhase", fontsize = 15)

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


wavesFittingTest()
#wavesDirectionTest()
# wave = WaveFitter("Waves_Data/test_spiral.csv", timeEnd = 40, l = 2, rmax = 15)
# wave.check_k_splitting(50, 1.297)


#wave.check_k_fitting(48)
# wave.check_k_splitting(48, 1.09)
#wave.k_splitting_Plot(48, 1.09)
#wave.splitting_Plot(45)
#plt.show()
#wave.check_k_fitting(48)

#wavesFittingTest()

''' Run this '''

# wave = WaveFitter("Test_density.csv", timeEnd = 100, l = 2, rmax = 10)
# wave.splitting_Plot()
