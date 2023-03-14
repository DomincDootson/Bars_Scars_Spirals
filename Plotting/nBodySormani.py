from generalPlottingFunctions import *
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm

def sormaniBackwardsFractionOmegaP():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	data = readingInRealCSV("Nbody_Sormani_Data/BO_omegaP_cold_large.csv")

	op, d = data[:,0],data[:,1:]
	mean, s = np.mean(d, axis = (1)), np.std(d, axis = (1))
	plt.errorbar(op, mean, yerr = s, color = 'royalblue', label = "0.35")

	data = readingInRealCSV("Nbody_Sormani_Data/BO_omegaP_warm_large.csv")

	op, d = data[:,0],data[:,1:]
	mean, s = np.mean(d, axis = (1)), np.std(d, axis = (1))
	plt.errorbar(op, mean, yerr = s, color = 'firebrick', label = "0.45")

	plt.ylabel("Reverse Fraction", fontsize = 15)
	plt.xlabel(r"$\Omega_{p}$", fontsize = 15)

	plt.legend(title = r"$\sigma_{R}/v_{c}$", title_fontsize =15, fontsize =15)

	plt.show()


def sormaniBackwardsFractionAmplitude():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	data = readingInRealCSV("Nbody_Sormani_Data/BO_Amp_realOmega_cold_large.csv")

	op, d = data[:,0],data[:,1:]
	mean, s = np.mean(d, axis = (1)), np.std(d, axis = (1))
	plt.errorbar(op, mean, yerr = s, color = 'royalblue', label = "0.35")

	data = readingInRealCSV("Nbody_Sormani_Data/BO_Amp_realOmega_warm_large.csv")

	op, d = data[:,0],data[:,1:]
	mean, s = np.mean(d, axis = (1)), np.std(d, axis = (1))
	plt.errorbar(op, mean, yerr = s, color = 'firebrick', label = "0.45")

	plt.ylabel("Reverse Fraction", fontsize = 15)
	plt.xlabel(r"$\epsilon$", fontsize = 15)

	plt.legend(title = r"$\sigma_{R}/v_{c}$", title_fontsize =15, fontsize =15)

	plt.show()

def angularMomentum(evolution):
	return evolution[:, 1] * evolution[:, 4] - evolution[:, 2] *  evolution[:, 3] 

def particleEvolution():
	data = readingInRealCSV("BoxOrbit_1.csv")
	cmap = ScalarMappable(cmap = 'plasma', norm = Normalize(vmin=-1.5, vmax=2))

	am = angularMomentum(data)
	am = am/ np.absolute(am)

	colors = cmap.to_rgba(am)
	
	for i in range(np.shape(data)[0]-1):
		color = 'firebrick' if am[i]<0 else 'royalblue'
		plt.plot([data[i, 1],data[i+1, 1]], [data[i, 2],data[i+1, 2]], color = color)

	plt.xlabel("x")
	plt.ylabel("y")

	plt.scatter([0],[0], color = 'firebrick')
	plt.show()

def phaseSpacePlots(filename, omegaP, ep):
	data = readingInRealCSV(filename)
	fig, axs = plt.subplots(ncols = 2)

	axs[0].hist2d(data[:,1],data[:,0], bins = 500)
	axs[1].hist2d(data[:,3],data[:,2], bins = 500)

	axs[0].set_xlabel("Angular Momentum")
	axs[1].set_xlabel("Angular Momentum Corotating")

	axs[0].set_ylabel("Energy")
	axs[1].set_ylabel("Jacobi Integral")
	axs[0].set_ylim([-0.5, 4])

	fig.suptitle(f"Pattern Speed {omegaP} epsilon {ep}")

	plt.show()


def changePhaseSpace(files, omegaP, epsilon):
	# plt.rc('text', usetex=True)
	# plt.rc('font', family='serif')
	# data = [readingInRealCSV(file) for file in files]
	# fig, axs = plt.subplots(ncols = len(data), nrows = 2, sharey = 'row', sharex = 'row')

	# # xlim = [np.amin(data[0][:,1]), np.amax(data[0][:,1])]
	# # ylim = [np.amin(data[0][:,0]), np.amax(data[0][:,0])]
	# ylim = [0, 4.55]
	# xlim = [0, 2.5]
	# jacobi =  1.28
	# for ax, d in zip(axs[0,:], data):
		
	# 	ax.hist2d(d[:,0],d[:,1], bins = [np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/100), np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/100)])
	# 	ax.set_xlim([0, 2.45])
	# 	ax.set_ylim([0, 4.5])
	# 	ax.set_xlabel(r"$H(t)$", fontsize = 15)
	# 	#ax.plot([0.09*L*2 + jacobi for L in np.linspace(0,4.5)], np.linspace(0,4.5), color = 'firebrick')
	
	# ylim = [0, 5]
	# xlim = [0, 3]
	# for ax, d in zip(axs[1,:], data):
		
	# 	ax.hist2d(d[:,2],d[:,1], bins = [np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/100), np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/100)])
	# 	ax.axvline(jacobi, color = 'firebrick')
	# 	ax.set_xlim([0, 2.9])
	# 	ax.set_ylim([0, 4.5])
	# 	#ax.set_xlabel(r"$H_{J}(t)$", fontsize = 15)

	# axs[0,0].set_ylabel(r"$L_{z}$", fontsize = 15)
	# axs[1,0].set_ylabel(r"$L_{z}$", fontsize = 15)
	# axs[0,0].set_title(r"$t = 0$", fontsize = 15)
	# axs[0,1].set_title(r"$t = 7.5 \times 2\pi$", fontsize = 15)
	# axs[0,2].set_title(r"$t = 15.0 \times 2\pi$", fontsize = 15)

	# plt.show()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	data = [readingInRealCSV(file) for file in files]
	fig, axs = plt.subplots(ncols = len(data), sharey = 'row', sharex = 'row')

	
	
	ylim = [0, 4.55]
	xlim = [0, 2.5]
	jacobi =  1.28
	for ax, d in zip(axs, data):
		
		ax.hist2d(d[:,0],d[:,1], bins = [np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/100), np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/100)])
		ax.set_xlim([0, 2.45])
		ax.set_ylim([0, 4.5])
		ax.set_xlabel(r"$H(t)$", fontsize = 15)
		#ax.plot([0.09*L*2 + jacobi for L in np.linspace(0,4.5)], np.linspace(0,4.5), color = 'firebrick')
	

	axs[0].set_ylabel(r"$L_{z}$", fontsize = 15)
	
	axs[0].set_title(r"$t = 0$", fontsize = 15)
	axs[1].set_title(r"$t = 7.5 \times 2\pi$", fontsize = 15)
	axs[2].set_title(r"$t = 15.0 \times 2\pi$", fontsize = 15)

	plt.show()
	

def varyingEpsilon(linearEp = 0.005): 
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(ncols = 2, sharey = True)


	epS, epD, epC = ['005', '01', '05', '1'], [0.005, 0.01, 0.05, 0.1], [0.005, 0.01, 0.05, 0.1]
	#epS, epD = ['005'], [0.005]
	cmap = ScalarMappable(cmap = 'viridis_r', norm = LogNorm(vmin=0.004, vmax=0.25))

	for s, d, c in zip(epS, epD, epC):
		stem = "evolution_" + s + "_18"
		mean, std = torqueMeanAndStd(stem, upperIndex = 4, file = "Nbody_Sormani_Data/Varying_Ep/")
		time = np.linspace(0,100, np.shape(mean)[0])
		axs[0].plot(time, mean[:,2]/(d**2), color = cmap.to_rgba(c), label =r"$\log_{10} \epsilon = $ " f"{log(d,10):.2f}")
		axs[0].fill_between(time, (mean[:,2]-std[:,2])/(d**2), (mean[:,2]+std[:,2])/(d**2), alpha = 0.5, color = cmap.to_rgba(c))
	linear = readingInRealCSV("Nbody_Sormani_Data/Varying_Ep/Linear_Test_Particle.csv")
	timeL = np.linspace(0, 100, np.shape(linear)[0])
	axs[0].plot(timeL[:-1], linear[1:,2]/(linearEp**2), color = 'firebrick', label = "Linear")


	linear = readingInRealCSV("Nbody_Sormani_Data/Varying_Ep/Linear_Consistent_Particle.csv")
	timeL = np.linspace(0, 125, np.shape(linear)[0])
	axs[1].plot(timeL[:-2], 1.25*linear[2:,2]/(linearEp**2), color = 'firebrick', label = "Linear")
	epS, epD, epC = ['001', '05', '1'], [0.001, 0.05, 0.1], [0.005, 0.05, 0.1]
	#epS, epD = ['001', '05'], [0.001, 0.05]
	for s, d, c in zip(epS, epD, epC):
		stem = "evolutionC_" + s + "_18"
		mean, std = torqueMeanAndStd(stem, upperIndex = 3, file = "Nbody_Sormani_Data/Varying_Ep/")
		time = np.linspace(0,100, np.shape(mean)[0])
		axs[1].plot(time, mean[:,2]/(d**2), color = cmap.to_rgba(c), label =r"$\log_{10} \epsilon = $ " f"{log(d,10):.2f}")
		axs[1].fill_between(time, (mean[:,2]-std[:,2])/(d**2), (mean[:,2]+std[:,2])/(d**2), alpha = 0.5, color = cmap.to_rgba(c))

	axs[1].set_xlim([0,100])
	
	
	
	xlabelPositions, xlabels = [2*pi*i for i in range(0,16, 3)], [0] + [str(i) + r"$\times2\pi$" for i in range(3,16, 3)]
	axs[0].set_xticks(xlabelPositions)
	axs[0].set_xticklabels(xlabels)
	axs[1].set_xticks(xlabelPositions)
	axs[1].set_xticklabels(xlabels)
	
	axs[0].set_xlabel("Time", fontsize = 15)
	axs[1].set_xlabel("Time", fontsize = 15)
	axs[0].set_ylabel(r"$\tau(t)/\epsilon^{2}$", fontsize = 15)
	axs[0].ticklabel_format(axis='y', style = 'sci', scilimits=(0,0))
	axs[0].legend(fontsize = 12, loc = "lower left")
	axs[0].axhline(0, linestyle = '--', color = 'k')
	axs[1].axhline(0, linestyle = '--', color = 'k')
	#axs[0].set_ylim([-70, 20])
	plt.show()


#sormaniBackwardsFractionOmegaP()
#sormaniBackwardsFractionAmplitude()

#particleEvolution()

#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_18_14.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_18_14.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_18_14.csv"], 0.18, 0.14)
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_1_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_1_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_1_05.csv"], 0.1, 0.05) # This good one for jacobi
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_1_01.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_1_01.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_1_01.csv"], 0.1, 0.01)
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_18_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_18_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_18_05.csv"], 0.1, 0.01)
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_18_005.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_18_005.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_18_005.csv"], 0.1, 0.005)
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning.csv","Nbody_Sormani_Data/Phase_Space_Data/Middle.csv","Nbody_Sormani_Data/Phase_Space_Data/End.csv"], 0.18, 0.005)
varyingEpsilon()


