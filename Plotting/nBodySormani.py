from generalPlottingFunctions import *

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
	data = [readingInRealCSV(file) for file in files]
	fig, axs = plt.subplots(ncols = len(data), nrows = 2, sharey = 'row', sharex = 'row')

	# xlim = [np.amin(data[0][:,1]), np.amax(data[0][:,1])]
	# ylim = [np.amin(data[0][:,0]), np.amax(data[0][:,0])]
	ylim = [0, 10]
	xlim = [0, 5]
	jacobi =  1.4
	for ax, d in zip(axs[0,:], data):
		
		ax.hist2d(d[:,0],d[:,1], bins = [np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/100), np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/100)])
		ax.set_xlim([0, 2.45])
		ax.set_ylim([0, 4.5])
		ax.set_xlabel("Energy")
		ax.plot([0.1*L + jacobi for L in np.linspace(0,4.5)], np.linspace(0,4.5))
	
	ylim = [0, 5]
	xlim = [0, 3]
	for ax, d in zip(axs[1,:], data):
		
		ax.hist2d(d[:,2],d[:,1], bins = [np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/100), np.arange(ylim[0], ylim[1], (ylim[1]-ylim[0])/100)])
		ax.axvline(jacobi)
		ax.set_xlim([0, 2.9])
		ax.set_ylim([0, 4.5])
		ax.set_xlabel("Jacobi")

	axs[0,0].set_ylabel("Angular Momentum")
	axs[1,0].set_ylabel("Angular Momentum")

	plt.show()


#sormaniBackwardsFractionOmegaP()
#sormaniBackwardsFractionAmplitude()
#particleEvolution()
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_18_14.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_18_14.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_18_14.csv"], 0.18, 0.14)
changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_1_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_1_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_1_05.csv"], 0.1, 0.05) # This good one for jacobi
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_1_01.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_1_01.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_1_01.csv"], 0.1, 0.01)
#changePhaseSpace(["Nbody_Sormani_Data/Phase_Space_Data/Beginning_18_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/Middle_18_05.csv", "Nbody_Sormani_Data/Phase_Space_Data/End_18_05.csv"], 0.1, 0.01)
