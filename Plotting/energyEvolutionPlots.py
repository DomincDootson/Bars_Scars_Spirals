from generalPlottingFunctions import *

test = EnergyEvolutionData("KalnajsEnergyEvolutionTest.csv")
selfConsistent = EnergyEvolutionData("KalnajsEnergyEvolution.csv")



def energyEvolutionGivenRadius(energyEvolution, energyEvolutionTest, radius):
	littleSigma = energyEvolution.little_sigma()
	angHarmonic = energyEvolution.ang_harmonic()
	radii = energyEvolution.radii()
	print(radii)
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(np.shape(littleSigma)[0], np.shape(angHarmonic)[0], sharex = True)
	time = np.linspace(0,50, energyEvolution.numb_Time_Steps())

	for i in range(np.shape(littleSigma)[0]):
		for j in range(np.shape(angHarmonic)[0]):
			
			axs[i,j].plot(time, energyEvolution.energy_evolution(littleSigma[i], angHarmonic[j], radius), label = 'Self Consistent', color = 'firebrick')
			axs[i,j].plot(time, energyEvolutionTest.energy_evolution(littleSigma[i], angHarmonic[j], radius), label = 'Test Particle', color = 'cornflowerblue')


	for l in range(np.shape(angHarmonic)[0]):
		axs[0,l].set_title(r"$\ell_{p}=$ " + str(angHarmonic[l]))

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r}=$ " + str(littleSigma[s]) + r"$v_{c}$")

	
	axs[0,0].legend()
	axs[-1, 1].set_xlabel(r"Time $[t_{0}]$")
	fig.suptitle(r"Energy Evolution $r_{n} = $ " + str(radius))

	plt.show()

## A function that varies the point where you ping it


def maxEnergyDifferentRadii(energyEvolution, energyEvolutionTest):
	littleSigma = energyEvolution.little_sigma()
	angHarmonic = energyEvolution.ang_harmonic()
	radii = energyEvolution.radii()
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(np.shape(littleSigma)[0], np.shape(angHarmonic)[0], sharex = True)
	

	for i in range(np.shape(littleSigma)[0]):
		for j in range(np.shape(angHarmonic)[0]):
			data = energyEvolution.max_Energy(littleSigma[i], angHarmonic[j])
			axs[i,j].plot(data[:,0], data[:,1], label = 'Self Consistent', color = 'firebrick')

			data = energyEvolutionTest.max_Energy(littleSigma[i], angHarmonic[j])
			axs[i,j].plot(data[:,0], data[:,1], label = 'Test Particle', color = 'cornflowerblue')


	for l in range(np.shape(angHarmonic)[0]):
		axs[0,l].set_title(r"$\ell_{p}=$ " + str(angHarmonic[l]))

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r}=$ " + str(littleSigma[s]) + r"$v_{c}$")

	
	axs[0,0].legend()
	axs[-1, 1].set_xlabel(r"Nudge Radius $[r_{0}]$")
	fig.suptitle("Maximum Energy")

	plt.show()

		

def timeMaxEnergyDifferentRadii(energyEvolution, energyEvolutionTest):	
	littleSigma = energyEvolution.little_sigma()
	angHarmonic = energyEvolution.ang_harmonic()
	radii = energyEvolution.radii()
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(np.shape(littleSigma)[0], np.shape(angHarmonic)[0], sharex = True)
	

	for i in range(np.shape(littleSigma)[0]):
		for j in range(np.shape(angHarmonic)[0]):
			data = energyEvolution.time_Max_Energy_All_Radii(littleSigma[i], angHarmonic[j])
			axs[i,j].plot(energyEvolution.radii(), data, label = 'Self Consistent', color = 'firebrick')

			data = energyEvolutionTest.time_Max_Energy_All_Radii(littleSigma[i], angHarmonic[j])
			axs[i,j].plot(energyEvolutionTest.radii(), data, label = 'Test Particle', color = 'cornflowerblue')


	for l in range(np.shape(angHarmonic)[0]):
		axs[0,l].set_title(r"$\ell_{p}=$ " + str(int(angHarmonic[l])))

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r}=$ " + str(littleSigma[s]) + r"$v_{c}$")

	
	axs[-1,-1].legend()
	axs[-1, 1].set_xlabel(r"Nudge Radius $[r_{0}]$")
	fig.suptitle("Time of Maximum Energy")

	plt.show()



energyEvolutionGivenRadius(selfConsistent, test,7)
