from generalPlottingFunctions import *

test = EnergyEvolutionData("GaussianLogEnergy_10_100.csv",0,0)
#selfConsistent = EnergyEvolutionData("KalnajsEnergyEvolution.csv")



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
			norm = np.amax(np.absolute(data[:25,1]))
			axs[i,j].plot(data[:,0], data[:,1], label = 'Self Consistent', color = 'firebrick')

			data = energyEvolutionTest.max_Energy(littleSigma[i], angHarmonic[j])
			norm = np.amax(np.absolute(data[:25,1]))
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



def taperExtension(rInner, rOuter, filename = "KalnajsEnergyEvolution", test = False):
	
	if rInner >= 1:
		filename = filename + "_" + str(int(10*rInner))+"_"+str(rOuter) 
	else:
		order = abs(floor(log(rInner,10)))
		filename = filename + "_"   + order*"0"+str(int((10**order)*rInner))+"_"+str(rOuter)

	if (test):
		return filename+"_Test.csv"
	else:
		return filename+".csv"


def maxEnergyDifferentTapers(rInner, rOuter):
	energyClasses = []
	for i in range(len(rInner)):
		energyClasses.append(EnergyEvolutionData(taperExtension(rInner[i], rOuter[i]), rInner[i], rOuter[i]))

	littleSigma = energyClasses[0].little_sigma()
	angHarmonic = energyClasses[0].ang_harmonic()
	radii = energyClasses[0].radii()
	
	fig, axs = plt.subplots(np.shape(littleSigma)[0], np.shape(angHarmonic)[0], sharex = True)
	
	for i in range(np.shape(littleSigma)[0]):
		for j in range(np.shape(angHarmonic)[0]):
			for energyEvolution in energyClasses:
				data = energyEvolution.max_Energy(littleSigma[i], angHarmonic[j])
				norm = np.amax(np.absolute(data[:25,1]))
							
				axs[i,j].plot(data[:,0], data[:,1], label = energyEvolution.m_rInner)
				#axs[i,j].plot(data[:,0], data[:,1], label = energyEvolution.m_rInner)



	for l in range(np.shape(angHarmonic)[0]):
		axs[0,l].set_title(r"$\ell_{p}=$ " + str(int(angHarmonic[l])))

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r}=$ " + str(littleSigma[s]) + r"$v_{c}$")

	
	axs[-1,-1].legend()
	axs[-1, 1].set_xlabel(r"Nudge Radius $[r_{0}]$")
	fig.suptitle("Maximum Energy")

	plt.show()


def maxEnergyDifferentTapersl1Harmonic(rInner, rOuter):
	energyClasses = []
	for i in range(len(rInner)):
		energyClasses.append(EnergyEvolutionData(taperExtension(rInner[i], rOuter[i]), rInner[i], rOuter[i]))

	littleSigma = energyClasses[0].little_sigma()
	angHarmonic = energyClasses[0].ang_harmonic()
	radii = energyClasses[0].radii()
	
	fig, axs = plt.subplots(np.shape(littleSigma)[0], np.shape(angHarmonic)[0], sharex = True)
	
	for i in range(np.shape(littleSigma)[0]):
		for j in range(np.shape(angHarmonic)[0]):
			for energyEvolution in energyClasses:
				data = energyEvolution.max_Energy(littleSigma[i], 1)
				if j ==0:
					norm = 1
				elif j ==2:
					norm = np.amax(np.absolute(data[25:,1]))
				elif j ==1:
					norm = np.amax(np.absolute(data[:25,1]))
							
				axs[i,j].plot(data[:,0], data[:,1]/norm, label = energyEvolution.m_rInner)
				#axs[i,j].plot(data[:,0], data[:,1], label = energyEvolution.m_rInner)



	normtile = ["No Normalisation", "Normalised on 1st Peak", "Normalised on 2nd Peak"]
	for l in range(np.shape(angHarmonic)[0]):
		axs[0,l].set_title(normtile[l])

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r}=$ " + str(littleSigma[s]) + r"$v_{c}$")

	
	axs[-1,-1].legend()
	axs[-1, 1].set_xlabel(r"Nudge Radius $[r_{0}]$")
	fig.suptitle(r"Maximum Energy for $\ell_{p}=1$")

	plt.show()





def timeMaxEnergyDifferentTapers(rInner, rOuter):
	energyClasses = []
	for i in range(len(rInner)):
		energyClasses.append(EnergyEvolutionData(taperExtension(rInner[i], rOuter[i]), rInner[i], rOuter[i]))

	littleSigma = energyClasses[0].little_sigma()
	angHarmonic = energyClasses[0].ang_harmonic()
	radii = energyClasses[0].radii()
	
	fig, axs = plt.subplots(np.shape(littleSigma)[0], np.shape(angHarmonic)[0], sharex = True)
	
	for i in range(np.shape(littleSigma)[0]):
		for j in range(np.shape(angHarmonic)[0]):
			for energyEvolution in energyClasses:
				data = energyEvolution.time_Max_Energy_All_Radii(littleSigma[i], angHarmonic[j])
				axs[i,j].plot(energyEvolution.radii(), data, label = energyEvolution.m_rInner)
				#axs[i,j].plot(data[:,0], data[:,1], label = 'Self Consistent')



	for l in range(np.shape(angHarmonic)[0]):
		axs[0,l].set_title(r"$\ell_{p}=$ " + str(int(angHarmonic[l])))

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r}=$ " + str(littleSigma[s]) + r"$v_{c}$")

	
	axs[-1,-1].legend()
	axs[-1, 1].set_xlabel(r"Nudge Radius $[r_{0}]$")
	fig.suptitle("Time of Maximum Energy")

	plt.show()

def taperLabel(rInner, rOuter):
	return r"$R_{i} = $ " +str(rInner) + r", $R_{o} = $ " + str(rOuter) 

def gaussianVaryingTapers(listOfEnergies):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig, axs = plt.subplots(3,3)
	radii = listOfEnergies[0].radii()
	nSigmas, littleSigma = np.shape(listOfEnergies[0].little_sigma())[0], listOfEnergies[0].little_sigma()
	angHarmonic = listOfEnergies[0].ang_harmonic()

	for s in range(nSigmas):
		for l in listOfEnergies[0].ang_harmonic():
			for i in range(len(listOfEnergies)):
				print(listOfEnergies[i].max_Energy(littleSigma[s], l))
				axs[s, int(l)].plot(radii, listOfEnergies[i].max_Energy(littleSigma[s], l), label = taperLabel(listOfEnergies[i].m_rInner, listOfEnergies[i].m_rOuter))

	
	axs[-1,-1].legend()	
	
	for l in range(np.shape(angHarmonic)[0]):
		axs[0, l].set_title(r"$\ell_{p} = $ " + str(angHarmonic[l]))

	for s in range(np.shape(littleSigma)[0]):
		axs[s, 0].set_ylabel(r"$\sigma_{r} = $ " + str(littleSigma[s]) + r"$v_{c}$")		
	
	plt.show()


def energyFromCoeff(data): #we assume that it is all orthogonal
	energies = []

	for time in range(np.shape(data)[0]):
		sum_ = 0
		for n in range(np.shape(data)[1]):
			sum_ += data[time, n].real**2 + data[time, n].imag**2

		energies.append(sum_)

	return energies


def varyingTaperPositition(files = ["General_Data/kalnajsComparison_10_1.out", "General_Data/kalnajsComparison_20_1.out", "General_Data/kalnajsComparison_10_2.out", "General_Data/kalnajsComparison_20_2.out"]):
	energies = [energyFromCoeff(readingInComplexCSV(file)) for file in files]	

	fig, axs = plt.subplots()
	axs.plot(energies[0], label = r"$N_{max} = 10$ $R_{i}=1$")
	axs.plot(energies[1], label = r"$N_{max} = 20$ $R_{i}=1$")
	axs.plot(energies[2], label = r"$N_{max} = 10$ $R_{i}=2$")
	axs.plot(energies[3], label = r"$N_{max} = 20$ $R_{i}=2$")

	axs.set_yscale('log')
	plt.legend()

	plt.show()
	


'''
listOfEnergies = [EnergyEvolutionData("GaussianLogEnergy_5_150.csv",0.5,15), EnergyEvolutionData("GaussianLogEnergy_10_150.csv",1,15), EnergyEvolutionData("GaussianLogEnergy_20_150.csv",2,15)]
gaussianVaryingTapers(listOfEnergies)
listOfEnergies = [EnergyEvolutionData("GaussianLogEnergy_10_100.csv",1,10), EnergyEvolutionData("GaussianLogEnergy_10_150.csv",1,15), EnergyEvolutionData("GaussianLogEnergy_10_175.csv",1,17.5)]
gaussianVaryingTapers(listOfEnergies)
'''



varyingTaperPositition()


