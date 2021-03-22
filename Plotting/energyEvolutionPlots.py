from generalPlottingFunctions import *

#test = EnergyEvolutionData("KalnajsEnergyEvolutionTest.csv")
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
			axs[i,j].plot(data[:,0], data[:,1]/norm, label = 'Self Consistent', color = 'firebrick')

			data = energyEvolutionTest.max_Energy(littleSigma[i], angHarmonic[j])
			norm = np.amax(np.absolute(data[:25,1]))
			axs[i,j].plot(data[:,0], data[:,1]/norm, label = 'Test Particle', color = 'cornflowerblue')


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



def taperExtension(rInner, rOuter, test = False):
	
	if rInner >= 1:
		filename = "KalnajsEnergyEvolution_"+str(int(10*rInner))+"_"+str(rOuter) 
	else:
		order = abs(floor(log(rInner,10)))
		filename = "KalnajsEnergyEvolution_" + order*"0"+str(int((10**order)*rInner))+"_"+str(rOuter)

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


rInner,rOuter = [.1, .5, 1.3, 1, 2], [10, 10, 10, 10, 10]
rInner,rOuter = [.1, .5, 1, 1.3, 2], [10, 10, 10, 10, 10]
maxEnergyDifferentTapers(rInner, rOuter)
#rInner, rOuter = [1], [10]
#maxEnergyDifferentTapersl1Harmonic(rInner, rOuter)
# #timeMaxEnergyDifferentTapers(rInner, rOuter)

selfConsistent = EnergyEvolutionData("KalnajsEnergyEvolution_10_15.csv", 1, 15) 
test = EnergyEvolutionData("KalnajsEnergyEvolution_10_15.csv", 1, 15)


maxEnergyDifferentRadii(selfConsistent,test)


# energyEvolutionGivenRadius(selfConsistent, test,1)
# energyEvolutionGivenRadius(selfConsistent, test,2)
# energyEvolutionGivenRadius(selfConsistent, test,3)
# energyEvolutionGivenRadius(selfConsistent, test,5)

# maxEnergyDifferentRadii(selfConsistent, test)

# timeMaxEnergyDifferentRadii(selfConsistent, test)