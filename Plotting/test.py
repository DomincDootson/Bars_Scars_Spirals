from generalPlottingFunctions import *




#new = readingInComplexCSV("kalnajsTest.csv")
'''

n = 0
old = readingInComplexCSV("evolution200.csv")
plt.plot(np.linspace(0,200, 200), np.absolute(old[:,n]), label = 'Linear calculation')
new = readingInComplexCSV("other.csv")
plt.plot(np.absolute(new[:,n]), label = 'New n body')

new = readingInComplexCSV("kalnajsTest.csv")
#plt.plot(np.absolute(9.7*new[:,n]), label = 'Old n body')
plt.legend()
plt.show()
'''

# 2  + 1.57
'''
n = 8

time = np.linspace(0,20,200)

data = readingInComplexCSV("../evolution200.csv")
#plt.plot(time, np.absolute(data[:,n]), label = 'Linear Perturbation')

data = readingInComplexCSV("../barCoeff.csv")
plt.plot(time, np.absolute(data[:,n]), label = 'Linear Bar')

data = 1.02*readingInComplexCSV("../nBody/coeffEvolutionN.csv")
plt.plot(time, np.absolute(-data[:,n]), label = 'N-body bar')
data = 1.02*readingInComplexCSV("../nBody/coeffEvolutionPerturbationN.csv")
#plt.plot(time, np.absolute(data[:,n]), label = 'N-body Perturbation')
'''

'''
data = readingInRealCSV("../nBody/barCoeffN.csv")
lst = []
for i in range(np.shape(data)[0]):
	lst.append(data[i,0] + 1j*data[i,1])
data = np.asarray(lst)

plt.plot(np.absolute(lst))

data = readingInComplexCSV("../barPerturbationFile.csv")
plt.plot(np.absolute(data[:,0]))
'''

data = readingInRealCSV("../barEvolution.csv")
time  = np.linspace(0,20, np.shape(data)[0])
plt.plot(time, data[:,2], color = 'firebrick', label = 'Linear Response')

data = readingInRealCSV("../nBody/barEvolutionN.csv")
plt.plot(time, -1.02*data[:,2], label = 'N-Body', color = 'cornflowerblue')


'''
data = readingInComplexCSV("../evolution200.csv")
plt.plot(time, np.absolute(data[1:,n]), label = 'Linear Perturbation')

data = readingInComplexCSV("../barCoeff.csv")
plt.plot(time, np.absolute(data[1:,n]), label = 'Linear Bar')

#data = 1.025*readingInComplexCSV("../nBody/coeffEvolutionN.csv")
#plt.plot(time, np.absolute(data[:,n]), label = 'N-body bar')
data = 1.025*readingInComplexCSV("../nBody/coeffEvolutionPerturbationN.csv")
plt.plot(time, np.absolute(data[1:,n]), label = 'N-body Perturbation')
'''
plt.title(r"Torque on slowly rotating gaussian Bar ($0.1t_{0}^{-1}$)")
plt.xlabel("Time")
plt.ylabel("Torque")
plt.legend()




plt.show()


