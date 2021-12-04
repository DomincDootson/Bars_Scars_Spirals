from generalPlottingFunctions import *
	#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
'''
n = 8
plt.subplot(131)

data = readingInComplexCSV("Coeff.csv")

plt.plot(1.035*np.real(data)[:,n])

data = -readingInComplexCSV("evolution.csv")

plt.plot(np.real(data)[:,n])



plt.subplot(132)
data = readingInComplexCSV("Coeff.csv")

plt.plot(1.035*np.imag(data)[:,n])

data = -readingInComplexCSV("evolution.csv")

plt.plot(np.imag(data)[:,n])


plt.subplot(133)
data = readingInComplexCSV("Coeff.csv")

plt.plot(1.035*np.absolute(data)[:,n])

data = -readingInComplexCSV("evolution.csv")

plt.plot(np.absolute(data)[:,n])
plt.show()
'''
fig, axs = plt.subplots(nrows=1, ncols=1)

data = readingInRealOUT("Kalnaj_0.out")


time = np.linspace(0,-80, np.shape(data)[0])

axs.plot(time, data[:,0])
axs.plot(time, data[:,1])

axs.set_xlabel(r"$t-t'$")



plt.show()