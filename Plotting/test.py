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

data = readingInRealCSV("fitGaussian.csv")
#data = np.log(np.absolute(data))
plt.plot(data)
data = readingInRealCSV("trueKalnajs.csv")
plt.plot(data)
plt.show()