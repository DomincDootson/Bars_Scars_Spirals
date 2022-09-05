from generalPlottingFunctions import *
from Density_Classes.OneDdensity import *

from numpy import linalg as LA

'''
data = readingInRealCSV("potential.csv")
r = np.linspace(0,20, np.shape(data[200,210:])[0])
plt.plot(r, 2*3.14/data[200,210:])
plt.plot(r, 2*3.14*r)'''
#plt.plot(r[1:], -r[1:])
'''
fig, axs = plt.subplots(ncols = 1)

data = readingInRealCSV("withScar.csv")
l = np.linspace(0, 3, np.shape(data)[0])

axs.plot(l,data, label = "Scarred Mestel")
axs.set_title("Grooved Mestel")
axs.set_xlabel(r"$J_{\phi}$")
axs.set_ylabel(r"$F_{m}(J_{\phi}, J_{r} = 0)$")

data = readingInRealCSV("withoutScar.csv")
axs.plot(l, data)
'''
'''
n = 3

data = readingInComplexCSV("Spiral_Data/SpiralAnalyticTest.csv")
plt.plot(np.absolute(data[:,n]))

data = readingInComplexCSV("Spiral_Data/SpiralQuasiTest.csv")
plt.plot(np.absolute(data[:,n]))

data = readingInComplexCSV("Spiral_Data/SpiralTrueTest.csv")
plt.plot(np.absolute(data[:,n]))

plt.show()'''
'''
data = readingInRealCSV("Spiral_Data/QuasiKernel.csv")
data1 = readingInRealCSV("Spiral_Data/AnalyticKernel.csv")

plotting = [(data[i,i]-data1[i,i])/data1[i,i] for i in range(1,11)]
print(plotting)

plt.plot(plotting)
plt.show()'''
'''
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, axs = plt.subplots(ncols=3)

data = readingInRealCSV("coeff.csv")

axs[0].plot(data[0:10,0], label = "Real N Body")
axs[0].plot(data[0:10,1], label = "Imag N Body")
axs[0].plot(data[0:10,2], label = "Real Spiral Coefficent")
axs[0].plot(data[0:10,3], label = "Imag Spiral Coefficent")

axs[0].set_xlabel(r"$n$")
axs[0].set_ylabel(r"$B_{n}$")
axs[0].legend()

#plt.axhline()

axs[1].plot(data[:,0])
axs[1].plot(data[:,1])
axs[1].plot(data[:,2])
axs[1].plot(data[:,3])
axs[1].set_ylabel(r"$B_{n}$")
axs[1].set_xlabel(r"$n$")

deltaPower = (np.sum(np.square(data[:,0]) + np.square(data[:,1]))- np.sum(np.square(data[:,2]) + np.square(data[:,3])))/np.sum(np.square(data[:,2]) + np.square(data[:,3]))
print(deltaPower)
n, delta = [20000, 50000, 100000, 500000], [0.0636, 0.0269, 0.0163, deltaPower]

axs[2].plot(n,delta)
axs[2].ticklabel_format(axis = 'x', style = 'sci', scilimits=(0,0))
axs[2].set_xlabel("Number of Particles")
axs[2].set_ylabel(r"Fractional Error")


fig.suptitle("Initial Coefficents for Spiral Sampling")

plt.show()'''
'''n = 0
fig, axs = plt.subplots(ncols=2)
data = readingInComplexCSV("SpiralCoef.csv")
time = np.linspace(0, 50, np.shape(data)[0])

lst = []
for i in range(np.shape(data)[0]):
	total =0
	for j in range(np.shape(data)[1]):
		total += data[i,j].real**2 + data[i,j].imag**2

	lst.append(total)
axs[0].plot(time, np.real(data[:,n]))
axs[1].plot(time, np.imag(data[:,n]))


data = readingInComplexCSV("Spiral_Data/nBodyEvolution.csv")
time = np.linspace(0, 50, np.shape(data)[0])

lst = []
for i in range(np.shape(data)[0]):
	total =0
	for j in range(np.shape(data)[1]):
		total += data[i,j].real**2 + data[i,j].imag**2

	lst.append(total)


axs[0].plot(time, np.real(data[:,n]))
axs[1].plot(time, np.imag(data[:,n]))

plt.show() '''



'''
oneDDensity = OneDdensity("Waves_Data/RadiusPull_0_5.csv")
oneDDensity.densityAnimation(write2file = "Waves_Plots/Waves_Videos/RadiusPull_0_5.mp4",remove_ic = True)

oneDDensity = OneDdensity("Waves_Data/RadiusPull_0_10.csv")
oneDDensity.densityAnimation(write2file = "Waves_Plots/Waves_Videos/RadiusPull_0_10.mp4",remove_ic = True)

oneDDensity = OneDdensity("Waves_Data/RadiusPull_0_13.csv")
oneDDensity.densityAnimation(write2file = "Waves_Plots/Waves_Videos/RadiusPull_0_13.mp4",remove_ic = True)
'''
'''
m1 = [1, 2, 3, 4, 5, 6, 7]

real = [0.0395109, 0.00481113, 0.000799331, 0.000185919, 0.000116288, 9.12043e-05, 7.28621e-05]

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(m1, real, color = 'firebrick', label = "Real")

plt.ylabel(r"$\left|\hat{\mathcal{M}}(\omega)-I\right|$")
plt.xlabel(r"Max $m_{1}$")
plt.yscale('log')
plt.title(r"Mestel Mode $\omega = 0.88 +i0.13$")
plt.show()




plt.plot(m1, real)
'''

data = readingInRealCSV("../wavetest.csv")
plt.plot(data[:,0], data[:,1])
plt.plot(data[:,0], data[:,2])
plt.show()

