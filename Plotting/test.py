from generalPlottingFunctions import *
from Density_Classes.OneDdensity import *
from Density_Classes.TwoDdensity import *

from numpy import linalg as LA

'''leading = TwoDdensity("Spiral_Data/spiralLeading.csv")
trailing = TwoDdensity("Spiral_Data/spiralTrailing.csv")

fig, axs = plt.subplots(nrows = 2, ncols = 2)

axs[0,0].imshow(leading[0])
axs[0,1].imshow(trailing[0])

data = readingInRealCSV("Spiralden1D.csv")


axs[1,0].plot(data[:,0], data[:,1])
axs[1,0].plot(np.linspace(0,10, 101), trailing[0][100:,100])


axs[1,1].plot(data[:,0], data[:,2])
axs[1,1].plot(np.linspace(0,10, 101), leading[0][100:,100])

plt.show()'''
'''
wave = TwoDdensity("Wave.csv")
wave.fourierAnimations(angHarmonic = 2, rMax = 10)
#wave.densityAnimation(rMax = 10, lines = [2.5, 4])
#wave.densityAnimationPolar(rMax = 5)
#plt.imshow(wave[1])
plt.show()'''
# data = OneDdensity("selfconsistent.csv")

# plt.plot(data.radii-8, data.radii*data[50]/(2*pi))

# data = OneDdensity("sheet_consistant.csv")
# plt.plot(data.radii-8, data[50])
# plt.show()
'''data = readingInRealCSV("test.csv")
plt.plot(data[:,0], data[:,1])
plt.plot(data[:,0], data[:,2])
plt.show()'''
'''
CA, CC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentAntiClockwise.csv"), readingInRealCSV("Swing_Data/Amplification_Fixed_5/ConsistentClockwise.csv")
TA, TC = readingInRealCSV("Swing_Data/Amplification_Fixed_5/TestAntiClockwise.csv"), readingInRealCSV("Swing_Data/Amplification_Fixed_5/TestClockwise.csv")

test = readingInRealCSV("swing_test_15.csv")
plt.scatter((4/CA[1:6,0]), CC[1:6,1], label = "Consistent Leading", color = 'firebrick')


#plt.scatter((4/CA[1:6,0]), CA[1:6,1], label = "Consistent Trailing", color = 'royalblue')

#plt.plot((CA[:,0]), TC[:,1], label = "Test Leading", color = 'firebrick', linestyle = '--')
#plt.plot((CA[:,0]), TA[:,1], label = "Test Trailing", color = 'royalblue', linestyle = '--')
plt.plot(4/test[:,0], test[:,1], color = 'firebrick')
test = readingInRealCSV("swing_test_12.csv")
plt.plot(4/test[:,0], test[:,1], color = 'royalblue')
#plt.plot(4/test[:,0], test[:,2])



plt.legend()
plt.ylabel("Max Density")
plt.xlabel(r"$\lambda/\lambda_{crit}$")

plt.yscale('log')
z
plt.show()'''
linear = readingInComplexCSV("KalnajsTorque/Sormani_Bar/coeffTest35.csv")
n = 0
fig,axs = plt.subplots(ncols = 2)
axs[0].plot(np.linspace(0, 100, np.shape(linear[1:100,n])[0]), np.absolute(linear[1:100,n]))
axs[1].plot(np.linspace(0, 100, np.shape(linear[1:100,n])[0]), np.angle(linear[1:100,n]), label = "Linear")

# linear_old = readingInComplexCSV("test.csv")



# axs[0].plot(np.linspace(0, 100, np.shape(linear_old[:100,n])[0]), np.absolute(linear_old[:100,n]))
# axs[1].plot(np.linspace(0, 100, np.shape(linear_old[:100,n])[0]), np.angle(linear_old[:100,n]), label = "Linear")


nbody = readingInComplexCSV("KalnajsTorque/Sormani_Bar/Coefficent_0.csv")

axs[0].plot(np.linspace(0, 100, np.shape(nbody[:100,n])[0]), np.absolute(nbody[:100,n]))
axs[1].plot(np.linspace(0, 100, np.shape(nbody[:100,n])[0]), np.angle(nbody[:100,n]), label = "Test Particle")
axs[0].set_ylabel("Absolute")
axs[1].set_ylabel("Phase")
axs[1].legend()
fig.suptitle(f"n = {n}")

# Coefficent_0.csv





plt.show()