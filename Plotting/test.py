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
data = OneDdensity("selfconsistent.csv")

plt.plot(data.radii-8, data.radii*data[50]/(2*pi))

data = OneDdensity("sheet_consistant.csv")
plt.plot(data.radii-8, data[50])
plt.show()

#data.densityAnimation()