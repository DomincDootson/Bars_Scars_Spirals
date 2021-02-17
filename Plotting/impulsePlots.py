import numpy as np
import matplotlib.pyplot as plt
import csv
from cmath import *



def impulseDensityName(littleSigma, radius, angHarmonic):
	return "../Disk_Kicking/littleSigma_" + str(round(littleSigma*100)) + "/Density" + str(round(radius*10)) +"_" + str(angHarmonic) + ".out" 


def readingInRealCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [float(i) for i in row]
			data.append(lst)

		return np.asarray(data)


def densityEvolutionPlot(littleSigma, timestep = 0.25):
	
	radii = [.1, 1, 5]
	angHarmonic = [0, 1, 2]
	times = [50, 100, 150]



	data = []

	for time in times:
		holdingAng = []
		for l in angHarmonic:
			holdingRadius = []
			for r in radii:
				holdingRadius.append(readingInRealCSV(impulseDensityName(littleSigma, r, l)))
			holdingAng.append(holdingRadius)
		data.append(holdingAng)

	fig, axs = plt.subplots(len(times),len(angHarmonic), sharex = True)

	for i in range(len(times)):
		for j in range(len(angHarmonic)):
			for r in range(len(radii)):
				axs[i,j].plot(data[i][j][r][0,:], data[i][j][r][times[i],:], label = str(radii[r])+r"$r_{0}$")

	

	for l in angHarmonic:
		axs[0,l].set_title(r"$\ell_{p} = $"+str(l))

	for t in range(len(times)):
		axs[t, 0].set_ylabel(str(times[t]*timestep)+r"$t_{0}$")
	

	axs[-1,1].set_xlabel(r"Radius")
	
	axs[0,0].legend()
	fig.suptitle(r"$\sigma=$"+str(littleSigma))



	plt.show()

	



## Could we do a really nice animation?!?!?!?!?!?!?!?!?!?!?!?!?!?!

densityEvolutionPlot(.35)