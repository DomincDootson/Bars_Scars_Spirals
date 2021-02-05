import numpy as np
import matplotlib.pyplot as plt
import csv
from cmath import *

def readingInComplexCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [complex(i) for i in row]
			
			data.append(lst)

		del data[0]
		return np.asarray(data)


def readingInRealCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [float(i) for i in row]
			data.append(lst)

		del data[0]
		return np.asarray(data)


def torqueAndSpeedPlot():
	omega = [r'$0.04t_{0}^{-1}$', r'$0.08t_{0}^{-1}$', r'$0.12t_{0}^{-1}$', r'$0.16t_{0}^{-1}$', r'$0.20t_{0}^{-1}$']
	evolution  = [readingInRealCSV("BarSlowing/evolution4.csv"), readingInRealCSV("BarSlowing/evolution8.csv"), readingInRealCSV("BarSlowing/evolution12.csv"), readingInRealCSV("BarSlowing/evolution16.csv"), readingInRealCSV("BarSlowing/evolution20.csv") ]




	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(2,1, sharex = True)

	time  = np.linspace(0,20, np.shape(evolution[0])[0])
	axs[0].plot(time, evolution[0][:,2], label = omega[0], color = '#191970')
	#plt.plot(time, evolution[1][:,2], label = omega[1], color = '#191970')
	axs[0].plot(time, evolution[2][:,2], label = omega[2], c = '#1034A6')
	#plt.plot(time, evolution[3][:,2], label = omega[3], c = '#1F75FE')
	axs[0].plot(time, evolution[4][:,2], label = omega[4], c = '#1F75FE')
	axs[0].axhline(linestyle = '--')
	axs[0].set_ylabel(r"Torque, $\tau_{b}$")
	axs[0].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))


	axs[1].plot(time, (evolution[0][:,1] - evolution[0][0,1]), label = omega[0], color = '#191970')
	#plt.plot(time, evolution[1][:,2], label = omega[1], color = '#191970')
	axs[1].plot(time, (evolution[2][:,1] - evolution[2][0,1]), label = omega[2], c = '#1034A6')
	#plt.plot(time, evolution[3][:,2], label = omega[3], c = '#1F75FE')
	axs[1].plot(time, (evolution[4][:,1] - evolution[4][0,1]), label = omega[4], c = '#1F75FE')
	axs[1].set_ylabel(r"$\omega-\omega_{0}$")
	axs[1].axhline(linestyle = '--')
	axs[1].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))

	axs[1].set_xlabel(r"$t_{0}$")
	axs[0].legend()
	plt.show()


#torqueAndSpeedPlot()
n = 10

old = readingInComplexCSV("comparison0.csv")
new = readingInComplexCSV("evolution.csv")

plt.plot(old[:,n])
plt.plot(9.5*new[:,n])

plt.show()