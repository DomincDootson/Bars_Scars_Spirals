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


def plottingVaringOmega():
	omega = [r'$0.04t_{0}^{-1}$', r'$0.08t_{0}^{-1}$', r'$0.12t_{0}^{-1}$', r'$0.16t_{0}^{-1}$', r'$0.20t_{0}^{-1}$']
	evolution  = [readingInRealCSV("BarEvolution/VaryingOmega/evolution4.csv"), readingInRealCSV("BarEvolution/VaryingOmega/evolution8.csv"), 
	readingInRealCSV("BarEvolution/VaryingOmega/evolution12.csv"), readingInRealCSV("BarEvolution/VaryingOmega/evolution16.csv"), readingInRealCSV("BarEvolution/VaryingOmega/evolution20.csv") ]




	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots(2,1, sharex = True)

	time  = np.linspace(0,50, np.shape(evolution[0])[0])
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



def plottingVaryingKka():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	Kka = [4, 5, 6, 7]
	Rmax = []
	labels =[r"$0.35 r_{Ka}$", r"$0.32 r_{Ka}$", r"$0.29 r_{Ka}$", r"$0.27 r_{Ka}$"] 

	evolution = [readingInRealCSV("BarEvolution/VaryingKka/evolution4.csv"), readingInRealCSV("BarEvolution/VaryingKka/evolution5.csv"), readingInRealCSV("BarEvolution/VaryingKka/evolution6.csv"), readingInRealCSV("BarEvolution/VaryingKka/evolution7.csv")]

	fig, axs = plt.subplots(2,1, sharex = True)

	time  = np.linspace(0,50, np.shape(evolution[0])[0])
	axs[0].plot(time, evolution[0][:,2], label = labels[0], color = '#191970')
	axs[0].plot(time, evolution[1][:,2], label = labels[1], c = '#1034A6')
	axs[0].plot(time, evolution[3][:,2], label = labels[3], c = '#1F75FE')
	axs[0].axhline(linestyle = '--')
	axs[0].set_ylabel(r"Torque, $\tau_{b}$")
	axs[0].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))


	axs[1].plot(time, (evolution[0][:,1] - evolution[0][0,1]), label = labels[0], color = '#191970')
	axs[1].plot(time, (evolution[1][:,1] - evolution[2][0,1]), label = labels[1], c = '#1034A6')
	axs[1].plot(time, (evolution[3][:,1] - evolution[3][0,1]), label = labels[3], c = '#1F75FE')
	axs[1].set_ylabel(r"$\omega-\omega_{0}$")
	axs[1].axhline(linestyle = '--')
	axs[1].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))

	axs[1].set_xlabel(r"$t_{0}$")
	axs[0].legend()
	plt.show()

def plottingVaryingRka():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	Kka = [5, 10, 15, 20]
	Rmax = []
	labels =[r"$5 r_{0}$", r"$10 r_{0}$", r"$15 r_{0}$", r"$20 r_{0}$"] 

	evolution = [readingInRealCSV("BarEvolution/VaryingRka/evolution5.csv"), readingInRealCSV("BarEvolution/VaryingRka/evolution10.csv"), readingInRealCSV("BarEvolution/VaryingRka/evolution15.csv"), readingInRealCSV("BarEvolution/VaryingRka/evolution20.csv")]

	fig, axs = plt.subplots(2,1, sharex = True)

	time  = np.linspace(0,50, np.shape(evolution[0])[0])
	axs[0].plot(time, evolution[1][:,2], label = labels[1], color = '#191970')
	axs[0].plot(time, evolution[2][:,2], label = labels[2], c = '#1034A6')
	axs[0].plot(time, evolution[3][:,2], label = labels[3], c = '#1F75FE')
	axs[0].axhline(linestyle = '--')
	axs[0].set_ylabel(r"Torque, $\tau_{b}$")
	axs[0].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))


	axs[1].plot(time, (evolution[1][:,1] - evolution[1][0,1]), label = labels[1], color = '#191970')
	axs[1].plot(time, (evolution[2][:,1] - evolution[2][0,1]), label = labels[2], c = '#1034A6')
	axs[1].plot(time, (evolution[3][:,1] - evolution[3][0,1]), label = labels[3], c = '#1F75FE')
	axs[1].set_ylabel(r"$\omega-\omega_{0}$")
	axs[1].axhline(linestyle = '--')
	axs[1].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))

	axs[1].set_xlabel(r"$t_{0}$")
	axs[0].legend()
	plt.show()



plottingVaringOmega()
plottingVaryingKka()
plottingVaryingRka()



