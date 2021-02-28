import numpy as np
import matplotlib.pyplot as plt
import csv
#from cmath import *
from math import *

def readingInRealCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [float(i) for i in row]
			data.append(lst)

		return np.asarray(data)

def readingInComplexCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [complex(i) for i in row]
			
			data.append(lst)

		return np.asarray(data)


ENERGY_DIR = "../Disk_Kicking/Energy_Evolution/"

class EnergyEvolutionData(): # This holds the data output by the C++ code
	
	def __init__(self, filename, timeStep = 0.25):
		self.m_data = readingInRealCSV(ENERGY_DIR + filename)
		self.m_timeStep = timeStep

	def check_Agreement(self, row, littleSigma, angHarmonic, radius):
		if (self.m_data[row, 0] == littleSigma) & (self.m_data[row, 1] == angHarmonic) & (self.m_data[row,2] == radius):
			return True
		else:
			return False

	def energy_evolution(self, littleSigma, angHarmonic, radius):
		for i in range(np.shape(self.m_data)[0]):
			if self.check_Agreement(i, littleSigma, angHarmonic, radius):
				return self.m_data[i, 3:]
		
		print("Incorrect values.")
		print("Values entered: ", littleSigma, angHarmonic, radius)
		exit(0)


	def little_sigma(self):
		return np.unique(self.m_data[:,0])

	def ang_harmonic(self):
		return np.unique(self.m_data[:,1])

	def radii(self):
		return np.unique(self.m_data[:,2])


	def max_Energy(self, littleSigma, angHarmonic):
		radii = self.radii()
		lst = []
		for radius in radii:
			lst.append([radius, np.amin(self.energy_evolution(littleSigma, angHarmonic, radius))])
		return np.asarray(lst)

	def numb_Time_Steps(self):
		return np.shape(self.m_data[0,3:])[0]

	def index_Max_Energy(self, littleSigma, angHarmonic, radius):
		maxValue = np.amin(self.energy_evolution(littleSigma, angHarmonic, radius))
		energy = self.energy_evolution(littleSigma, angHarmonic, radius)
		for i in range(self.numb_Time_Steps()):
			if energy[i] == maxValue:
				return i

	def time_Max_Energy(self, littleSigma, angHarmonic, radius):
		return self.m_timeStep * self.index_Max_Energy(littleSigma, angHarmonic, radius)

	def time_Max_Energy_All_Radii(self, littleSigma, angHarmonic):
		return [self.time_Max_Energy(littleSigma, angHarmonic, radius) for radius in self.radii()]