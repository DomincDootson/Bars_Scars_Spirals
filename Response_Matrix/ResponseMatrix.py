import numpy as np
import csv
from cmath import *

from numpy import linalg as LA


def readingInComplexCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [complex(i) for i in row]
			#lst = [print(row) for i in row] 
			
			data.append(lst)

		return np.asarray(data, dtype=np.cdouble)

class ResponseMatrix():
	def __init__(self, filename):
		self.rm = readingInComplexCSV(filename)

	def getEigenVectors(self):
		self.eigenValues, self.eigenVectors = LA.eig(self.rm)

	def findEigenMode(self):
		self.getEigenVectors()

		closeest = abs(self.eigenValues[0])
		self.eigenMode = self.eigenVectors[:, 0]

		for i, eigenValue in enumerate(self.eigenValues): # Can we check that the eigen mode is always the largest? 
			if ((abs(eigenValue - 1)) < abs(closeest - 1)):
				closeest = eigenValue
				self.eigenMode = self.eigenVectors[:, i]

		print("Closest Eigenvalue ", closeest, abs(closeest))
				


	def saveEigenMode(self, filename):
		self.findEigenMode()

		f = open(filename, "w")
		f.write(str(np.shape(self.rm)[0]-1) + '\n')
		for coeff in self.eigenMode:
			f.write(f"{coeff.real} {coeff.imag} ")

		f.close() 
				





		
for n in ['12', '14', '16', '18', '20']:
	rm = ResponseMatrix(f"RM/Single_Scarred_RM/RM_{n}_25_-95.csv")
	rm.findEigenMode()
	
	rm.saveEigenMode(f"RM/Single_Scarred_RM/EV_{n}_25_-95.out")
		