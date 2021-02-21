import numpy as np
import matplotlib.pyplot as plt
import csv
from cmath import *


def readingInRealCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []

		for row in csv_reader:
			lst = [float(i) for i in row]
			data.append(lst)

		return np.asarray(data)