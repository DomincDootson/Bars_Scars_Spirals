import csv
from math import *
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt 
import numpy as np


def inhomogenity(r, r_q, v_0 = 220, A = 0.2):
	return (A* exp(2)/(4 * pi * r_q**2))*np.exp(-2*r/r_q)

def fun(r, y, r_q):
	return np.vstack((y[1], -2 * y[1]/r + 6 * y[0]/r**2 +inhomogenity(r, r_q = r_q)))


def bc(ya, yb):
	return np.array([ya[0], yb[0]])

def write_Phi_2_2_file(filename, r_q = 1.5):
	r = np.linspace(0.1, 10,3)
	phi = np.zeros((2, r.size))
	phi[0,0] = 0.01
	r[1], phi[0,1] = 1.25,  -0.12
	
	def fun_with_param(r, y):
		return fun(r, y, r_q)
	
	sol = solve_bvp(fun_with_param, bc, r, phi)

	x_plot = np.linspace(0.10, 10, 400)
	y_plot = sol.sol(x_plot)[0]

	out = np.zeros((x_plot.size-1 ,2))
	out[:,0] = x_plot[1:]
	out[:,1] = y_plot[1:]/abs(np.min(y_plot[1:])) # Normalise such that the largest amplitudue is 1
	#out[:,1] = y_plot[1:]


	with open(filename, 'w', encoding='UTF8', newline='') as f:
		writer = csv.writer(f, delimiter=',')
		#writer.writerow([np.shape(out)[0], out[1,0]-out[0,0]])
		writer.writerows(out)


	
# write_Phi_2_2_file("Bar_Potentials/Sormani_Large.out", r_q = 2)
# write_Phi_2_2_file("Bar_Potentials/Sormani_Medium.out", r_q = 1.5)
# write_Phi_2_2_file("Bar_Potentials/Sormani_Small.out", r_q = 1.0)

write_Phi_2_2_file("../Plotting/Bar_Data/Sormani_Fitting/Sormani_Large.csv", r_q = 2)
write_Phi_2_2_file("../Plotting/Bar_Data/Sormani_Fitting/Sormani_Medium.csv", r_q = 1.5)
write_Phi_2_2_file("../Plotting/Bar_Data/Sormani_Fitting/Sormani_Small.csv", r_q = 1.0)