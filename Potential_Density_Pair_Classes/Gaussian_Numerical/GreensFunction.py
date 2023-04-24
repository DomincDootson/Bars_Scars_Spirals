import matplotlib.pyplot as plt 
import numpy as np
import time 
from math import *

from scipy.special import kn, iv 


class GreenFunction():
	def __init__(self, ell=2, ymin=0.01, ymax=50, vstep=2000):
		self.ell = ell 
		
		self.v = np.linspace(np.cbrt(ymin-1) +1, np.cbrt(ymax-1) +1,vstep)
		self.y = (self.v-1)**3 +1
		self.integral = np.fromiter(map(self.greens_Integral, [1/y for y in self.y]), dtype = float)
		

		self.v_jacobi = 2 * np.square(self.v-1)
	
	## Bessel Integral ##
	## --------------- ##

	def integrand(self, x, r_over_rp):
		if (x>50) and (r_over_rp*x >50):
			return np.exp(-abs(r_over_rp-1)*x)/(2*x*sqrt(r_over_rp))
		else:
			return iv(self.ell, x * r_over_rp)*kn(self.ell, x) if (r_over_rp <=1) else iv(self.ell, x)*kn(self.ell, x*r_over_rp)

	def greens_Integral(self, r_over_rp, nstep = 2000, upper = 20): 
		## This calculates the integral using linear spacing in u, so that x = u*u 
		## This integral doesn't converge for r_over_rp = 1
		if r_over_rp == 1: # So we don't get an error
			return 0
		upper_x = sqrt(upper/abs(r_over_rp-1)) 
		u = np.linspace(0.01, upper_x, nstep)
		x_lst = np.square(u)


		integral = list(map(self.integrand, x_lst, (r_over_rp for _ in x_lst)))
		integral[0] *= 0.5
		integral[-1] *= 0.5

		total = 0
		for igrand, jaco in zip(integral, u):
			total += igrand * 2 * jaco

		return total * (u[1]-u[0])


	## Calculating density ## 
	## ------------------- ##

	def potential_R(self, density_func, R):
		if R ==0:
			return 0 
		density = density_func(self.y * R)

		integrand = (self.y[1]-self.y[0]) * density * self.v_jacobi * self.integral
		return -4*R*np.sum(integrand)

	def potential(self, density_func, R_array):
		potential_with_func = lambda R : self.potential_R(density_func, R)
		return np.fromiter(map(potential_with_func, (R for R in R_array)), dtype = float)
	

	## Include a diganositc plot in the Guassian Container, somthing like ## 

	# rho = lambda R: (1/sqrt(2*pi*(1.10)**2))*np.exp(-np.square(R-15)/(0.5*(1.10**2)))
	# gf = GreenFunction(vstep = 2000)
	# r = np.linspace(0.0, 20, 500)
	# s = time.time()
	# pot = gf.potential(rho, r)
	# e = time.time()
	# print(e-s)
	# plt.plot(r, pot)
	# plt.show()


# rho = lambda R: (1/sqrt(2*pi*(0.0367)**2))*np.exp(-np.square(R-0.5)/(0.5*(0.0367**2)))
# gf = GreenFunction(vstep = 2000)
# r = np.linspace(0.0, 20, 1001)
# s = time.time()
# pot = gf.potential(rho, r)
# e = time.time()
# print(e-s)
# plt.plot(r, pot)
# plt.show()

