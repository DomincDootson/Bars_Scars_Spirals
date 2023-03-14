import csv
from math import *
import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import solve_bvp

from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm


class Gaussian_BF():
	"""docstring for Gaussian_BF"""
	def __init__(self, r_0, sigma, l):
		self.r_0, self.sigma = r_0, sigma
		self.l = l


	def density(self, r):
		return (1/sqrt(2*pi*(self.sigma)**2))*np.exp(-np.square(r-self.r_0)/(2 *self.sigma)**2)

	# Then some potentenal solver 

	def fun(self, r, y): # assumed l(l+1) = 6 (i.e. l = 2)
		return np.vstack((y[1], -2 * y[1]/r + self.l*(self.l+1) * y[0]/r**2 +4*pi*self.density(r)))


	def bc(self, ya, yb):
		return np.array([ya[0], yb[0]])


	def obtain_potential(self, r_lst):
		r = np.linspace(0.1, 50,3)
		phi = np.zeros((2, r.size))
		phi[0,0] = 0.01
		r[1], phi[0,1] = self.r_0,  self.density(self.r_0)
		
		sol = solve_bvp(self.fun, self.bc, r, phi)
		pot = sol.sol(r_lst)[0]

		self.potential = [min(0, p) for p in pot]
		
		return self.potential
		

	def laplacian_check(self, r_lst): # Code to test that we are solving Poissions equation properly, we are. 
		plt.plot(r_lst, 4*pi* self.density(r_lst))

		pot = self.obtain_potential(r_lst)

		delta = abs(r_lst[1]-r_lst[0])
		r = [r_lst[i] for i in range(len(pot[1:-1]))]
		phi_0 = [pot[i]/(r_lst[i]**2) for i in range(len(pot[1:-1]))]
		phi_1 = [(2/(2*delta*r_lst[i]))*(pot[i+1] - pot[i-1]) for i in range(len(pot[1:-1]))]
		phi_2 = [(1/((delta**2)))*(pot[i+1] + pot[i-1] - 2 * pot[i]) for i in range(len(pot[1:-1]))]

		den = np.asarray(phi_2) + np.asarray(phi_1) - self.l*(self.l+1) * np.asarray(phi_0)
		plt.plot(r, den)

		plt.show()



class Guassian_Container():
	
	def __init__(self, r_I = 0.5, r_O = 15, n = 49, l = 2):
		
		self.r_0s = [10**log_r for log_r in np.linspace(log(r_I,10), log(r_O,10), num = n)]
		
		
		new_upper = log(self.r_0s[-1],10) + log(self.r_0s[-1],10) - log(self.r_0s[-2],10)
		self.sigmas = [10**r_u-r_l for r_u, r_l in zip(np.linspace(log(self.r_0s[1],10), new_upper, num = n), self.r_0s)]

		self.bf = [Gaussian_BF(r, s, l) for r, s in zip(self.r_0s, self.sigmas)]

	def __getitem__(self, index):
		return self.bf[index]
	
	def __iter__(self):
		for bf in self.bf:
			yield bf 
	def __len__(self):
		return len(self.bf)


	def calculate_potential():
		for bf in self:
			bf.obtain_potential()

	## Some Plotting Functions ## 
	## ----------------------- ##

	def plot_bf(self):
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		fig, axs = plt.subplots(nrows =2)
		r = np.linspace(0,20, 1001)
		cmap_r = ScalarMappable(cmap = 'plasma_r', norm = LogNorm(vmin=0.5, vmax=15))
		cmap_s = ScalarMappable(cmap = 'plasma_r', norm = LogNorm(vmin=self.sigmas[0], vmax=self.sigmas[-1]))

		for bf in self[3::5]:
			axs[0].plot(r, bf.density(r)*bf.sigma, color = cmap_r.to_rgba(bf.r_0))
			axs[1].plot(r, bf.obtain_potential(r)/bf.sigma, color = cmap_s.to_rgba(bf.sigma))

		
		cbar = fig.colorbar(cmap_r, ax=axs[0]).set_label(label = r"$R_{n}$", fontsize = 15) 
		cbar = fig.colorbar(cmap_s, ax=axs[1]).set_label(label = r"$\sigma_{n}$", fontsize = 15) 


		axs[1].set_xlabel("Radius", fontsize = 15)
		axs[0].set_ylabel(r"$\rho_{n}(R)\sigma_{n}$", fontsize = 15)
		axs[1].set_ylabel(r"$\psi_{n}(R)/\sigma_{n}$", fontsize = 15)
		axs[0].set_title(r"Gaussian Potential Density Pairs", fontsize = 15)

		plt.show()

	## Saving Functions ##
	## ---------------- ##


	def write_Gaussian_2_file(self, filename, r_inner = 0, r_outer = 20):
		r = np.linspace(r_inner, r_outer, num = 1001)

		with open(filename, 'w', encoding='UTF8', newline='') as f:
			writer = csv.writer(f, delimiter=' ')
			writer.writerow([len(self)-1, self[0].l, r[1]-r[0], r[-1]])

			for bf in self:
				pot = bf.obtain_potential(r)
				writer.writerow(pot)
				writer.writerow(bf.density(r))
		




gc = Guassian_Container()
gc.write_Gaussian_2_file('Gaussian_Numerical.dat')
