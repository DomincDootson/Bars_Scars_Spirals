import csv
from math import *
import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import solve_bvp

from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm

from GreensFunction import * 


class Gaussian_BF():
	"""docstring for Gaussian_BF"""
	def __init__(self, r_0, sigma, l):
		self.r_0, self.sigma = r_0, sigma
		self.l = l


	def density(self, r):
		return (1/sqrt(2*pi*(self.sigma)**2))*np.exp(-0.5*np.square(r-self.r_0)/(self.sigma)**2)


	# Then some potentenal solver 

	def obtain_potential(self, potential):
		self.potential = potential

class Guassian_Container():
	
	def __init__(self, r_I = 0.5, r_O = 15, n = 49, l = 2):
		
		self.r_0s = [10**log_r for log_r in np.linspace(log(r_I,10), log(r_O,10), num = n)]
		self.l = l 
		
		self.gf = GreenFunction(ell = self.l)

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

		

	def calculate_potential(self, r_array):
		for bf in self:
			potential = self.gf.potential(bf.density, r_array) 
			bf.obtain_potential(potential)

	## Some Plotting Functions ## 
	## ----------------------- ##

	def plot_bf(self):
		plt.rc('text', usetex=True)
		plt.rc('font', family='serif')
		fig, axs = plt.subplots(nrows =2, sharex = True)
		r = np.linspace(0,20, 1001)
		cmap_r = ScalarMappable(cmap = 'plasma_r', norm = LogNorm(vmin=0.5, vmax=15))
		cmap_s = ScalarMappable(cmap = 'plasma_r', norm = LogNorm(vmin=self.sigmas[0], vmax=self.sigmas[-1]))

		self.calculate_potential(r)
		for bf in self[3::5]:
			axs[0].plot(r, bf.density(r)*bf.sigma, color = cmap_r.to_rgba(bf.r_0))
			axs[1].plot(r, bf.potential, color = cmap_s.to_rgba(bf.sigma))

		
		cbar = fig.colorbar(cmap_r, ax=axs[0]).set_label(label = r"$R_{n}$", rotation = 270, fontsize = 15) 
		cbar = fig.colorbar(cmap_s, ax=axs[1]).set_label(label = r"$\sigma_{n}$", rotation = 270, fontsize = 15) 


		axs[1].set_xlabel("Radius", fontsize = 15)
		axs[0].set_ylabel(r"$\rho_{n}(R)\sigma_{n}$", fontsize = 15)
		axs[1].set_ylabel(r"$\psi_{n}(R)$", fontsize = 15)
		axs[0].set_title(r"Gaussian Potential Density Pairs", fontsize = 15)

		plt.show()

	## Saving Functions ##
	## ---------------- ##


	def write_Gaussian_2_file(self, filename, r_inner = 0, r_outer = 20):
		r = np.linspace(r_inner, r_outer, num = 1001)

		with open(filename, 'w', encoding='UTF8', newline='') as f:
			writer = csv.writer(f, delimiter=' ')
			writer.writerow([len(self)-1, self[0].l, (1)/(len(r)-1), r[-1]])

			self.calculate_potential(r)
			for bf in self:
				#pot = bf.obtain_potential(r)
				
				writer.writerow(bf.potential)
				
				den = [(1/sqrt(2*pi*(bf.sigma)**2))*(exp(-0.5*((r_v-bf.r_0)/1)**2)) for r_v in r]
				writer.writerow(np.asarray(den))
				#writer.writerow(np.asarray(bf.density(r)))

				#plt.plot(r, bf.density(r))
		

		print(f"Written Gaussian BF to: {filename}")
		#plt.show()

gc = Guassian_Container()

gc.plot_bf()



#gc.write_Gaussian_2_file('Gaussian_Numerical.dat')

# r = np.linspace(0,15,1001)
# plt.plot(r, gc[0].density(r))
# plt.show()
