## This Class needs to be split up into two, one that does indivudal k and one that combines different k 

from math import *
import matplotlib.pyplot as plt
import numpy as np

class ShearingSheet():
	
	def __init__(self, ky, Q, omega = 1, kappa = sqrt(2), xi = 1):	# ky is in units of k_crit
		
		self.omega, self.kappa, self.sigma_0 = omega, omega * kappa, xi /(2*pi)
		self.k_crit = (self.kappa**2)/(2 * pi * self.sigma_0) 

		self.ky, self.Q = ky * self.k_crit, Q
		
		self.oort_A = self.omega * (1 - 0.25 * (self.kappa/self.omega)**2)

		self.time_step = 0.1 * (1/self.kappa)

	## Kernel Functions ## 
	## ---------------- ##

	def S(self, t):
		return sin(self.kappa * t)

	def C(self, t):
		return cos(self.kappa * t)

	def bx(self, t, tp): # Absorb the k/k_crit into def
		return (self.ky/self.k_crit) * (self.oort_A * (tp*self.S(tp) - t * self.S(t)) + (self.omega/self.kappa) * (self.C(tp) - self.C(t)))

	def by(self, t, tp):
		return (self.ky/self.k_crit) * (self.oort_A * (tp*self.C(tp) - t * self.C(t)) - (self.omega/self.kappa) * (self.S(tp) - self.S(t)))

	def mag_b(self, t, tp):
		return self.bx(t, tp)**2 + self.by(t, tp)**2

	def cx(self, tp):
		return -self.oort_A * tp * self.C(tp) + (self.omega/self.kappa) * (self.S(tp))

	def cy(self, tp):
		return  self.oort_A * tp * self.S(tp) + (self.omega/self.kappa) * (self.C(tp))

	def b_dot_c(self, t, tp):
		return self.bx(t, tp) * self.cx(tp) + self.by(t, tp) * self.cy(tp)

	def kernel(self, t, tp):
		return 4 * self.b_dot_c(t, tp) * exp(-0.572 * self.Q**2 * self.mag_b(t, tp)) / sqrt(1 + 4 * self.oort_A**2 * tp**2)

	def kernel_grid(self, t_array):
		kernels = self.kappa*np.fromiter(map(self.kernel, (t_array[-1] for _ in t_array), t_array), dtype = float)
		kernels[0] *= 0.5
		kernels[-1] *= 0.5
		return kernels


	## Delta Function Evolution ##
	## ------------------------ ##

	def delta_evolution(self, ti): 
		t_end = 5
		t_array = np.arange(ti * pi /self.kappa, t_end * pi /self.kappa, self.time_step)
		density_test_ptle, density_consistent = np.zeros_like(t_array), np.zeros_like(t_array)
		
		for t_index in range(1, t_array.shape[0]-1): 
			kernels = self.kernel_grid(t_array[:t_index+1])
			
			pert = kernels[0] * (2/self.kappa)  
			density_int = np.sum(kernels[:t_index] * density_consistent[:t_index]) * self.time_step

			density_consistent[t_index] = (pert + density_int) / (1 - kernels[-1])
			density_test_ptle[t_index]  = pert

		return (self.kappa/pi) * t_array[:-2], density_test_ptle[:-2],  density_consistent[:-2]

	def delta_amplification_ti(self, ti):
	 response  = self.delta_evolution(ti)
	 return np.amax(response[2])/np.amax(response[1])

	def delta_amplification(self):
		ti = np.linspace(-1.5, -0., 20)
		list_max = [sheet.delta_amplification_ti(t) for t in ti]
		return max(list_max)





	## Cloud Evolution ## 
	## --------------- ## 

	def mag_k_sq(self, tp):
		return self.ky * (1 + (2 * self.oort_A * tp)**2)

	def cloud_pert(self, tp, M = 1, delta = 0):
		return  M * exp(-0.5 * self.mag_k_sq(tp) * (delta **2)) 

	def cloud_evolution(self, ti, M = 1, delta = 0): 
		t_end = 5
		t_array = np.arange(ti * pi /self.kappa, t_end * pi /self.kappa, self.time_step)
		density_test_ptle, density_consistent, pert = np.zeros_like(t_array), np.zeros_like(t_array), np.asarray(list(map(self.cloud_pert, t_array, (M for _ in t_array), (delta for _ in t_array))))
		
		for t_index in range(1, t_array.shape[0]-1): 
			kernels = self.kernel_grid(t_array[:t_index+1])
			
			pert_int    = np.sum(kernels * pert[:t_index+1])    * self.time_step
			density_int = np.sum(kernels[:t_index] * density_consistent[:t_index]) * self.time_step

			density_consistent[t_index] = (pert_int + density_int) / (1 - kernels[-1])
			density_test_ptle[t_index]  = pert_int

		return (self.kappa/pi) * t_array[:-2], density_test_ptle[:-2],  density_consistent[:-2]









