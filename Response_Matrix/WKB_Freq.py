from math import *

import numpy as np
import matplotlib.pyplot as plt


# Do we want to define a mestel class?

class Mestel_Disc():

	def __init__(self, sigmaR):

		self.sigmaR = sigmaR

		self.vc = 1
		self.Sigma0 = 1/(2*pi)
		self.epsilon = 0.125 # For the time being we will ignore softening 
		self.m = 2

	## General Functions ##
	## ----------------- ##

	def Omega(self, radius):
		return self.vc/radius

	def kappa(self, radius):
		return self.Omega(radius) * sqrt(2)

	def Sigma(self, radius, activeFraction = 1/2):
		return activeFraction * self.Sigma0 * (1/radius)

	def k_crit(self, radius):
		return self.kappa(radius)**2 /(2*pi*self.Sigma(radius))

	def CR(self, omega, m = 2):
		return m / omega

	def ILR(self, omega, m = 2):
		return self.CR(omega, m) * (2-sqrt(2))*0.5

	def OLR(self, omega, m =2):
		return self.CR(omega, m) * (2+sqrt(2))*0.5
	

	## Reduction Factor ##
	## ---------------- ##

	def reduction_factor(self, s, chi, stepSize = 0.1):
		if chi ==0:
			return 1
		if s==1:
			s=0.9999999 # The two lines tend to each other 
    
		tau = np.arange(0, pi, stepSize)
	    
		integrand = np.exp(-chi *(1 + np.cos(tau)))*np.sin(s * tau)*np.sin(tau)
		return stepSize*((1-s**2)/sin(pi * s)) * np.sum(integrand)

	def partialRFpartialS(self, s, chi, deltaS = 0.0001):
		return (1/(2*deltaS)) * (self.reduction_factor(s+deltaS, chi) - self.reduction_factor(s-deltaS, chi))

	def partialRFpartialChi(self, s, chi, deltaChi = 0.0001):
		return (1/(2*deltaChi)) * (self.reduction_factor(s, chi+deltaChi) - self.reduction_factor(s, chi-deltaChi))

	def chi(self, k, radius):
		return ((self.sigmaR*k)/self.kappa(radius))**2

	def s(self, omega0, radius):
		return (omega0 - self.m *self.Omega(radius))/self.kappa(radius)

	def omegaFromS(self, s, radius):
		return s * self.kappa(radius) + self.m * self.Omega(radius)

	## Dispersion relation ## 
	## ------------------- ##
		

	def LHS_disp(self, s, k, radius):
		chi = self.chi(k, radius)
		return s**2 + 2 * pi * self.Sigma(radius)/(self.kappa(radius)**2) * self.reduction_factor(s, chi)*abs(k) * exp(-self.epsilon * abs(k))-1

	def s_from_k(self, k, radius, nstep = 100): # returns s
	    upperVal, middleVal, lowerVal = .99999, 0.5, 0.000001
	    upper, lower = self.LHS_disp(upperVal, k, radius), self.LHS_disp(lowerVal, k, radius)
	    if k ==0: # by inspection this must be true
	        return 1
	    for n in range(nstep):
	        middle = self.LHS_disp(middleVal, k, radius)
	        if (upper*middle < 0): 
	            lowerVal = middleVal
	            middleVal = (lowerVal + upperVal)*0.5
	            lower = middle
	        
	        else:
	            upperVal = middleVal
	            middleVal = (lowerVal + upperVal)*0.5
	            upper = middle
	    
	    return middleVal

	def k_from_omega(self, omega, r, nstep = 300, upper = 10, onlyShort = False):
	    sVal, kc = self.s(omega, r), self.k_crit(r)
	    k_list = np.linspace(0,upper*kc, nstep)
	    
	    
	    RHS = lambda k : sVal**2 + 2 * pi * self.Sigma(r)/(self.kappa(r)**2) * self.reduction_factor(sVal, self.chi(k, r))*abs(k)*exp(-abs(k) *self.epsilon)-1
	    lst = list(map(RHS, k_list))
	    
	    k = [k_list[i] for i in range(1, len(lst)) if lst[i]*lst[i-1] < 0]

	    if (len(k) ==1) & (onlyShort == True): # Returns message if we only find the shorter wave length
	        print("Only found short wavelength")
	    
	    return k

	def forbidden_radius(self, omega0, nstep=100): # Returns the inner forbidden radius/CR 
	    cr = self.CR(omega0)
	    radii = np.linspace(1.01*self.ILR(omega0), 0.99*cr, nstep)
	    k = list(map(self.k_from_omega, [omega0 for e in radii], radii))
	    r_plot = [radii[i]/cr for i, klst in enumerate(k) for void in klst]


	    return max(r_plot)


	def k_func_r(self, omega0, forbidden_radius = 0.99): # Use this for plotting r vs k
	    cr = self.CR(omega0)
	    radii = np.linspace(1.01*self.ILR(omega0), forbidden_radius*cr, 100)
	    k = list(map(self.k_from_omega, [omega0 for e in radii], radii, [1000 for e in radii]))
	    r_plot, k_plot = [radii[i]/cr for i, klst in enumerate(k) for void in klst], [kvalue for klst in k for kvalue in klst]
	    
	    return r_plot, k_plot

	def plotting_k_vs_r(self, omega0):
		fr = self.forbidden_radius(omega0)
		r, k = self.k_func_r(omega0, fr)
		plt.scatter(r, [k[i]/self.k_crit(r[i]) for i in range(len(r))], s =1)

		plt.axvline(fr, linestyle = '--', color = 'firebrick')
		plt.show()

	## Group Velocity ##
	## -------------- ##

	def stellar_vg(self, k, radius, deltaX = 0.001):
	    kc = self.k_crit(radius)
	    x0, x1, x2 = abs(k/kc)- deltaX, abs(k/kc), abs(k/kc) + deltaX
	    s0, s1, s2 = self.s_from_k(x0*kc, radius), self.s_from_k(x1*kc, radius), self.s_from_k(x2*kc, radius)
	    grad = -((1/(x2-x1)) * (abs(s2)-abs(s1)))
	    

	    
	    return (self.kappa(radius)/kc)* np.sign(k * s1) * grad

	

	## Groove Search for omega_0 ##
	## ------------------------- ##

	def integrate_k(self, omega0, r_inf, nstep = 100):
	    r_sup = self.forbidden_radius(omega0)
	    radii = self.CR(omega0) * np.linspace(r_inf, r_sup, nstep)
	    k = list(map(self.k_from_omega, [omega0 for e in radii], radii))
	    k_unpacked = [kvalue for klst in k for kvalue in klst]

	    return (sum(k_unpacked))*(radii[1] - radii[0]) # We need to integrate over both branches 
        

	def find_mode_omega0(self, inital_guess, r_inf, nstep = 20): 
	    
	    omegaL, omegaM, omegaU = inital_guess[0], sum(inital_guess)/len(inital_guess), inital_guess[1]
	    integralL = self.integrate_k(omegaL, r_inf/self.CR(omegaL)) - pi
	    integralU = self.integrate_k(omegaU, r_inf/self.CR(omegaU)) - pi
	    
	    for n in range(nstep):
	        integralM = self.integrate_k(omegaM, r_inf/self.CR(omegaM)) - pi
	        #print(omegaL, omegaU)
	        if (integralM*integralL < 0):
	            omegaU, integralU = omegaM, integralM
	            omegaM = 0.5 * (omegaL + omegaU)
	        else:
	            omegaL, integralL = omegaM, integralM
	            omegaM = 0.5 * (omegaL + omegaU)
	    
	    return omegaM;

	## Groove Search for eta ##
	## --------------------- ##

	def etaFromOmega0(self, omega0, scarRadius, deltaU = 0.01): # I think we need to be careful with this integration 
		rIn, rOut = self.innerRadius(omega0, scarRadius), self.outerRadius(omega0)

		integral = 0
		r = lambda u, rIn, rOut : rIn + (rOut - rIn) * (sin(u) **2)

		for u in np.arange(0, 0.5 * pi - deltaU, deltaU):
			R = r(u, rIn, rOut)
			integral += (1/self.groupVelocity(omega0, R, True) + 1/self.groupVelocity(omega0, R, False)) * (rOut - rIn) * sin(2*u) * deltaU
			print(R, (1/self.groupVelocity(omega0, R, False) + 1/self.groupVelocity(omega0, R, False)) * (rOut - rIn) * sin(2*u) * deltaU)

		return integral * (2 / log(2))
		

	def find_mode_eta(self, omega0, r_inf, nstep = 100):
	    r_sup = self.forbidden_radius(omega0)
	    radii = np.linspace(r_inf, self.CR(omega0) * r_sup, nstep)
	    k = list(map(self.k_from_omega, [omega0 for e in radii], radii, [750 for e in radii]))
	    r, k_plot   = [radii[i] for i, klst in enumerate(k) for void in klst], [kvalue for klst in k for kvalue in klst]
	    
	    vg_2_invert = list(map(self.stellar_vg, k_plot, r))
	    vg_inverted = [1/abs(v_g) for v_g in vg_2_invert]
	    '''plt.scatter(r, vg_2_invert, s = 1)
	    	    	    	    plt.show()'''
	    #return (2/log(2))*(sum(vg_inverted) - 0.5 * (max(vg_inverted) + min(vg_inverted))) * (radii[1] - radii[0])
	    return 0.25/((sum(vg_inverted) - 0.5 * (max(vg_inverted) + min(vg_inverted))) * (radii[1] - radii[0]))

	def correct_k(self, omega0, radius, bool_lst): # bool_lst = [isLong, isTrailing]
		sgn = 1 if bool_lst[1] == True else -1
		k_lst = self.k_from_omega(omega0, radius)
		#print(radius/self.ILR(omega0), k_lst, bool_lst)

		return sgn * min(k_lst) if bool_lst[0] == True else sgn * max(k_lst)

	def update_region(self, omega0, radius, bool_lst, fr):
		if radius > 0.999 * fr *self.CR(omega0): 
			if ((bool_lst[0] == bool_lst[1])):
				bool_lst[0] = not (bool_lst[0])
			

		if ((bool_lst[1] == False) & (bool_lst[0] == True) & (radius < 1.05 * self.ILR(omega0))): 
			bool_lst[1] = not bool_lst[1]

		return bool_lst




	def motion_of_wave(self, omega0): # Figure out why it is quite jumpu
		time, region, fr = np.linspace(0, 40, 300), [False, False], self.forbidden_radius(omega0)
		radius, k, vg, deltaT = np.zeros_like(time), np.zeros_like(time), np.zeros_like(time), time[1]
		radius[0] = 1.1*self.ILR(omega0) 
		k[0]  = self.correct_k(omega0, radius[0], region)
		vg[0] = self.stellar_vg(k[0], radius[0])
		print(radius[0], k[0], vg[0])
		
		for t in range(1, np.shape(time)[0]):
			radius[t] = radius[t-1] + vg[t-1] * deltaT 
			region = self.update_region(omega0, radius[t], region, fr)
			k[t]  = self.correct_k(omega0, radius[t], region)
			vg[t] = self.stellar_vg(k[t], radius[t])
			
			print(region, radius[t], k[t], np.sign(vg[t]))
			
		plt.plot([k[i]/self.k_crit(r) for i, r in enumerate(radius)], [r/self.CR(omega0) for r in radius])
		plt.axhline(fr)
		plt.axhline(self.ILR(omega0)/self.CR(omega0))
		plt.show()
	

## We need a function to select the correct k 

disc = Mestel_Disc(1/sqrt(12.4))

'''
disc.motion_of_wave(0.6)
k = disc.k_from_omega(0.6, 1.0739418023159923)
print(k, disc.stellar_vg(-k[1], 1.0739418023159923),disc.stellar_vg(k[1], 1.0739418023159923) )
'''


omega0 = disc.find_mode_omega0([0.5, 0.7], 1.2)
eta = disc.find_mode_eta(omega0, 1.2, 100) ## Find convergence in nstep for eta integral 

print(omega0, eta)


#disc.plotting_k_vs_r(0.547)
## Can we come up with a better way to figure out forbidden_radius so that if we use more points 

