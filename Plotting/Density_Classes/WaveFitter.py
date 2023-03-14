from Density_Classes.TwoDdensity import *
from statistics import median

from scipy.special import jn_zeros, jv, hankel1, hankel2 # Bessel Functions 
from scipy.fft import fft2, ifft2, fftfreq
from operator import gt, le

def grad(radius, phase, multiple = 3):
	rad, grad = [r1 for r1 in radius[1:]], [(ph1-ph0)/(r1-r0) for ph1, ph0, r1, r0 in zip(phase[1:], phase, radius[1:], radius)]
	r_o, g_o =[], []
	med = abs(median(grad))
	
	for r, g in zip(rad, grad):
		if abs(abs(g) - med)/med < multiple:
			r_o.append(r)
			g_o.append( g)

	return r_o, g_o

class WaveFitter():
	
	def __init__(self, filename, timeEnd, l = 2, rmax = 10): 
	
		self.density = TwoDdensity(filename)
		self.rmax, self.l = rmax, l
		self.time = np.linspace(0, timeEnd, self.density.nSteps)

		self.radii, self.amp, self.phase = self.density.fouierEvolution(self.l, self.rmax)
		self.r_step = self.radii[1] - self.radii[0]
		self.n_time_step = self.time.shape[0]

		#self.nmax = 100
		

	## Splitting Functions ## 
	## ------------------- ##

	def list_k_fitting(self):
		n_x, x_min, x_max = self.radii.shape[0], self.radii[0], self.radii[-1]
		return [ 2*pi*(n-n_x//2)/x_max for n in range(n_x) if (n!=n_x//2)] # We can't have k = 0 

	def forward_transform_k_time(self, k, timeIndex):
		integrand = np.exp(-1j * self.radii * k)*np.sqrt(np.abs(k)*self.radii) * self.amp[timeIndex, :] * np.exp(1j * self.phase[timeIndex,:]) ## Why do we have that factor of l? 
		return np.sum(integrand) * ((self.radii[1]-self.radii[0]))

	def forward_transform(self):
		k = self.list_k_fitting()
		k_projected = np.zeros((self.n_time_step, len(k)), dtype = complex)
		for t in range(self.n_time_step):
			k_projected[t, :] = np.fromiter(map(self.forward_transform_k_time, k, (t for _ in k)), dtype = complex)

		self.fTransformed_pos_l = fft2(k_projected, axes = (0))

	def split_components(self): # In-Clock, Out-Clock, Out-Anti, In-Anti ; anti has omega_p > 0
		self.fTransformed_pos_l_split = [np.zeros_like(self.fTransformed_pos_l) for _ in range(4)]
		oMid, kMid = self.fTransformed_pos_l.shape[0]//2, self.fTransformed_pos_l.shape[1]//2

		self.fTransformed_pos_l_split[0][:oMid,kMid:] = self.fTransformed_pos_l[:oMid,kMid:]
		self.fTransformed_pos_l_split[1][:oMid,:kMid] = self.fTransformed_pos_l[:oMid,:kMid]
		self.fTransformed_pos_l_split[2][oMid:,kMid:] = self.fTransformed_pos_l[oMid:, kMid:]
		self.fTransformed_pos_l_split[3][oMid:,:kMid] = self.fTransformed_pos_l[oMid:,:kMid]
	
	## Backwards Transforms - No K split ##
	## --------------------------------- ##

	def backward_transfrom_R_time(self, R, toTransform):
		k = np.asarray(self.list_k_fitting(), dtype = complex)
		integrand = np.exp(1j * k * R)*(1/np.sqrt(np.abs(k)*R)) * toTransform 
		return np.sum(integrand) * (abs(k[1]-k[0])/(2*pi))

	def backward_transform_split(self, toTransform):
		holding = np.zeros((self.n_time_step, self.radii.shape[0]), dtype = complex) 
		for time in range(self.n_time_step):
			holding[time, :] = np.fromiter(map(self.backward_transfrom_R_time, self.radii, (toTransform[time, :] for _ in self.radii)), dtype = complex)
		
		return ifft2(holding, axes = (0))

	def backward_transform(self):
		self.real_space_seperated = [self.backward_transform_split(each) for each in self.fTransformed_pos_l_split]

	def transform_and_split(self):
		self.forward_transform()
		self.split_components()
	
	def split_waves(self):
		self.transform_and_split()
		self.backward_transform()

	## In-going & Out-going Test ## 
	## --------------------------##

	def outgoing_T(self, time : int):
		return self.real_space_seperated[1][time, :]+ self.real_space_seperated[2][time, :]

	def  ingoing_T(self, time : int):
		return self.real_space_seperated[0][time, :]+ self.real_space_seperated[3][time, :]

	def total_T(self, time : int):
		total = self.ingoing_T(time) + self.outgoing_T(time)
		return np.absolute(total), np.angle(total)

	def check_k_fitting(self, time : int):
		fig,axs = plt.subplots(ncols = 2)
		axs[0].plot(self.radii, self.amp[time, :], color = 'firebrick', linestyle = "--")
		axs[1].plot(self.radii, self.phase[time, :], color = 'firebrick', linestyle = "--")

		self.split_waves()

		outgoing = self.outgoing_T(time)
		amp, phase = np.absolute(outgoing), np.angle(outgoing)
		#plt.plot(self.radii, amp)

		ingoing = self.ingoing_T(time)
		amp, phase = np.absolute(ingoing), np.angle(ingoing)
		#plt.plot(self.radii, amp)

		total = ingoing + outgoing 
		amp, phase = np.absolute(total), np.angle(total)
		axs[0].plot(self.radii, amp)
		axs[1].plot(self.radii, phase)
			#plt.plot(self.radii, amp[time, :]*np.cos(phase[time,:]), linestyle = '--')

		plt.xlabel("Radius")
		plt.show()	

	def splitting_Plot(self, time):
		self.split_waves()

		fig,axs = plt.subplots(ncols = 2, nrows = 2, sharex = True, sharey = 'row')
		axs[0,0].set_xlim([1.5, self.radii[-1]])
		inInd = 10

		axs[0,0].plot(self.radii[inInd:], np.absolute(self.ingoing_T(time))[inInd:], linestyle = '--')
		axs[0,0].plot(self.radii[inInd:], np.real(self.ingoing_T(time))[inInd:])
		axs[0,0].axhline(color = 'k', linestyle = 'dotted')

		axs[0,1].plot(self.radii[inInd:], np.absolute(self.outgoing_T(time))[inInd:], linestyle = '--')
		axs[0,1].plot(self.radii[inInd:], np.real(self.outgoing_T(time))[inInd:])
		axs[0,1].axhline(color = 'k', linestyle = 'dotted')


		axs[1,0].plot(*grad(self.radii[inInd:], np.angle(self.ingoing_T(time))[inInd:]))
		axs[1,0].axhline(color = 'k', linestyle = 'dotted')


		axs[1,1].plot(*grad(self.radii[inInd:], np.angle(self.outgoing_T(time))[inInd:]))
		axs[1,1].axhline(color = 'k', linestyle = 'dotted')
		
		axs[0,0].set_title("In-going")
		axs[1,0].set_title("Out-going")
	
		plt.show()


	## Backwards Transform - K Split ##
	## ----------------------------- ##
	def k_crit(self, radius, activeFraction = 0.5):
		return 2/ (activeFraction * radius)	

	def backward_transfrom_R_time_k(self, R, toTransform, k_tp, comp_func):
		kc = self.k_crit(R) * k_tp
		k_lst  = np.asarray(self.list_k_fitting(), dtype = complex)
		integrand = [np.exp(1j * k * R)*(1/np.sqrt(np.abs(k)*R)) * f_k for f_k, k in zip(toTransform, k_lst) if comp_func(k, kc)]
		return sum(integrand) * (abs(k_lst[1]-k_lst[0])/(2*pi))

	def backward_transform_split_k(self, toTransform, k_tp, comp_func):
		holding = np.zeros((self.n_time_step, self.radii.shape[0]), dtype = complex) 
		for time in range(self.n_time_step):
			holding[time, :] = np.fromiter(map(self.backward_transfrom_R_time_k, self.radii, (toTransform[time, :] for _ in self.radii), (k_tp for _ in self.radii), (comp_func for _ in self.radii)), dtype = complex)
		
		return ifft2(holding, axes = (0))

	def backward_transform_k(self, k_tp):
		fTransformed = [self.fTransformed_pos_l_split[0] + self.fTransformed_pos_l_split[3], self.fTransformed_pos_l_split[1] + self.fTransformed_pos_l_split[2]]
		comp_funcs = [gt, le]
		self.real_space_seperated = [self.backward_transform_split_k(wave, k_tp, cf) for wave in fTransformed for cf in comp_funcs] # The last for is looped through first
	

	def split_waves_k(self, k_tp):
		self.transform_and_split() 
		self.backward_transform_k(k_tp) # [IS, IL, OS, OL]

	## In-going & Out_going Splitting K Test ## 
	## ------------------------------------- ## 

	def ingoing_S(self, time : int):
		return self.real_space_seperated[0][time, :]

	def ingoing_L(self, time : int):
		return self.real_space_seperated[1][time, :]

	def outgoing_S(self, time : int):
		return self.real_space_seperated[2][time, :]

	def outgoing_L(self, time : int):
		return self.real_space_seperated[3][time, :]

	def total_K(self, time : int):
		total = self.ingoing_S(time) + self.ingoing_L(time) + self.outgoing_S(time) + self.outgoing_L(time)
		return np.absolute(total), np.angle(total)

	def check_k_splitting(self, time, k_tp):
		fig,axs = plt.subplots(ncols = 2)

		axs[0].plot(self.radii, self.amp[time, :], color = 'firebrick', linestyle = "--")
		axs[1].plot(self.radii, self.phase[time, :], color = 'firebrick', linestyle = "--")
		
		self.split_waves_k(k_tp)
		
		tot_amp, tot_phase = self.total_K(time)
		axs[0].plot(self.radii, tot_amp)	
		axs[0].set_xlabel("Radius")

		axs[1].plot(self.radii, tot_phase)
	

		plt.show()

	def k_splitting_Plot(self, time, k_tp):
		self.split_waves_k(k_tp)

		fig,axs = plt.subplots(ncols = 2, nrows = 4, sharex = True, sharey = 'row')
		axs[0,0].set_xlim([1.5, self.radii[-1]])
		inInd = 10

		axs[0,0].plot(self.radii[inInd:], self.radii[inInd:]*np.absolute(self.ingoing_S(time))[inInd:], linestyle = '--')
		axs[0,0].plot(self.radii[inInd:], self.radii[inInd:]*np.real(self.ingoing_S(time))[inInd:])
		axs[0,0].axhline(color = 'k', linestyle = 'dotted')

		axs[0,1].plot(self.radii[inInd:], self.radii[inInd:]*np.absolute(self.ingoing_L(time))[inInd:], linestyle = '--')
		axs[0,1].plot(self.radii[inInd:], self.radii[inInd:]*np.real(self.ingoing_L(time))[inInd:])
		axs[0,1].axhline(color = 'k', linestyle = 'dotted')

		axs[1,0].plot(self.radii[inInd:], self.radii[inInd:]*np.absolute(self.outgoing_S(time))[inInd:], linestyle = '--')
		axs[1,0].plot(self.radii[inInd:], self.radii[inInd:]*np.real(self.outgoing_S(time))[inInd:])
		axs[1,0].axhline(color = 'k', linestyle = 'dotted')

		axs[1,1].plot(self.radii[inInd:], self.radii[inInd:]*np.absolute(self.outgoing_L(time))[inInd:], linestyle = '--')
		axs[1,1].plot(self.radii[inInd:], self.radii[inInd:]*np.real(self.outgoing_L(time))[inInd:])
		axs[1,1].axhline(color = 'k', linestyle = 'dotted')

		axs[2,0].plot(*grad(self.radii[inInd:], np.angle(self.ingoing_S(time))[inInd:]))
		axs[2,0].axhline(color = 'k', linestyle = 'dotted')


		axs[2,1].plot(*grad(self.radii[inInd:], np.angle(self.ingoing_L(time))[inInd:]))
		axs[2,1].axhline(color = 'k', linestyle = 'dotted')

		
		axs[3,0].plot(*grad(self.radii[inInd:], np.angle(self.outgoing_S(time))[inInd:]))
		axs[3,0].axhline(color = 'k', linestyle = 'dotted')

		
		axs[3,1].plot(*grad(self.radii[inInd:], np.angle(self.outgoing_L(time))[inInd:]))
		axs[3,1].axhline(color = 'k', linestyle = 'dotted')


		
		
		axs[0,0].set_title("Short")
		axs[0,1].set_title("Long")
		
		axs[0,0].set_ylabel("In-going")
		axs[1,0].set_ylabel("Out-going")
		axs[2,0].set_ylabel("In-going")
		axs[3,0].set_ylabel("Out-going")
		plt.show()

		
	## Animation Functions ##
	## ------------------- ## 

	def phase_evolution(self, radius, filename = None):
		index = int(radius/self.r_step)

		fig, axs = plt.subplots(ncols = 2, sharex = True)

		for i, each in enumerate(self.real_space_seperated[2:]):
			axs[0].plot(self.time, np.angle(each[:, index]))
			axs[1].plot(self.time, np.abs(each[:, index]))

		axs[0].set_xlabel("Time")

		axs[0].set_ylabel("Phase")
		axs[1].set_ylabel("Amplitude")

		plt.ylabel("Phase")
		plt.title(f"Phase Evolution at $R = {radius}$")
		if filename == None:
			plt.show()
		else:
			ani.save(filename, writer = writer)
			print(f"Animation saved to: {filename}")


	def density_animation(self, filename = None):
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=20, metadata=dict(artist='Me'))

		fig, axs1 = plt.subplots(1,1)
		ims = []
		colors = ['firebrick', 'royalblue', 'royalblue', 'firebrick']
		style = ['-', '-', '--', '--']

		for time in range(self.n_time_step):
			items = []
			for each, c, s in zip(self.real_space_seperated, colors, style):
				item, = axs1.plot(self.radii[5:], np.real(each[time, 5:]), linestyle = s, color = c, label = "Magnitude", animated = True)
				items.append(item)
			
			items.append(fig.text(.4,.9,(r"Time: " +str(time))))
			ims.append([*items])


		ani = animation.ArtistAnimation(fig, ims, interval=30)

		if filename == None:
			plt.show()
		else:
			ani.save(filename, writer = writer)
			print(f"Animation saved to: {filename}")


	## Density Reconstruction ##
	## ---------------------- ## 

	def component_index(self, is_clockwise, is_ingoing):
		if is_clockwise ==True:
			if is_ingoing == True:
				return 0
			else: 
				return 1
		else:
			if is_ingoing:
				return 3
			else:
				return 2 

	def density_reconstuction(self, time, is_clockwise, is_ingoing): ## Check this with the orignal, then change to the correct array
		index = self.component_index(is_clockwise, is_ingoing)

		xA, yA = np.linspace(-self.rmax, self.rmax, 2*np.shape(self.radii)[0]), np.linspace(-self.rmax, self.rmax, 2*np.shape(self.radii)[0])
		array = np.zeros((xA.shape[0], yA.shape[0])) 

		for i, x in enumerate(xA):
			for j, y in enumerate(yA):
				R, phi = sqrt(x**2 + y**2), atan2(y, x)
				if R < (self.rmax-self.r_step): 
					index = R/self.r_step
					amp = (index - (index//1)) * self.amp[time][int(index +1)] + ((index+1)//1 - index) * self.amp[time][int(index)]
					phase = (index - (index//1)) * self.phase[time][int(index +1)] + ((index+1)//1 - index) * self.phase[time][int(index)]
					array[i,j] = amp * cos(2 * phi + phase)

		return array

	## Recombination Testing ## 
	## --------------------- ## 


	def fitted_at_time(self, time : int):
		return (self.real_space_seperated[2])[time, :], (self.real_space_seperated[3])[time, :]

	def reconstruction_at_time(self, time : int):
		a, b = self.fitted_at_time(time)
		amp_1, amp_2, phase_1, phase_2 = np.absolute(a), np.absolute(b), np.angle(a), np.angle(b)

		eq_1 = amp_1 * np.cos(phase_1) + amp_2 * np.cos(phase_2)
		eq_2 = amp_1 * np.sin(phase_1) + amp_2 * np.sin(phase_2)

		return np.sqrt(eq_1**2 + eq_2**2), np.arctan2(eq_2, eq_1)

	

	# Some function that does the reconsturction 
	# Some function that plots this along side the orginal 




