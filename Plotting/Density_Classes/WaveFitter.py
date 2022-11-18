from Density_Classes.TwoDdensity import *

from scipy.special import jn_zeros, jv, hankel1, hankel2 # Bessel Functions 
from scipy.fft import fft2, ifft2, fftfreq

class WaveFitter():
	
	def __init__(self, filename, timeEnd):
	
		self.density = TwoDdensity(filename)
		self.rmax, self.l = 10, 2
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
		integrand = np.exp(-1j * self.radii * k)*np.sqrt(np.abs(k)*self.radii) * self.amp[timeIndex, :] * np.exp(1j * self.l * self.phase[timeIndex,:])
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

	def split_singular_k(self, k_index): # In-Clock, Out-Clock, Out-Anti, In-Anti ; anti has omega_p > 0
		self.fTransformed_pos_l_split = [np.zeros_like(self.fTransformed_pos_l) for _ in range(4)]
		oMid, kMid = self.fTransformed_pos_l.shape[0]//2, self.fTransformed_pos_l.shape[1]//2

		for k in k_index:
			if k < kMid:
				self.fTransformed_pos_l_split[0][:oMid,k] = self.fTransformed_pos_l[:oMid,k]
				self.fTransformed_pos_l_split[2][oMid:,k] = self.fTransformed_pos_l[oMid:,k]
				
			if k >= kMid:	
				self.fTransformed_pos_l_split[1][:oMid,k] = self.fTransformed_pos_l[:oMid,k]
				self.fTransformed_pos_l_split[3][oMid:,k] = self.fTransformed_pos_l[oMid:,k]
		

	def backward_transfrom_R_time(self, R, toTransform):
		k = np.asarray(self.list_k_fitting(), dtype = complex)
		integrand = np.exp(1j * k * R)*(1/np.sqrt(np.abs(k)*R)) * toTransform 
		return np.sum(integrand) *(abs(k[1]-k[0])/(2*pi))

	def backward_transform_split(self, toTransform):
		holding = np.zeros((self.n_time_step, self.radii.shape[0]), dtype = complex) 
		for time in range(self.n_time_step):
			holding[time, :] = np.fromiter(map(self.backward_transfrom_R_time, self.radii, (toTransform[time, :] for _ in self.radii)), dtype = complex)
		
		return ifft2(holding, axes = (0))#

	def backward_transform(self):
		self.real_space_seperated = [self.backward_transform_split(each) for each in self.fTransformed_pos_l_split]
		

	def split_waves(self):
		self.forward_transform()
		self.split_components()
		self.backward_transform()

	def split_waves_singular_k(self, k_index : list):
		self.forward_transform()
		self.split_singular_k(k_index) 
		# fig, axs = plt.subplots(ncols = 4)
		# axs[0].imshow(np.abs(self.fTransformed_pos_l_split[0]))
		# axs[1].imshow(np.abs(self.fTransformed_pos_l_split[1]))
		# axs[2].imshow(np.abs(self.fTransformed_pos_l_split[2]))
		# axs[3].imshow(np.abs(self.fTransformed_pos_l_split[3]))
		# plt.show()
		self.backward_transform() 
		
		
	## Animation Functions ##
	## ------------------- ## 

	def phase_evolution(self, radius, filename = None):
		index = int(radius/self.r_step)

		fig, axs = plt.subplots(ncols = 2)

		for i, each in enumerate(self.real_space_seperated):
			axs[0].plot(self.time, np.angle(each[:, index]))
			axs[1].plot(self.time, np.abs(each[:, index]))

		plt.xlabel("Time")
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