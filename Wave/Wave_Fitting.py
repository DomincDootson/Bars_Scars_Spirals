## This script decomposes a wave into a left and right moving part ##  

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy.fft import fft2,ifft2, rfftfreq, fftfreq, ifft
from math import *

def wave_packet(centre, x):
	return (1/sqrt(2*pi))*2*np.exp(-np.square((x-centre)/0.5))

def wave_animation(wave_packets1, wave_packets2, wave_packets3, x, t, c = 1):
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=20, metadata=dict(artist='Me'))

	fig, axs1 = plt.subplots(1,1)
	ims = []

	for i, time in enumerate(t):
		waves1, = axs1.plot(x, wave_packets1[i, :], color = 'firebrick', label = "Total Wave")
		waves2, = axs1.plot(x, wave_packets2[i, :], color = 'royalblue', label = "Left Moving Wave")
		waves3, = axs1.plot(x, wave_packets3[i, :], color = 'royalblue', linestyle = '--', label = "Right Moving Wave") 
		title = fig.text(.4,.9,f"Centre: {time*c:.2f}")
		ims.append([waves1, waves2, waves3,title])

	ani = animation.ArtistAnimation(fig, ims, interval=30)

	#ani.save("Waves_on_a_string.mp4", writer = writer) 
	plt.show()


def fft(wave_packets): 
	
	ft = fft2(wave_packets) #, axes = (1,0)
	return ft

def fft_plotting(wave_packets, x, t):
	ft = fft(wave_packets)
	ift = ifft(ft)
	
	xf, tf = fftfreq(np.shape(x)[0], abs(x[1]-x[0])), fftfreq(np.shape(t)[0],  abs(t[1]-t[0]))
	xv, tv = np.meshgrid(xf, tf)

	plt.contourf(xv, tv, np.log(np.abs(ft)), levels = 100)
	plt.colorbar()
	
	plt.show()

def splitting_freq(ft_wave):
	left = np.zeros_like(ft_wave)
	ft_wave_shifted = np.fft.fftshift(ft_wave)

	midO, midK = int(ft_wave.shape[0]//2), int(ft_wave.shape[1]//2)
	left[0:midO,0:midK] = ft_wave_shifted[0:midO,0:midK]
	left[midO:, midK:]  = ft_wave_shifted[midO:,midK:]

	right = ft_wave_shifted-left

	return np.fft.fftshift(left), np.fft.fftshift(right)
	

def inverting_waves(ft_left, ft_right):
	return np.real(ifft2(ft_left)), np.real(ifft2(ft_right))


def splitting_waves():
	x = np.linspace(-5, 5,1000)
	t = np.linspace(0,10,1000)
	c = 1

	waves = np.zeros((np.shape(t)[0], np.shape(x)[0]))

	for i, time in enumerate(t):
		waves[i, :] = exp(0.2*time)*(wave_packet(c * time-5, x) + wave_packet(-c * time+4, x))#wave_packet(-4+c*time, x) + wave_packet(+4-c*time , x)


	left, right = inverting_waves(*splitting_freq(fft(waves)))
	wave_animation(waves, left, right, x, t)


