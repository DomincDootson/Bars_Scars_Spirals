## This script impliments projecting onto the plane wave cylindrical basis ## 
from math import *
import numpy as np
import matplotlib.pyplot as plt


R_max = 10
x = np.linspace(0.011,R_max, 101)
y = np.exp(-0.5*np.square((x-R_max*0.4)/(R_max*0.1)))

def generate_basis(x, k):
	return np.exp(1j * x * k)/np.sqrt(abs(x*k))

## Generate list of ks
def list_k(x):
	n_x, x_min, x_max = x.shape[0], x[0], x[-1]

	return [ 2*pi*(n-n_x//2)/x_max for n in range(n_x) if (n!=n_x//2)] # We can't have k = 0 

## Do forward Projection 
def forward_transform_k(k, y, x):
	integrand = np.exp(-1j * x * k)*np.sqrt(np.abs(k)*x) * y
	return np.sum(integrand) * ((x[1]-x[0]))

def forward_transform(y, x):
	k = list_k(x)
	return np.fromiter(map(forward_transform_k, k, (y for _ in k), (x for _ in k)), dtype = complex)
	

## Backwards Projection ##

def backward_transfrom_R(R, y, k):
	integrand = np.exp(1j * k * R)*(1/np.sqrt(np.abs(k)*R)) * y
	return np.sum(integrand) *(abs(k[1]-k[0])/(2*pi))

def backward_transform(y, x):
	k = np.asarray(list_k(x))


	return np.fromiter(map(backward_transfrom_R, x, (y for _ in x), (k for _ in x)), dtype = complex)


k1, k2 = list_k(x)[45], list_k(x)[62]

f = y#generate_basis(x, k1) + generate_basis(x, k2) 

plt.plot(x, np.real(f))
plt.plot(x, np.imag(f))


f_k = forward_transform(f, x)
f_o = backward_transform(f_k, x)

plt.plot(x, np.real(f_o))
plt.plot(x, np.imag(f_o))


plt.show()