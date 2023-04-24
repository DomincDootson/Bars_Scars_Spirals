import matplotlib.pyplot as plt 
import numpy as np
from math import *

Omega, kappa, rT = 1, 1*sqrt(2), 1
B = -(kappa)**2/(4 * Omega)

def H_x(p_x, v_phi):
	return 0.5 * (p_x*p_x - (Omega/B)*v_phi*v_phi)


def x_bar(p_y):
	return 2 * (Omega/kappa**2)*(p_y - rT * Omega)

def v_phi(x, p_y):
	return 2 * B * (x - x_bar(p_y))


def moment(p_x, p_y, order):
	p_y = np.linspace(-5,5,1000)
	p_x = np.linspace(-5,5, 1000)
	PX, PY = np.meshgrid(p_x,p_y)

	return np.sum((PY**order)*np.exp((-H_x(PX,v_phi(x, PY))/(sig*sig))))/np.sum(np.exp((-H_x(PX,v_phi(x, PY))/(sig*sig))))


def potential(radius):
	

x, sig = 0, 0.1
p_y = np.linspace(-5,5,1000)
p_x = np.linspace(-5,5, 1000)

# PX, PY = np.meshgrid(p_x,p_y)
# plt.contourf(PY,PX, np.exp((-H_x(PX,v_phi(x, PY))/(2*sig*sig))), levels = 50)
# plt.xlabel(r"$p_{y}$")
# plt.ylabel(r"$p_{x}$")
# plt.colorbar()

print(sqrt(Omega/(-B))*sqrt(moment(p_x, p_y, 2) - moment(p_x, p_y, 1)**2))


plt.show()