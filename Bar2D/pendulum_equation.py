import numpy as np
import matplotlib.pyplot as plt
from math import *

from scipy.integrate import odeint

plt.rc('text', usetex=True)
plt.rc('font', family='serif')



def pendulum_func(y, t, omega = 1):
	return np.array([y[1], -sin(y[0])])


def solve_pedulum(time, IC):
	return odeint(pendulum_func, np.array(IC), time)
# Some function that gives the seperatirix

#plottings 
time = np.linspace(0,50, 500)
delta_t = time[1]



ICs = np.linspace(0.2, 2.7, 11)

for ic in ICs:
	evolution = solve_pedulum(time, [0, ic])
	color = 'royalblue' if ic < 2 else 'firebrick'
	if ic == ICs[0]:
		plt.plot(evolution[:,0]/pi,evolution[:,1], color = color, label = "Libration")
	if ic == ICs[-1]:
		plt.plot(evolution[:,0]/pi,evolution[:,1], color = color, label = "Rotation")
	plt.plot(evolution[:,0]/pi,evolution[:,1], color = color)
	plt.plot(-evolution[:,0]/pi,evolution[:,1], color = color)

	dx_arrow, dy_arrow = pendulum_func(evolution[0], 0)* delta_t
	x_arrow, y_arrow = evolution[0,0] - 0.5*dx_arrow, evolution[0,1]
	plt.arrow(x_arrow, y_arrow+0.01, dx_arrow, dy_arrow, width = 0.02, color = 'k', head_length = 0.02) 

	evolution = solve_pedulum(time, [0, -ic])
	
	plt.plot(evolution[:,0]/pi,evolution[:,1], color = color)
	plt.plot(-evolution[:,0]/pi,evolution[:,1], color = color)
	plt.arrow(-x_arrow, -y_arrow-0.01, -dx_arrow, -dy_arrow, width = 0.02, color = 'k', head_length = 0.02) 
	if ic == ICs[1]:
		evolution = solve_pedulum(time, [0, 1.999999])
		plt.plot(evolution[:,0]/pi,evolution[:,1], color = 'k', linestyle = '--', label = "Separatrix")	


plt.xlim([-1,1])
plt.ylim([-2.6,2.6])
plt.xlabel(r"$\theta/\pi$", fontsize = 15)
plt.ylabel(r"$J_{\theta}$", fontsize = 15)
plt.legend(fontsize = 12, loc = 'lower right')



plt.show()