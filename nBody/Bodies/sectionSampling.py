import matplotlib.pyplot as plt

# we assume that vc =1 for this code

from math import *
VC = 1

def potential(r):
	return log(r)*VC

def energy(vr, angMom, radius):
	return 0.5*(vr**2 + (angMom**2)/(radius**2)) + potential(radius)

def jacobiIntegral(energy, angMom, patternSpeed):
	return energy - angMom*patternSpeed

def circularOrbitParams(radius):
	angMom = radius * VC
	return angMom, energy(0, angMom, radius) 



#First find the set of ICs for the circular orbit --> from this get the energy 
# take steps to decrease the angular momentum 


def vR(r, E, J):
	if (2*(E-potential(r)) - (J/r)**2 < 0): 
		print("Unphysical Value of vr")
	return sqrt(abs(2*(E-potential(r)) - (J/r)**2))

def radialMomentum(r, E, J):
	return (2 * (E - potential(r)) - (J/r)**2)

def rMin(E, J):
	x1, x2, nStep  = 0, exp(E-0.5), 20

	for i in range(nStep):
		holding = x1 + 0.5*(x2-x1)
		if radialMomentum(holding, E, J) <0:
			x1 = holding
		else:
			x2 = holding

	return x2 


def rMax(E, J):
	x1, x2, nStep  = exp(E-0.5), exp(E),  20

	for i in range(nStep):
		holding = x1 + 0.5*(x2-x1)
		if radialMomentum(holding, E, J)  > 0:
			x1 = holding
		else:
			x2 = holding

	return x1

def pow(x, y):
	x**y

def energyFuncRad(rApo, rPer):
	return (rApo*rApo*potential(rApo) -(rPer*rPer)*potential(rPer)) / (rApo*rApo - rPer*rPer)


def angMomFuncRad(rApo, rPer):
	return sqrt(2*(potential(rApo) - potential(rPer))/ (1/(rPer*rPer) - 1/(rApo*rApo)))


def lstEnergies(rCirc, patternSpeed):
	jCirc, eCirc = circularOrbitParams(rCirc)
	jacobi = 1.2*(eCirc - jCirc*patternSpeed)
	print("The jacobi Integral is: ", jacobi)

	print(jCirc, eCirc, jacobi)
	step, nOrbits = 0.02, 25
	j = [jCirc-(i+1)*step for i in range(nOrbits)]

	f = open("particleSamplesSections.out", "w+")
	f.write(str(nOrbits))
	x, y = [], []
	for i in range(nOrbits):
		#print(rMin(eCirc, j[i]), rMax(eCirc, j[i]))
		#print(rCirc, 0, vR(rCirc, eCirc, j[i]), j[i]/rCirc, rMin(eCirc, j[i]), rMax(eCirc, j[i])
		f.flush()
		#f.write('\n' +"9.82208e-05" + " " + str(rCirc) + " " + str(0) + " " +  str(vR(rCirc, eCirc, j[i])) + " " + str(j[i]/rCirc))
		
		#jacobi = eCirc - j[i]*patternSpeed
		vRadial = vR(rCirc, jacobi + j[i]*patternSpeed, j[i])

		f.write('\n' +"9.82208e-05" + " " + str(rCirc) + " " + str(0) + " " +  str(vRadial) + " " + str(j[i]/rCirc))

		e = energy(vRadial, j[i], rCirc)
		#e = energy(vR(rCirc, eCirc, j[i]), j[i], rCirc)
		y.append(jacobiIntegral(e, j[i], patternSpeed))
		x.append(j[i])
		print(rMin(jacobi + j[i]*patternSpeed, j[i]), rMax(jacobi + j[i]*patternSpeed, j[i]))


	f.close()
	plt.plot(x, y)
	plt.xlabel("Ang Mom")
	plt.ylabel("Jacobi Integral")

	plt.show()
'''
nOrbits, step = 5, .1
rCirc, patternSpeed = 2.1, 0.5
jCirc, energy = circularOrbitParams(rCirc)

j = [jCirc-(i+1)*step for i in range(nOrbits)]

for angMom in j:
	#energy = circularOrbitParams(rCirc)
	jacobi = energy - patternSpeed * angMom



	vrE, vrJ = vR(rCirc, energy, angMom), vR(rCirc, jacobi + patternSpeed * angMom, angMom)

	print(vrE, vrJ)
'''



lstEnergies(2.1, 0.5)
