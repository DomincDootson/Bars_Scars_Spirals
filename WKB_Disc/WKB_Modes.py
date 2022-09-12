from WKB_Disc import *

def readingInRealCSV(filename): # Reads in the basis functions
	with open(filename) as csv_file:
		csv_reader = csv.reader(csv_file, delimiter = ',')
		data = []
		for row in csv_reader:
			lst = [float(i) for i in row]
			data.append(lst)

		return np.asarray(data, dtype=np.float64)

def deRijkeMode():
	disc = WKB_Disc(1/sqrt(12.4))

	print(disc.modeFinder([0.4, 0.9], 1)) 




def varyingSoftening(softeningValues, innerPosition = 1.2):
	omega0 = []

	for ep in softeningValues:
		print(f"Finding mode for {ep:.2f}")
		disc = WKB_Disc(1/sqrt(12.4), epsilon = ep)
		omega0.append(disc.modeFinder([0.4, 0.7], innerPosition))
		print(omega0)

	fig, axs = plt.subplots(ncols = 2)

	axs[0].plot(softeningValues, [each.real for each in omega0])
	
	axs[0].set_xlabel(r"Softening, $\epsilon$")
	axs[0].set_ylabel(r"$\omega_{0}$")

	axs[1].plot(softeningValues, [each.imag for each in omega0]) 
	axs[1].set_xlabel(r"Softening, $\epsilon$")
	axs[1].set_ylabel(r"$\eta$")
	plt.show()
	# ep = softeningValues[-1]
	# disc = WKB_Disc(1/sqrt(12.4))
	# print(disc.find_mode_eta(0.5657197, 1.2, 200)) 
	# print(disc.find_mode_eta(0.5657197, 1.2, 500)) 
	# print(disc.find_mode_eta(0.5657197, 1.2, 1000)) 
	# print(disc.find_mode_eta(0.5657197, 1.2, 2000)) 

	# (0.5657197475433351+0.007992813861886824j)
	# 0.017550837684988975
	
	
def varyingInnerPositon(innerPosition, filename = "Cavity_Modes/VaryingInnerPosition.csv"):
	omega0, rInner = [], []
	f = open(filename, 'w')
	for rIn in innerPosition:
		print(f"Finding mode for {rIn:.2f}")
		densityFile = f"Disc_Density/Tapered_R_{str(int(rIn*10))}_W_25_D_-95_G.csv"
		print("Density file: " + densityFile)
		disc = WKB_Disc(1/sqrt(12.4), epsilon = 0, densityFile = densityFile)
		omega0.append(disc.modeFinder([0.4, 0.9], rScar = rIn))
		rInner.append(disc.forbidden_radius(omega0[-1].real) * disc.CR(omega0[-1].real))
		f.write(f"{rIn},{omega0[-1].real},{omega0[-1].imag}\n")
		print(f"{rIn},{omega0[-1].real},{omega0[-1].imag}\n")


	f.close()

	
def plottingInnerPosition(files = ["Cavity_Modes/VaryingInnerPosition.csv", "Cavity_Modes/VaryingInnerPosition_Tapered.csv", "Cavity_Modes/VaryingInnerPosition_Tapered_Scarred.csv"], scatterFiles = []):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	
	fig, axs = plt.subplots(ncols = 1)
	labels, colors, ls = ["Untapered","Tapered", "Tapered Scarred"], ["royalblue", "firebrick", "firebrick"], ['-', '--', '-']
	for file, l,c,s in zip(files, labels, colors, ls):
		data = readingInRealCSV(file)
		axs.plot(data[:,0]+0.125, data[:,1], label = l, color = c, linestyle = s)
		#axs[1].plot(data[:,0], data[:,2]) 

	
	for file in scatterFiles:
		data = readingInRealCSV(file)
		axs.scatter(data[1:,0], data[1:, 1], color = 'firebrick')
		#axs[1].scatter(data[:,0], data[:, 2]) 		

	axs.legend(fontsize = 12, title = "Density Model", title_fontsize=12)
	axs.set_xlabel(r"Inner Scar Radius, $R_{i}$", fontsize = 12)
	axs.set_ylabel(r"$\omega_{0}$", fontsize = 12)

	
	#axs[1].set_xlabel(r"Inner Scar Radius")
	#axs[1].set_ylabel(r"$\eta$")

	plt.show() 

def varyingOuterPositon(outerPosition, innerPosition = 1.2):
	omega0, rOuter = [], []

	for rO in outerPosition:
		print(f"Finding mode for {rO:.2f}")
		disc = WKB_Disc(1/sqrt(12.4), epsilon = 0)
		omega0.append(disc.modeFinder([0.2, 2], rScar = innerPosition, rUpperScar = rO))
		rOuter.append(disc.forbidden_radius(omega0[-1].real) * disc.CR(omega0[-1].real))
		print(omega0[-1])

	fig, axs = plt.subplots(ncols = 3)

	axs[0].plot(outerPosition, [each.real for each in omega0])
	axs[0].set_xlabel(r"Outer Scar Radius")
	axs[0].set_ylabel(r"$\omega_{0}$")

	axs[1].plot(outerPosition, [each.imag for each in omega0]) 
	axs[1].set_xlabel(r"Outer Scar Radius")
	axs[1].set_ylabel(r"$\eta$")

	axs[2].plot(outerPosition, rOuter)
	axs[2].plot(outerPosition, [rOuter[i] / disc.CR(omega0[i].real) for i in range(len(omega0))]) 
	axs[2].set_xlabel(r"Outer Scar Radius")
	axs[2].set_ylabel(r"Outer Position")

	plt.show()

def WKBapproximation():
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	fig, axs = plt.subplots()
	disc = WKB_Disc(1/sqrt(12.4), densityFile = "Disc_Density/Tapered_R_20_W_25_D_-95_G.csv", epsilon = 0)
	innerScar = disc.k_func_r(0.62)
	outerScar = disc.k_func_r(0.39)

	innerR = [innerScar[0][i] for i in range(0,len(innerScar[0]), 2)] + [innerScar[0][i] for i in range(len(innerScar[0])-1, 0, -2)]
	innerk = [innerScar[1][i] for i in range(0,len(innerScar[1]), 2)] + [innerScar[1][i] for i in range(len(innerScar[1])-1, 0, -2)]

	axs.plot([i -1 for i in innerR], [k*r for k, r in zip(innerR, innerk)], color = "royalblue", label = r"$1.0$")
	

	outerR = [outerScar[0][i] for i in range(0,len(outerScar[0]), 2)] + [outerScar[0][i] for i in range(len(outerScar[0])-1, 0, -2)]
	outerk = [outerScar[1][i] for i in range(0,len(outerScar[1]), 2)] + [outerScar[1][i] for i in range(len(outerScar[1])-1, 0, -2)]

	axs.plot([i-2 for i in outerR], [k*r for k, r in zip(outerR, outerk)], color ="firebrick", label = r"$2.0$")
	

	axs.legend(title = r"Scar Radius, $R_{i}$", fontsize = 12)

	axs.set_xlabel(r"$R-R_{i}$", fontsize = 12)
	axs.set_ylabel(r"$kR$", fontsize = 12)
	axs.set_ylim([0, 17.5])
	axs.axvline(0, linestyle ='--', color ="black")
	
	plt.show()

def forbiddenRadii(rScar, omega0 = []):
	
	for o, rIn in zip(omega0, rScar):
		densityFile = f"Disc_Density/Tapered_R_{str(int(rIn*10))}_W_25_D_-95_G.csv"
		disc =  WKB_Disc(1/sqrt(12.4), epsilon = 0, densityFile = densityFile)
		print(disc.CR(o) * disc.forbidden_radius(o))


#deRijkeMode()



#varyingSoftening(np.linspace(0.08, 0.10, 3))
#varyingInnerPositon(np.linspace(1, 2, 11), "Cavity_Modes/VaryingInnerPosition_Tapered_Scarred.csv")
#varyingInnerPositon(np.linspace(1, 2, 11))
#varyingOuterPositon(np.linspace(1.5, 2, 1))

#plottingInnerPosition(scatterFiles = ["Cavity_Modes/VaryingInnerPositionResponse_25_-95.csv"])
forbiddenRadii(rScar = [1.2, 1.4, 1.6, 1.8, 2.0], omega0 = [0.564407, 0.510169, 0.469492, 0.428814, 0.394915])
#deRijkeMode()

'''disc = WKB_Disc(1/sqrt(12.4), epsilon = 0)
disc.plotting_k_vs_r(0.40453948974609377, scars = [2])'''


# disc = WKB_Disc(1/sqrt(12.4), densityFile = "Disc_Density/Tapered_R_20_W_25_D_-95_G.csv", epsilon = 0)
# disc.modeFinder([0.35, 0.55], rScar = 2)
