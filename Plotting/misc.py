from generalPlottingFunctions import *

def plotting_kalnajs_BF():	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	data = readingInRealCSV("BF_Comparison/KalnajsFunctions.csv")

	fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = 'row')

	axs[0,0].plot(data[:,0], data[:,1], label = r"$\mathcal{U}^{\ell}_{0}(R)$", color = "navy")
	axs[0,0].plot(data[:,0], data[:,3], label = r"$\mathcal{U}^{\ell}_{1}(R)$", color = "royalblue")
	axs[0,0].plot(data[:,0], data[:,5], label = r"$\mathcal{U}^{\ell}_{2}(R)$", color = "cornflowerblue")

	axs[1,0].plot(data[:,0], data[:,2], label = r"$\mathcal{D}^{\ell}_{0}(R)$", color = "firebrick")
	axs[1,0].plot(data[:,0], data[:,4], label = r"$\mathcal{D}^{\ell}_{1}(R)$", color = "red")
	axs[1,0].plot(data[:,0], data[:,6], label = r"$\mathcal{D}^{\ell}_{2}(R)$", color = 'lightcoral')

	axs[0,1].plot(data[:,0], data[:,7], color = "navy")
	axs[0,1].plot(data[:,0], data[:,9], color = "royalblue")
	axs[0,1].plot(data[:,0], data[:,11], color = "cornflowerblue")

	axs[1,1].plot(data[:,0], data[:,8], color = "firebrick")
	axs[1,1].plot(data[:,0], data[:,10], color = "red")
	axs[1,1].plot(data[:,0], data[:,12], color = "lightcoral")

	axs[0,2].plot(data[:,0], data[:,13], color = "navy")
	axs[0,2].plot(data[:,0], data[:,15], color = "royalblue")
	axs[0,2].plot(data[:,0], data[:,17], color = "cornflowerblue")

	axs[1,2].plot(data[:,0], data[:,14], color = "firebrick")
	axs[1,2].plot(data[:,0], data[:,16], color = "red")
	axs[1,2].plot(data[:,0], data[:,18], color = "lightcoral")


	axs[0,0].set_ylabel("Potential", fontsize = 15)
	axs[1,0].set_ylabel("Density", fontsize = 15)

	axs[1,0].set_xlabel(r"$R/R_{Ka}$", fontsize = 15)
	axs[1,1].set_xlabel(r"$R/R_{Ka}$", fontsize = 15)
	axs[1,2].set_xlabel(r"$R/R_{Ka}$", fontsize = 15)

	axs[0,0].set_title(r"$\ell=0$", fontsize = 15)
	axs[0,1].set_title(r"$\ell=1$", fontsize = 15)
	axs[0,2].set_title(r"$\ell=2$", fontsize = 15)

	axs[0,0].legend(fontsize = 12)
	axs[1,0].legend(fontsize = 12)


	plt.show()


plotting_kalnajs_BF()