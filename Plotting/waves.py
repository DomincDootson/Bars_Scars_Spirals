from generalPlottingFunctions import *
from Density_Classes.OneDdensity import *
from Density_Classes.TwoDdensity import *
from Density_Classes.WaveFitter import * 
from CoefficientClass import *

import matplotlib.animation as animation

def modesComparison(oneDdensity, filename = None, index2Normalise = -1):
	modePower = [each.modePower() for each in oneDdensity]

	if index2Normalise != -1:
		for each in modePower:
			each[:, index2Normalise] -= each[0,index2Normalise]

	bins, minValue, maxValue = [range(each.nMax) for each in oneDdensity], min([np.amin(each) for each in modePower]), max([np.max(each) for each in modePower])

	fig, axs = plt.subplots()
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=20, metadata=dict(artist='Me'))

	color, alpha = ['firebrick', "royalblue"], [1, 0.8]


	def animate(time):
		axs.clear()
		axs.set_title("Time = " + str(time*0.2))
		axs.set_ylim([minValue, maxValue])
		
		for i in range(len(modePower)):
			axs.bar(bins[i], modePower[i][time, :], color = color[i], alpha = alpha[i], label =len(bins[i])-1)

		axs.legend()
		return axs

	ani = animation.FuncAnimation(fig, animate, frames = 250)
	if (filename):
		ani.save(filename, writer = writer)
		print("Animation saved to: " + filename) 
	else:
		plt.show()


# coeff = [CoeffficientClass("Waves_Data/Self_Consistent_0_30.csv"), CoeffficientClass("Waves_Data/Self_Consistent_0_Small_30.csv")]
# modesComparison(coeff, filename = "Waves_Plots/Waves_Videos/SelfConsistent_Comparison_30.mp4", index2Normalise = 30)

# coeff = [CoeffficientClass("Waves_Data/Self_Consistent_0_10.csv"), CoeffficientClass("Waves_Data/Self_Consistent_0_Small_10.csv")]
# modesComparison(coeff, filename = "Waves_Plots/Waves_Videos/SelfConsistent_Comparison_10.mp4", index2Normalise = 10)



def saveAnimations(): 
	oneDdensity = OneDdensity("Waves_Data/RadiusPull_Density_0_2.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPull_Density_0_2.mp4")
	oneDdensity = OneDdensity("Waves_Data/RadiusPull_Density_0_5.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPull_Density_0_5.mp4")
	oneDdensity = OneDdensity("Waves_Data/RadiusPull_Density_0_10.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPull_Density_0_10.mp4")
	oneDdensity = OneDdensity("Waves_Data/RadiusPull_Density_0_13.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPull_Density_0_13.mp4")

	oneDdensity = OneDdensity("Waves_Data/RadiusPing_Density_0_2.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPing_Density_0_2.mp4")
	oneDdensity = OneDdensity("Waves_Data/RadiusPing_Density_0_5.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPing_Density_0_5.mp4")
	oneDdensity = OneDdensity("Waves_Data/RadiusPing_Density_0_10.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPing_Density_0_10.mp4")
	oneDdensity = OneDdensity("Waves_Data/RadiusPing_Density_0_13.csv")
	oneDdensity.densityAnimation(remove_ic = True, write2file = "Waves_Plots/Waves_Videos/RadiusPing_Density_0_13.mp4")

def decomposed1Danimations(readStem, writeStem, file2read):
	for f2r in file2read:
		wavefitter = WaveFitter(readStem + f2r + '.csv', 10)
		wavefitter.split_waves()
		
		write2 = writeStem + f2r +".mp4"
		wavefitter.density_animation(filename = write2)


#saveAnimations()

#decomposed1Danimations("Waves_Data/Stirring/", "Waves_Plots/Waves_Videos/Stirring/", [f"CR_5_OLR_{w}_{d}" for w, d  in zip(['W', 'W', 'R', 'R'], ['P', 'N', 'P', 'N'])])
#twoD.fourierAnimations(2, rMax = 10, filename ="Another_vid_4_ali.mp4 2")

#den = TwoDdensity("Waves_Data/Stirring/CR_5_OLR_W_N.csv")


fitter = WaveFitter("Waves_Data/Stirring/CR_5_OLR_W_N.csv", 10)
fitter.split_waves_singular_k([49])
fitter.density_animation()




#fitter.density_animation()
#den.fourierAnimations(2)