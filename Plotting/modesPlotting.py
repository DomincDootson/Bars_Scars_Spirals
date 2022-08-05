from Density_Classes.TwoDdensity import *
from generalPlottingFunctions import *
from CoefficientClass import *

def eigenModeComparison(filename, patternSpeed):
	mode = TwoDdensity(filename)

	plt.imshow(mode.densityAtTime(-1))
	plt.show()


#eigenModeComparison("Modes_Data/JB_mode.csv", 1)

coef = CoeffficientClass("Modes_Data/JB_mode_Coeff.csv")
coef.powerPlot(100)

#Modes_Data/JB_mode_Coeff.csv