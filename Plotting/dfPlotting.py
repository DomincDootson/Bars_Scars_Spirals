from generalPlottingFunctions import *
from Density_Classes.TwoDdensity import *


density= twoDdensity("densityEvolution.csv")
density.densityEvolution(10)


df = DFClass("dfApoPer.csv")
df.apoDFEvolution()

