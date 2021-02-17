CXX = g++
CXXFLAGS = -O3 -std=c++2a -Wall 

HEADERFILES = Action_Angle_Basis_Functions/ActionAngleBasisContainer.h Potential_Density_Pair_Classes/KalnajsBasis.h Potential_Density_Pair_Classes/GaussianLogBasis.h Potential_Density_Pair_Classes/PotentialDensityPairContainer.h DF_Class/Mestel.h \
		Volterra_Solver/VolterraSolver.h Volterra_Solver/EvolutionKernels.h Volterra_Solver/ExpansionCoeff.h DF_Class/DFClass.h Action_Angle_Basis_Functions/ActionAngleBasisFunction.h Potential_Density_Pair_Classes/PotentialDensityPair.h Bar2D/Bar2D.h


all : BasisFunctions 

BasisFunctions : main.o physics.o
	$(CXX) $(CXXFLAGS) -o main main.o physics.o





main.o : main.cpp Physics_Functions/barEvolutionFunctions.h Physics_Functions/densityEvolutionFunctions.h Physics_Functions/generalFunctions.h
	$(CXX) $(CXXFLAGS) -c main.cpp


#physics.o : physics.cpp $(HEADERFILES)
#	$(CXX) $(CXXFLAGS) -c physics.cpp

generalFunctions.o : Physics_Functions/generalFunctions.cpp $(HEADERFILES)
	$(CXX) $(CXXFLAGS) -c Physics_Functions/generalFunctions.cpp

densityEvolutionFunctions.o : Physics_Functions/densityEvolutionFunctions.cpp $(HEADERFILES)
	$(CXX) $(CXXFLAGS) -c Physics_Functions/densityEvolutionFunctions.cpp

barEvolutionFunctions.o : Physics_Functions/barEvolutionFunctions.cpp $(HEADERFILES)
	$(CXX) $(CXXFLAGS) -c Physics_Functions/barEvolutionFunctions.cpp

clean : 
	@rm -f *.o