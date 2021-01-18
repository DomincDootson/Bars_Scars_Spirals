CXX = g++
CXXFLAGS = -O3 -std=c++2a -Wall 

HEADERFILES = Action_Angle_Basis_Functions/ActionAngleBasisContainer.h Potential_Density_Pair_Classes/KalnajsBasis.h Potential_Density_Pair_Classes/GaussianLogBasis.h Potential_Density_Pair_Classes/PotentialDensityPairContainer.h DF_Class/Mestel.h \
		Volterra_Solver/VolterraSolver.h Volterra_Solver/EvolutionKernels.h Volterra_Solver/ExpansionCoeff.h 

all : BasisFunctions 

BasisFunctions : main.o physics.o
	$(CXX) $(CXXFLAGS) -o main main.o physics.o





main.o : main.cpp $(HEADERFILES)
	$(CXX) $(CXXFLAGS) -c main.cpp

physics.o : physics.cpp $(HEADERFILES)
	$(CXX) $(CXXFLAGS) -c physics.cpp

clean : 
	@rm -f *.o