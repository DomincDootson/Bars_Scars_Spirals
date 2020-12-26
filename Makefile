CXX = g++
CXXFLAGS = -O3 -std=c++2a -Wall 

all : BasisFunctions

BasisFunctions : main.o Potential_Density_Pair_Classes/KalnajsBasis.o Potential_Density_Pair_Classes/GaussianLogBasis.o 
	$(CXX) $(CXXFLAGS) -o main main.o Potential_Density_Pair_Classes/KalnajsBasis.o Potential_Density_Pair_Classes/GaussianLogBasis.o 





main.o : main.cpp Action_Angle_Basis_Functions/ActionAngleBasisContainer.h Potential_Density_Pair_Classes/KalnajsBasis.h Potential_Density_Pair_Classes/PotentialDensityPairContainer.h DF_Class/Mestel.h \
		Volterra_Solver/VolterraSolver.h Volterra_Solver/EvolutionKernels.h Volterra_Solver/ExpansionCoeff.h 
	$(CXX) $(CXXFLAGS) -c main.cpp

KalnajsBasis.o : Potential_Density_Pair_Classes/KalnajsBasis.cpp Potential_Density_Pair_Classes/KalnajsBasis.h 
	$(CXX) $(CXXFLAGS) -c Potential_Density_Pair_Classes/KalnajsBasis.cpp

GaussianLogBasis.o : Potential_Density_Pair_Classes/GaussianLogBasis.cpp Potential_Density_Pair_Classes/GaussianLogBasis.h 
	$(CXX) $(CXXFLAGS) -c Potential_Density_Pair_Classes/GaussianLogBasis.cpp

clean : 
	@rm -f *.o