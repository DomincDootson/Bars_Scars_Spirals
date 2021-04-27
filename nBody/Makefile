# First give all the flags for the complier

CXX = g++
CXXFLAGS = -O3 -std=c++2a -Wall 
INCLUDES =  -I/../DF_Class -I/../Potential_Density_Pair_Classes -I/../Volterra_Solver -I/../Physics_Functions \
 -I/Bodies -I/Box -I/NBody_Classes -I/Perturbation_Grids

# Label all the different files that we want to compile and how they link


NBODY = nBody
NBODY_SRC = nBody.cpp nBodyFunctions.cpp
NBODY_OBJ = nBody.o nBodyFunctions.o 


# Define all the rules

all : $(NBODY) particleSampling

$(NBODY): $(NBODY_OBJ)
	$(CXX) $(CXXFLAGS) -lfftw3 -lm  -o $@ $^ 

particleSampling : particleSampling.o 
	$(CXX) $(CXXFLAGS) -o particleSampling particleSampling.o

particleSampling.o : Bodies/particleSampling.cpp Bodies/Bodies.h ../DF_Class/*.h  
	 $(CXX) $(CXXFLAGS) $(INCLUDES) -c Bodies/particleSampling.cpp

# Define the rule for making .o from .cpp files 
%.o: %.cpp ../Bar2D/*.h ../DF_Class/*.h ../Potential_Density_Pair_Classes/*.h ../Volterra_Solver/ExpansionCoeff.h \
Bodies/*.h Box/*.h NBody_Classes/*.h  Perturbation_Grids/*.h nBodyFunctions.h 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 

#nBody.o: nBody.cpp ../Bar2D/*.h ../DF_Class/*.h ../Potential_Density_Pair_Classes/*.h ../Volterra_Solver/ExpansionCoeff.h nBody_Classes/*.h 
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -c nBody.cpp 