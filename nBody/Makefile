# First give all the flags for the complier

CXX = g++
CXXFLAGS = -O3 -std=c++2a -Wall 
INCLUDES =  -I/../DF_Class -I/../Potential_Density_Pair_Classes -I/../Volterra_Solver -I/../Physics_Functions -I/nBody_Classes

# Label all the different files that we want to compile and how they link


NBODY = nBody
NBODY_SRC = nBody.cpp 
NBODY_OBJ = nBody.o 


# Define all the rules

all : $(NBODY)


$(NBODY): $(NBODY_OBJ)
	$(CXX) $(CXXFLAGS) -lfftw3 -lm  -o $@ $^ 


# Define the rule for making .o from .cpp files 
%.o: %.cpp ../Bar2D/*.h ../DF_Class/*.h ../Potential_Density_Pair_Classes/*.h ../Volterra_Solver/ExpansionCoeff.h nBody_Classes/*.h 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 

#nBody.o: nBody.cpp ../Bar2D/*.h ../DF_Class/*.h ../Potential_Density_Pair_Classes/*.h ../Volterra_Solver/ExpansionCoeff.h nBody_Classes/*.h 
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -c nBody.cpp 