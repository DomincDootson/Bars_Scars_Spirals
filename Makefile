# First give all the flags for the complier

CXX = g++
CXXFLAGS = -O3 -std=c++2a -Wall 
INCLUDES = -I/Action_Angle_Basis_Functions -I/DF_Class -I/Potential_Density_Pair_Classes -I/Volterra_Solver -I/Physics_Functions -I/DF_Function

# Label all the different files that we want to compile and how they link

FILE_DIR = Physics_Functions

BAR = bar
BAR_SRC = $(FILE_DIR)/barEvolution.cpp $(FILE_DIR)/barEvolutionFunctions.cpp 
BAR_OBJ = $(FILE_DIR)/barEvolution.o $(FILE_DIR)/barEvolutionFunctions.o 


DENSITY = density
DENSITY_SRC = $(FILE_DIR)/densityEvolution.cpp $(FILE_DIR)/densityEvolutionFunctions.cpp 
DENSITY_OBJ = $(FILE_DIR)/densityEvolution.o $(FILE_DIR)/densityEvolutionFunctions.o

GENERAL = general
GENERAL_SRC = $(FILE_DIR)/general.cpp $(FILE_DIR)/generalFunctions.cpp 
GENERAL_OBJ = $(FILE_DIR)/generalFunctions.o  $(FILE_DIR)/general.o

STABILITY = stability
STABILITY_SRC = $(FILE_DIR)/parameterStability.cpp
STABILITY_OBJ = $(FILE_DIR)/parameterStability.o


# Define all the rules

all : $(BAR) $(DENSITY) $(GENERAL) $(STABILITY)

$(BAR): $(BAR_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(DENSITY): $(DENSITY_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(GENERAL): $(GENERAL_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(STABILITY): $(STABILITY_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 


# Define the rule for making .o from .cpp files 
%.o: %.cpp Action_Angle_Basis_Functions/*.h Bar2D/*.h DF_Class/*.h DF_Function/*.h Potential_Density_Pair_Classes/*.h Volterra_Solver/*.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 