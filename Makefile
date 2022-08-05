# First give all the flags for the complier

CXX = g++
CXXFLAGS = -O3 -std=c++2a -Wall
INCLUDES = -I/usr/include/eigen3  -I/Action_Angle_Basis_Functions -I/DF_Class -I/Potential_Density_Pair_Classes -I/Volterra_Solver -I/Physics_Functions -I/DF_Function -I/nBody/Box -I/usr/local/include/eigen3 -I/usr/local/include/boost/math/special_functions -I/Response_Matrix

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

MODES = modes
MODES_SRC = $(FILE_DIR)/modes.cpp $(FILE_DIR)/generalFunctions.cpp 
MODES_OBJ = $(FILE_DIR)/modes.o  $(FILE_DIR)/modesFunctions.o





SCARRED = scarred
SCARRED_SRC = $(FILE_DIR)/scarredEvolution.cpp $(FILE_DIR)/scarredEvolutionFunctions.cpp 
SCARRED_OBJ = $(FILE_DIR)/scarredEvolution.o  $(FILE_DIR)/scarredEvolutionFunctions.o

SPIRAL = spiral
SPIRAL_SRC = $(FILE_DIR)/spiralEvolution.cpp $(FILE_DIR)/spiralEvolutionFunctions.cpp 
SPIRAL_OBJ = $(FILE_DIR)/spiralEvolution.o  $(FILE_DIR)/spiralEvolutionFunctions.o

STABILITY = stability
STABILITY_SRC = $(FILE_DIR)/parameterStability.cpp
STABILITY_OBJ = $(FILE_DIR)/parameterStability.o

WAVES = waves
WAVES_SRC = $(FILE_DIR)/waves.cpp $(FILE_DIR)/wavesFunctions.cpp 
WAVES_OBJ = $(FILE_DIR)/wavesFunctions.o  $(FILE_DIR)/waves.o


# Define all the rules

all : $(BAR) $(DENSITY) $(GENERAL) $(MODES) $(SCARRED) $(SPIRAL) $(STABILITY) $(WAVES)

$(BAR): $(BAR_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(DENSITY): $(DENSITY_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(GENERAL): $(GENERAL_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(MODES): $(MODES_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(SCARRED): $(SCARRED_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(SPIRAL): $(SPIRAL_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(STABILITY): $(STABILITY_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

$(WAVES): $(WAVES_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ 


# Define the rule for making .o from .cpp files 
%.o: %.cpp Action_Angle_Basis_Functions/*.h Bar2D/*.h DF_Class/*.h DF_Function/*.h nBody/Box/Box.h Potential_Density_Pair_Classes/*.h Response_Matrix/*.h Spiral2d/*.h Volterra_Solver/*.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $< 