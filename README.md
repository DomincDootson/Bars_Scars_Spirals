# Bars, Scars and Spirals 
This folder contains the main back end for my PhD work studying the linear evolution of disc galaxies. For more information about the mathematical framework behind the code, please see my thesis that I have included in this document (`Bars_Stars_Spirals.pdf`).  

In this readme, I will run over the different objects in this file.


## Code Overview

### Potential-Density Pairs
An important part of my approach was to project onto a set of basis functions, the so-called 'potential-density' pairs. This folder contains their virtual base classes and the implementations of two specific sets of basis functions: the Kalnajs BF and the Gaussian BF. 

This folder contains two abstract classes 
- `PotentialDensityPair`: This class contains the potential density pair for a given index. 
- `PotentialDensityPairContainer`: This class wraps together a set of basis functions, with different indices (usually the radial index). It also contains functions for outputting the density field given a set of coefficients. 

The other files define the specific implementations of these virtual classes. 

### Action Angle Basis Functions 
The next step is to convert the BF, which are defined in real space, to action angled (AA) coordinates. The files in this folder do that. Again, they contain a specific class and a container class
- `ActionAngleBasisFunction`: A class that contains all the indice information and a 2D grid
- `ActionAngleBasisContainer`: A class that wraps together all the AA BF and defines the methods that do the transformations

### DF Class
This file contains a virtual distribution function and examples that inherit from it (e.g. Mestel Disc). As the DF sources the potential, it also contains a potential method as well as other methods that are specific to the DF.

### Volterra Solver
This folder contains the key methods for generating the response:
- `ExpansionCoeff`: This class contains the basis coefficients at each point in time. 
- `EvolutionKernel`: Contains the methods to calculate the evolution kernel at each point in time
- `VolterraSolver`: Pulls together the expansion coefficients, the evolution kernel and the perturbation and solves the Volterra equation to generate the response. 

### Bar2d
This file defined a perturbation that rotated at a fixed speed, that we used to represent a bar.

### Response Matrix
This file allowed us to generate the response matrix:

-`ResponseMatrix.h`: A class that had all the methods to calculate the response matrix either via direct integrate or from the time domain
-`ResponseMatrixThread.h`: A class that allowed for the RM to be evaluated at different omega via multithreading
-`ResponseMatrix.py`: Code that took a response matrix and solve for its eigenvalues (and hence the modes of the discs)

### nBody
Contains code that run N body simulations for comparison 