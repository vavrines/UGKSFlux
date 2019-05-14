# Description
## Flux function of unified gas-kinetic scheme for 1D and 2D cases
The unified gas-kinetic scheme (UGKS) is a numerical solver of kinetic-type (Boltzmann) equations. Despite the demo program given in https://github.com/vavrines/UGKS, sometimes you just want to omit the details of implementation and directly take advantage of the algorithm. This repo serves that purpose where two Fortran90 modules are provided, with 1D/2D UGKS flux functions implemented. Thus, direct call of the subroutines is enabled. 

## For the source file
1. UGKSFlux1D.f90 is for 1D flux evaluation  
2. UGKSFlux2D.f90 is for 2D flux evaluation  

Although only 1D and 2D solvers are provided here, based on the same numerical principle, its extension to 3D case is straightforward.

# Requirement
A Fortran compiler is need which supports Fortran90 or later versions, e.g. ifort or gfortran.

# Usage
To use the module for 1D case, you need first to use it,  
`use UGKSFlux1D`  
and call the subroutine,  
`call flux_ugks1d(wL, hL, bL, shL, sbL, lenL, fluxw, fluxh, fluxb, wR, hR, bR, shR, sbR, lenR, unum, uspace, weight, ink, gamma, muref, omega, prandtl, dt)` 

For 2D case, the things are quite similar.  
`use UGKSFlux2D`  
`call flux_ugks2d(wL, hL, bL, shL, sbL, lenL, fluxw, fluxh, fluxb, lenFace, wR, hR, bR, shR, sbR, lenR, unum, vnum, uspace, vspace, weight, ink, gamma, muref, omega, prandtl, dt, dirc, cosa, sina)`

The parameters required stand for the conservative variables, particle distribution function and its slopes inside the cells next to a interface, their geometric features, interface flux functions, velocity quadrature points and weights, gas properties and time step.
The detailed meanings and usage of the variables required can be found in the annotations of the source code.

# Example
Two examples for using this program are provided, i.e. numerical simulations of 1D normal shock structure and 2D lid-driven cavity. To run them, first create a working directory, and then copy the files ./example/1D shock structure/ or ./example/2D lid-driven cavity/ into the new directory.  
Make the program through  
`$ Make`  
and run the program  
`$./bin/CloudCFD$`

The Makefile is provided for the users who have gnu make installed. If you are using some IDEs (e.g. Visual Studio), please use the compiling function provided by the IDE.

# License
Copyright (C) 2019 Tianbai Xiao tianbaixiao@gmail.com

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. To use this program, you are supposed to receive a copy of the GNU General Public License. See the GNU General Public License for more details. http://www.gnu.org/licenses/
