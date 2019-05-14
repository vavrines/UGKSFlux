# Description
## Flux function of unified gas-kinetic scheme for 1D and 2D cases
The unified gas-kinetic scheme (UGKS) is a numerical solver of kinetic-type (Boltzmann) equations. Despite the demo program given in https://github.com/vavrines/UGKS, sometimes you just want to omit the details of implementation and directly take advantage of the algorithm. This repo serves that purpose where two Fortran90 modules are provided, with 1D/2D UGKS flux functions implemented. Thus, direct call of the subroutines is enabled. 

## For the source file
UGKSFlux1D.f90 is for 1D flux evaluation  
UGKSFlux2D.f90 is for 2D flux evaluation  

Although only 1D and 2D solvers are provided here, based on the same numerical principle, its extension to 3D case is straightforward.

# Pre-requirements
A Fortran compiler is need which supports Fortran90 or later versions, e.g. ifort or gfortran.

# Usage

`use UGKSFlux1D`
`call flux_ugks1d(wL, hL, bL, shL, sbL, lenL, fluxw, fluxh, fluxb, wR, hR, bR, shR, sbR, lenR, unum, uspace, weight, ink, gamma, muref, omega, prandtl, dt)`

The Makefile is provided for those who have gnu make installed. If you are using any IDE (e.g. Visual Studio), use the compiling function provided by the IDE.

Note: openmp and Intel Fortran compiler is used by default. DO NOT type the prompt symbol $


# License
Copyright (C) 2019 Tianbai Xiao tianbaixiao@gmail.com

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. To use this program, you are supposed to receive a copy of the GNU General Public License. See the GNU General Public License for more details. http://www.gnu.org/licenses/
