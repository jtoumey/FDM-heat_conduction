# Finite Difference Solution of a 1-D Boundary Value Problem for Heat Conduction
## Julian M. Toumey
## March 2015
* I wrote this code originally in `Matlab`, this version is a translation to `FORTRAN90`
* Tested with `gcc 4.6.3`

## Instructions:
1. Compile with: `gfortran HEATCONDUCT.f90 thomas.f90 calc_jacobian.f90 calc_lineareqn.f90 write_results.f90 -o heatconduct.o`
2. Run with `./heatconduct.o`
