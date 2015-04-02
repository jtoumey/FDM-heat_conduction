heatconduct.o: HEATCONDUCT.o solve_tdma.o calc_jacobian.o calc_lineareqn.o write_results.o
	gfortran HEATCONDUCT.o solve_tdma.o calc_jacobian.o calc_lineareqn.o write_results.o -o heatconduct.o
HEATCONDUCT.o: HEATCONDUCT.f90 solve_tdma.f90 calc_jacobian.f90 calc_lineareqn.f90 write_results.f90
	gfortran -c HEATCONDUCT.f90 solve_tdma.f90 calc_jacobian.f90 calc_lineareqn.f90 write_results.f90
