all : 
	gfortran -g -O3 -fopenmp gif_util.f90 octtree.f90 snaptool.f90 -o snaptool