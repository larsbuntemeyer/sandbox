# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

F90 = gfortran
F90 = mpif90

#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source

OBJECTS = 

CFLAGS   = -I/home/buntemeyer/Downloads/netcdf/netcdf-fortran-4.2/f90

LFLAGS   = -L/home/buntemeyer/Downloads/netcdf/netcdf-fortran-4.2/f90 -lnetcdff 

simple_xy_par_wr:  simple_xy_par_wr.o 
	$(F90) simple_xy_par_wr.o -o simple_xy_par_wr $(LFLAGS)

simple_xy_wr:  simple_xy_wr.o 
	$(F90) simple_xy_wr.o -o simple_xy_wr $(LFLAGS)

simple_xy_rd:  simple_xy_rd.o 
	$(F90) simple_xy_rd.o -o simple_xy_rd $(LFLAGS)

main.o : main.F90
	$(F90) -c main.F90

%.o: %.f90 
	$(F90) -c $(CFLAGS) $<

clean:
	rm $(OBJECTS) 
