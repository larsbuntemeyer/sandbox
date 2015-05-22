# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

#F90 = gfortran
F90 = mpiifort
EXE = nonblocking_write

#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source

OBJECTS = 

#CFLAGS   = -I/home/buntemeyer/Downloads/netcdf/netcdf-fortran-4.2/f90
#CFLAGS   = -I${NETCDFF}/include

CFLAGS = -I/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.0-impi-intel14/include

#LFLAGS   = -L/home/buntemeyer/Downloads/netcdf/netcdf-fortran-4.2/f90 -lnetcdff 
#LFLAGS   = -L${NETCDFF}/lib -lnetcdff 
#LFLAGS   =  -lnetcdff 
LFLAGS = -L/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.0-impi-intel14/lib -lpnetcdf

$(EXE):  $(EXE).o 
	$(F90) $(EXE).o -o $(EXE) $(LFLAGS)

%.o: %.f90 
	$(F90) -c $(CFLAGS) $<

clean:
	rm $(OBJECTS) 
