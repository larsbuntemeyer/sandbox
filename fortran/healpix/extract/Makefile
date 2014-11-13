# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

F90 = ifort 
CC = gcc

#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source


OBJECTS = 

CFLAGS   = #-I/usr/include/x86_64-linux-gnu

LFLAGS   = 

FFLAGS  = -O3 

incpix = convert_inplace_1d_inc.f90 \
	convert_inplace_nd_inc.f90 \
	convert_nest2ring_nd_inc.f90 \
	convert_nest2ring_1d_inc.f90 \
	convert_ring2nest_nd_inc.f90 \
	convert_ring2nest_1d_inc.f90 \
	remove_dipole_inc.f90 \
	pixel_routines.F90

incind = indmed_part1.f90 \
         indmed_part2.f90 

obj0  = healpix_types.o indmed.o cgetEnvironment.o utilities.o 
obj1  = bit_manipulation.o extension.o long_intrinsic.o 
obj2  = misc_utils.o
obj3  = num_rec.o statistics.o

obj   = $(obj0) $(obj1) $(obj2) $(obj3) pix_tools.o

main:   main.o 
	$(F90) main.o $(obj) -o main

main.o : pix_tools.o

pix_tools.o: pix_tools.F90 $(incpix) $(obj)

indmed.o: indmed.f90 $(incind)

%.o: %.F90 
	$(F90) -c $(FFLAGS) $<

%.o: %.f90 
	$(F90) -c $(FFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<

clean:
	rm *.o *.mod 

