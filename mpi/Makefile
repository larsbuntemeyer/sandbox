# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

#F90 = mpif90-vt
F90 = mpiifort
#F90 = vtfort -vt:fc mpif90
#F90 = vtf90 -vt:f90 mpif90
#F90 = /opt/openmpi-1.7.2/bin/mpif90
#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source

OBJECTS = 

CFLAGS   = -r8 -c

FFLAGS_OPT   =  -c -r8 -i4 -m64 -mcmodel=large -shared-intel 

LFLAGS_OPT   = -i4 -r8 -mcmodel=large -shared-intel -L/usr/lib64
#LFLAGS_OPT   = -i4 -r8 -mcmodel=large -i-dynamic -L/usr/lib64 

#LFLAGS_OPT   =-L/usr/lib/x86_64-redhat-linux5E/lib64 -L/usr/lib64 -L-m64 -lpthread -lrt -lssl -lcrypto -lappllib -lnsl -lpthread -lresolv -lrt -lm -ldl -lmp

main:   main.o 
	$(F90) $(LFLAGS_OPT) main.o -o main

main.o : main.F90
	$(F90) $(FFLAGS_OPT) main.F90

%.o: %.F90 
	$(F90) $(CFLAGS) $<

clean:
	rm $(OBJECTS) main main.o 
