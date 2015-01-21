# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

#F90 = ifort 
F90 = gfortran

#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source

OBJECTS = test_module_interface.o testA.o code11.o 

CFLAGS   = -c 

LFLAGS   = 

main:   main.o $(OBJECTS) 
	$(F90) main.o $(OBJECTS) -o main

main.o : main.F90 $(OBJECTS)
	$(F90) -c main.F90

%.o: %.F90 
	$(F90) $(CFLAGS) $<

%.o: %.f
	$(F90) $(CFLAGS) $<

clean:
	rm $(OBJECTS) main main.o *.mod 
