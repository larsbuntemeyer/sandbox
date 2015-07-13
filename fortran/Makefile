# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

F90 = ifort 

#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source

OBJECTS = 

CFLAGS   = 

LFLAGS   = 

main:   main.o 
	$(F90) main.o -o main

main.o : main.F90
	$(F90) -c main.F90

%.o: %.F90 
	$(F90) -c $<

clean:
	rm $(OBJECTS) main main.o 
