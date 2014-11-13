# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

F90 = ifort 

#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source

OBJECTS = 

CFLAGS   = -c -check all

LFLAGS   = 

main:   main.o 
	$(F90) main.o -o main

%.o: %.F90 
	$(F90) $(CFLAGS) $<

clean:
	rm $(OBJECTS) main main.o 
