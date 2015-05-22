# makefile for the raytracer
# from quickcf tutorial:  g++ *.cpp -lSDL -O3 

#F90 = gfortran
F90 = mpiifort
#vpath %.cpp ../source
#vpath %.h ../source
#VPATH = ../source

OBJECTS = 

#CFLAGS   = -I/home/buntemeyer/Downloads/netcdf/netcdf-fortran-4.2/f90
#CFLAGS   = -I${NETCDFF}/include

CFLAGS = -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-parallel-impi-intel14/include -I/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.0-impi-intel14/include -DgFortran -I/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-impi-intel14/include -I/sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-impi-intel14/include -I/sw/rhel6-x64/sys/libaec-0.3.2-intel14/include

#LFLAGS   = -L/home/buntemeyer/Downloads/netcdf/netcdf-fortran-4.2/f90 -lnetcdff 
#LFLAGS   = -L${NETCDFF}/lib -lnetcdff 
#LFLAGS   =  -lnetcdff 
LFLAGS = -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-parallel-impi-intel14/lib -lnetcdff -L/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.0-impi-intel14/lib -L/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-impi-intel14/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-impi-intel14/lib -lnetcdf -L/sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-impi-intel14/lib -Wl,-rpath,/sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-impi-intel14/lib -lhdf5 -lhdf5_hl -L/sw/rhel6-x64/sys/libaec-0.3.2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/sys/libaec-0.3.2-intel14/lib -lsz -lz -lcurl -Wl,-rpath,/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.0-impi-intel14/lib -lnetcdf

all: simple_xy_par_wr simple_xy_par_rd simple_xy_wr simple_xy_rd \
     pres_temp_4D_wr pres_temp_4D_rd \
     nc4_pres_temp_4D_wr

simple_xy_par_wr:  simple_xy_par_wr.o 
	$(F90) simple_xy_par_wr.o -o simple_xy_par_wr $(LFLAGS)

simple_xy_par_rd:  simple_xy_par_rd.o 
	$(F90) simple_xy_par_rd.o -o simple_xy_par_rd $(LFLAGS)

simple_xy_wr:  simple_xy_wr.o 
	$(F90) simple_xy_wr.o -o simple_xy_wr $(LFLAGS)

simple_xy_rd:  simple_xy_rd.o 
	$(F90) simple_xy_rd.o -o simple_xy_rd $(LFLAGS)

pres_temp_4D_wr:  pres_temp_4D_wr.o 
	$(F90) pres_temp_4D_wr.o -o pres_temp_4D_wr $(LFLAGS)

pres_temp_4D_rd:  pres_temp_4D_rd.o 
	$(F90) pres_temp_4D_rd.o -o pres_temp_4D_rd $(LFLAGS)

nc4_pres_temp_4D_wr:  nc4_pres_temp_4D_wr.o 
	$(F90) nc4_pres_temp_4D_wr.o -o nc4_pres_temp_4D_wr $(LFLAGS)


main.o : main.F90
	$(F90) -c main.F90

%.o: %.f90 
	$(F90) -c $(CFLAGS) $<

clean:
	rm *.o
