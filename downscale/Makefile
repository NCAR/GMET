# ============= Makefile for downscaling program ======================
# AWW-2016:  split mega-file coefficients.f90 into several smaller ones

# NOTE: first link one numrec.f file needed into this directory if not present
#  eg, ln -s ../from_NR/numrec.f

EXE = downscale.exe

FC = gfortran
#FC = pgf90
#FC = ifort

# flags for debugging or for maximum performance, comment as necessary

ifeq ($(FC),pgf90)
#	FCFLAGS = -Mbackslash -fPIC

	FCFLAGS = -g -Mbackslash -fPIC -Mlarge_arrays
#	FCFLAGS = -g -Mbackslash -fPIC -Ktrap=fp,inexact,inv,unf,ovf -Kieee
#	FCFLAGS = -O2 -Mbackslash -fPIC

        # flags forall (e.g. look for system .mod files, required in gfortran)
#	FCFLAGS += -I/usr/local/netcdf4-pgi/include
#	FCFLAGS += -I/opt/netcdf4/include
        FCFLAGS += -I/opt/netcdf4-pgi/include
#	FCFLAGS += -I/opt/netcdf-4.3.3.1/pgi-15.1/include

        # libraries needed for linking, unused in the examples
#	LDLOCAL = -rpath,/usr/local/netcdf-4.1.1/pgi-v11-7/lib
#	LDLOCAL = -rpath,/opt/netcdf4/lib
	LDLOCAL = -rpath,/opt/netcdf4-pgi/lib

#	LDFLAGS = $(LDLOCAL) -L/usr/local/netcdf4-pgi/lib -lnetcdff -lnetcdf
#	LDFLAGS = -L/usr/local/netcdf/lib -lnetcdff
#	LDFLAGS = $(LDLOCAL) -L/opt/netcdf4/lib -lnetcdff -lnetcdf
	LDFLAGS = $(LDLOCAL) -L/opt/netcdf4-pgi/lib -lnetcdff -lnetcdf
endif

ifeq ($(FC),gfortran)
#	FCFLAGS = -g
        #FCFLAGS = -O2

        # flags forall (e.g. look for system .mod files, required in gfortran)
#	FCFLAGS += -I/usr/local.net/netcdf-4.3.0/gcc-4.7.2/include
        FCFLAGS += -I/opt/netcdf4-gnu/include


        # libraries needed for linking, unused in the examples
        LDLOCAL = -Wl,-rpath,/opt/netcdf-4.7.0/gcc-6.3.0/lib
#     LDFLAGS = $(LDLOCAL) -L/usr/local.net/netcdf-4.3.0/gcc-4.7.2/lib -lnetcdf -lnetcdff
	LDFLAGS = $(LDLOCAL) -L/opt/netcdf4-gnu/lib -lnetcdf -lnetcdff
        #LDFLAGS = -L/usr/local/netcdf/lib -lnetcdff
endif

ifeq ($(FC),ifort)
#        FCFLAGS = -debug -warn all -check all -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0 -heap_arrays
        FCFLAGS = -O2 -warn all -FR -auto -WB -fltconsistency -fpe1 -heap_arrays
        

#       FCFLAGS = -g -Mbackslash -fPIC -Ktrap=fp,inexact,inv,unf,ovf -Kieee
        #FCFLAGS = -O2 -Mbackslash -fPIC

        # flags forall (e.g. look for system .mod files, required in gfortran)
#       FCFLAGS += -I/usr/local/netcdf4-pgi/include
#       FCFLAGS += -I/opt/netcdf4-ifort/include

        NCDF_PATH = /opt/netcdf4-intel
        LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf  -L/opt/hdf5/lib 
        INCNETCDF = -I$(NCDF_PATH)/include

        # libraries needed for linking, unused in the examples
#       LDLOCAL = -rpath,/usr/local/netcdf-4.1.1/pgi-v11-7/lib
#       LDFLAGS = $(LDLOCAL) -L/usr/local/netcdf4-pgi/lib -lnetcdff -lnetcdf

        LDFLAGS = $(LIBNETCDF) $(INCNETCDF) 

endif


# "make" builds all
all: $(EXE)

$(EXE):	type.o \
	precmod.o \
	string_mod.o \
	numrec.o \
	utim.o \
	read_config.o \
	read_control.o \
	read_nc.o \
	read_domain_grid.o \
	write_nc.o \
        stats_routines.o \
        regression_routines.o \
        compute_station_weights.o \
        read_station_weights.o \
        estimate_coefficients.o \
        estimate_forcing_regression.o \
	main.o

	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^ 

%.o: %.f90
	$(FC) $(FCFLAGS) $(LDFLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD *genmod.f90
