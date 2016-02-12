# ======= Makefile for downscaling program ============

EXE = downscale.exe

# first link one numrec.f file needed into this directory
#  ln -s ../num_rec/numrec.f

#FC = gfortran
FC = pgf90

# flags for debugging or for maximum performance, comment as necessary

ifeq ($(FC),pgf90)
	FCFLAGS = -Mbackslash -fPIC
#	FCFLAGS = -g -Mbackslash -fPIC
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

	# libraries needed for linking, unused in the examples
        LDLOCAL = -Wl,-rpath,/usr/local.net/netcdf-4.3.0/gcc-4.7.2/lib
	LDFLAGS = $(LDLOCAL) -L/usr/local.net/netcdf-4.3.0/gcc-4.7.2/lib -lnetcdf -lnetcdff
	#LDFLAGS = -L/usr/local/netcdf/lib -lnetcdff
endif


# "make" builds all
all: $(EXE)

$(EXE):	type.o \
	precmod.o \
	stringmod.o \
	numrec.o \
	utim.o \
	read_config.o \
	read_control.o \
	read_nc.o \
	read_nc_grid.o \
	write_nc.o \
	coefficients.o \
	main.o

	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^ 

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD
