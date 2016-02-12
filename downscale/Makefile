#
#
#

FC = pgf90

# flags for debugging or for maximum performance, comment as necessary

ifeq ($(FC),pgf90)
	FCFLAGS = -g -Mbackslash -fPIC


	# flags forall (e.g. look for system .mod files, required in gfortran)
	FCFLAGS += -I/opt/netcdf4-pgi/include

	# libraries needed for linking, unused in the examples
	LDLOCAL = -rpath,/opt/netcdf4-pgi/lib
	LDFLAGS = $(LDLOCAL) -L/opt/netcdf4-pgi/lib -lnetcdff -lnetcdf

endif

ifeq ($(FC),gfortran)

	# flags forall (e.g. look for system .mod files, required in gfortran)

	# libraries needed for linking, unused in the examples
        LDLOCAL = -Wl,-rpath,/usr/local.net/netcdf-4.3.0/gcc-4.7.2/lib
	LDFLAGS = $(LDLOCAL) -L/usr/local.net/netcdf-4.3.0/gcc-4.7.2/lib -lnetcdf -lnetcdff
endif


# "make" builds all
all: downscale_test.exe


downscale_test.exe: 	type.o \
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
	mv *exe ../bin

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod *.MOD
