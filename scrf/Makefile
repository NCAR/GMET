#========================================================================
# Makefile to compile the SCRF Code from Martyn Clark
#========================================================================
#  Original Makefile from Martyn's FUSE Model
#  Modified:
#  Andy Newman
#
#
#========================================================================
# PART 0: Define directory paths
#========================================================================

# Define core directory below which everything resides
F_MASTER = 

# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)/scrf

# Location of the compiled modules
MOD_PATH = $(F_MASTER)/scrf

# Define the directory for the executables
EXE_PATH = $(F_MASTER)/bin
#
#========================================================================
# PART 1: Assemble all of the various sub-routines
#========================================================================

# Define directories
DRIVER_DIR = $(F_KORE_DIR)/driver
NUMREC_DIR = $(F_MASTER)/num_rec
MOD_DIR    = $(F_KORE_DIR)/share

#
# Numerical Recipes utilities
SCRF_NRUTIL= \
		nrtype.f90 \
		nrutil.f90 \
		nr.f90 \
		ran_state.f90
NRUTIL = $(patsubst %, $(NUMREC_DIR)/%, $(SCRF_NRUTIL))

#
# Numerical recipes procedures
SCRF_NRPROC= \
		ran1.f90 \
		gser.f90 \
		gcf.f90 \
		gamdev.f90 \
		gammln.f90 \
		gammp.f90 \
		gammq.f90 \
		lubksb.f90 \
		ludcmp.f90 \
		erf.f90 \
		erfc.f90 \
		erfcc.f90\
		gasdev.f90
NRPROC = $(patsubst %, $(NUMREC_DIR)/%, $(SCRF_NRPROC))


# Data modules
SCRF_DATAMS= \
		gridweight.f90 \
		linkstruct.f90 \
		dat_2dgrid.f90 \
		inputdat2d.f90 \
		qpe_namelist.f90
DATAMS = $(patsubst %, $(MOD_DIR)/%, $(SCRF_DATAMS))

# Routines taken from or modifed versions of Jason Craig's Prob QPE code
prob_qpe= \
		precmod.f90 \
		stringmod.f90 \
		save_precip.f90 \
		save_vars.f90 \
		read_nc_grid.f90 \
		read_grid_list.f90
qpe = $(patsubst %, $(F_KORE_DIR)/%, $(prob_qpe))

# Define routines for SCRF code
SCRF_GEN = \
		trig_degrees.f90 \
		erfinv.f90 \
		exit_scrf.f90 \
		field_rand.f90 \
		spcorr_grd.f90
GEN = $(patsubst %, $(F_KORE_DIR)/%, $(SCRF_GEN))


# ... and stitch it all together...
SCRF_ALL = $(NRUTIL) $(NRPROC) $(DATAMS) $(qpe) $(GEN)


#========================================================================
# PART 2: Define the libraries, driver programs, flags, and executables
#========================================================================

# Define the Fortran Compiler
FC  = pgf90

# Define the NetCDF libraries and path to include files
ifeq "$(FC)" "gfortran"
 LIBNETCDF = -L/opt/netcdf-4.1.3/gcc-4.4.5/lib -lnetcdff -lnetcdf -L/opt/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/opt/lib -lsz -lz
 INCNETCDF = -I/opt/netcdf-4.1.3/gcc-4.4.5/include
endif
ifeq "$(FC)" "ifort"
 NCDF_PATH = /usr/local/netcdf-4.1.3/ifort-12.0.5
 LOCL_PATH = /usr/local
 LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf -L$(LOCL_PATH)/lib -lsz -lz
 INCNETCDF = -I$(NCDF_PATH)/include
endif
ifeq "$(FC)" "pgf90"
 NCDF_PATH = /opt/netcdf4-pgi
 LOCL_PATH = /usr/local.net/
 LIBNETCDF = -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf -L$(LOCL_PATH)/lib -lz
 INCNETCDF = -I$(NCDF_PATH)/include
endif

# Define the driver routine
SCRF_DRIVER= \
	generate_ensembles.f90
DRIVER = $(patsubst %, $(DRIVER_DIR)/%, $(SCRF_DRIVER))

# Define the executable
DRIVER__EX = generate_ensemble.exe

# Define flags
ifeq "$(FC)" "gfortran"
 FLAGS = -static -Wall -g -ffree-line-length-none
endif
ifeq "$(FC)" "ifort"
 FLAGS = -static -debug -warn all -check all -FR -O0 -auto -WB -traceback -g -fltconsistency -fpe0
endif
ifeq "$(FC)" "pgf90"
 FLAGS_PROF = -Bstatic -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv -Mprof=lines,time
 FLAGS = -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv -rpath,/opt/netcdf4-pgi/lib
endif

#========================================================================
# PART 3: Checks
#======================================================================
# make sure that the paths are defined. These are just some high level checks
ifndef F_MASTER
 $(error F_MASTER is undefined)
endif
ifndef FC
 $(error FC is undefined: Specify your compiler)
endif
ifndef FLAGS
 $(error Specify flags for your compiler: $(FC))
endif
ifndef NCDF_PATH
 $(error NCDF_PATH is undefined)
endif

#========================================================================
# PART 4: Compile the puppy
#========================================================================

#FLAGS = -O3 -W -v

#.SUFFIXES: .f .o .f90

# Compile
all: compile_scrf link clean install

check:
	echo test
	echo $(FC)

# compile SCRF code
compile_scrf:
	$(FC) $(FLAGS) -c $(SCRF_ALL) $(DRIVER) \
	$(INCNETCDF)

# link routines
link:
#	$(FC) -fPIC -Bstatic_pgi -rpath,/usr/local/netcdf4-pgi/lib *.o -I./ $(LIBNETCDF) -o $(DRIVER__EX)
#	$(FC) -fPIC -rpath,/usr/local/netcdf4-pgi/lib -I./ $(LIBNETCDF) -o $(DRIVER__EX) *.o 
	$(FC) -fPIC -rpath,/opt/netcdf4-pgi/lib -I./ $(LIBNETCDF) -o $(DRIVER__EX) *.o 

# Remove object files
clean:
						rm -f *.o
						rm -f *.mod
#
# Copy the executable to the bin directory
install:
		mv $(DRIVER__EX) $(EXE_PATH)



