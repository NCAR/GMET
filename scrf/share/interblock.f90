MODULE interblock
 IMPLICIT NONE
 ! ---------------------------------------------------------------------------------------
 ! DEFINES EXPLICIT INTERFACES BETWEEN SUB-PROGRAMS
 ! ---------------------------------------------------------------------------------------
 ! interface for get_2dgrid -- used to get a generic 2-d grid
 INTERFACE
  SUBROUTINE GET_2DGRID(TSERIES_FILE,STANDARDNAME,CELL_METHODS,GEN_2DG,NENSM,NSPL1,NSPL2,NTIME,&
                        IDAT0,IDAT1,ISIM0,ISIM1,UNITSTR,VEXIST,TEXIST,EEXIST)
   USE nrtype                                                   ! variable types etc.
   USE dat_2dgrid                                               ! generic 2-d grid
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   CHARACTER(LEN=120),INTENT(IN)               :: STANDARDNAME  ! standard_name of data variable
   CHARACTER(LEN=120),INTENT(IN)               :: CELL_METHODS  ! cell methods for data variable
   TYPE(GENDAT),POINTER                        :: GEN_2DG       ! 2-d data structure
   INTEGER(I4B),INTENT(IN)                     :: NENSM         ! # ensemble members
   INTEGER(I4B),INTENT(OUT)                    :: NSPL1         ! # points 1st spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: NSPL2         ! # points 2nd spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: NTIME         ! number of time steps
   INTEGER(I4B),INTENT(OUT)                    :: IDAT0         ! first index in data array (time vector)
   INTEGER(I4B),INTENT(OUT)                    :: IDAT1         ! last index in data array (time vector)
   INTEGER(I4B),INTENT(OUT)                    :: ISIM0         ! first index of simulation array
   INTEGER(I4B),INTENT(OUT)                    :: ISIM1         ! last index of simulation array
   CHARACTER(LEN=50),INTENT(OUT)               :: UNITSTR       ! units of data variable
   LOGICAL(LGT),INTENT(OUT)                    :: VEXIST        ! flag if the variable exists
   LOGICAL(LGT),INTENT(OUT)                    :: TEXIST        ! flag if variable depends on time
   LOGICAL(LGT),INTENT(OUT)                    :: EEXIST        ! flag if ensembles exist
  END SUBROUTINE GET_2DGRID
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_station -- used to get lat-lon from station data
 INTERFACE
  SUBROUTINE XY_STATION(TSERIES_FILE,YVAR_ID,XVAR_ID,XY_LAT,XY_LON,&
                        SPL1_START,SPL1_COUNT,SPL2_START,SPL2_COUNT,VAR_ID)
   USE nrtype                                                   ! variable types etc.
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   INTEGER(I4B),INTENT(IN)                     :: YVAR_ID       ! NetCDF ID for latitude
   INTEGER(I4B),INTENT(IN)                     :: XVAR_ID       ! NetCDF ID for longitude
   INTEGER(I4B),INTENT(IN)                     :: VAR_ID        ! NetCDF ID for variable
   REAL(DP),DIMENSION(:,:),POINTER             :: XY_LAT        ! 2-dimensional grid of latitudes
   REAL(DP),DIMENSION(:,:),POINTER             :: XY_LON        ! 2-dimensional grid of longitudes
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL1_START    ! start index for 1st spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: SPL1_COUNT    ! count for 1st spatial dimension
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL2_START    ! start index for 2nd spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: SPL2_COUNT    ! count for 2nd spatial dimension
  END SUBROUTINE XY_STATION
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_lat_lon -- used to get lat-lon from a regular grid
 INTERFACE
  SUBROUTINE XY_LAT_LON(TSERIES_FILE,YVAR_ID,XVAR_ID,XY_LAT,XY_LON,&
                        SPL1_START,SPL1_COUNT,SPL2_START,SPL2_COUNT)
   USE nrtype                                                   ! variable types etc.
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   INTEGER(I4B),INTENT(IN)                     :: YVAR_ID       ! NetCDF ID for latitude
   INTEGER(I4B),INTENT(IN)                     :: XVAR_ID       ! NetCDF ID for longitude
   REAL(DP),DIMENSION(:,:),POINTER             :: XY_LAT        ! 2-dimensional grid of latitudes
   REAL(DP),DIMENSION(:,:),POINTER             :: XY_LON        ! 2-dimensional grid of longitudes
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL1_START    ! start index for 1st spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: SPL1_COUNT    ! count for 1st spatial dimension
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL2_START    ! start index for 2nd spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: SPL2_COUNT    ! count for 2nd spatial dimension
  END SUBROUTINE XY_LAT_LON
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_rotated -- used to get lat-lon from a rotated grid
 INTERFACE
  SUBROUTINE XY_ROTATED(TSERIES_FILE,YVAR_ID,XVAR_ID,XY_LAT,XY_LON,&
                        SPL1_START,SPL1_COUNT,SPL2_START,SPL2_COUNT)
   USE nrtype                                                   ! variable types etc.
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   INTEGER(I4B),INTENT(IN)                     :: YVAR_ID       ! NetCDF ID for latitude
   INTEGER(I4B),INTENT(IN)                     :: XVAR_ID       ! NetCDF ID for longitude
   REAL(DP),DIMENSION(:,:),POINTER             :: XY_LAT        ! 2-dimensional grid of latitudes
   REAL(DP),DIMENSION(:,:),POINTER             :: XY_LON        ! 2-dimensional grid of longitudes
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL1_START    ! start index for 1st spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: SPL1_COUNT    ! count for 1st spatial dimension
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL2_START    ! start index for 2nd spatial dimension
   INTEGER(I4B),INTENT(OUT)                    :: SPL2_COUNT    ! count for 2nd spatial dimension
  END SUBROUTINE XY_ROTATED
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for get_elevtn -- used to get elevation data
 INTERFACE
  SUBROUTINE GET_ELEVTN(TSERIES_FILE,CSYSTEM,XY_ELV,&
                        SPL1_START,SPL1_COUNT,SPL2_START,SPL2_COUNT)
   USE nrtype                                                   ! variable types etc.
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   CHARACTER(LEN=12),INTENT(IN)                :: CSYSTEM       ! Coordinate system of raw data
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL1_START    ! start index for 1st spatial dimension
   INTEGER(I4B),INTENT(IN)                     :: SPL1_COUNT    ! count for 1st spatial dimension
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL2_START    ! start index for 2nd spatial dimension
   INTEGER(I4B),INTENT(IN)                     :: SPL2_COUNT    ! count for 2nd spatial dimension
   REAL(DP),DIMENSION(:,:),POINTER             :: XY_ELV        ! 2-dimensional grid of elevation
  END SUBROUTINE GET_ELEVTN
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for get_timdat -- used to get full time arrays
 INTERFACE
  SUBROUTINE GET_TIMDAT(TSERIES_FILE,TVAR_ID,TB0,TB1,TIM,BEX,NTIM_DATA)
   USE nrtype                                                   ! variable types etc.
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   INTEGER(I4B),INTENT(IN)                     :: TVAR_ID       ! NetCDF ID for time
   REAL(DP),DIMENSION(:),POINTER               :: TB0           ! start of time interval (in data file)
   REAL(DP),DIMENSION(:),POINTER               :: TB1           ! end of time interval (in data file)
   REAL(DP),DIMENSION(:),POINTER               :: TIM           ! time stamp (in data file)
   LOGICAL(LGT),INTENT(OUT)                    :: BEX           ! .TRUE. if time bounds exist
   INTEGER(I4B),INTENT(OUT)                    :: NTIM_DATA     ! number of time intervals in the data file
  END SUBROUTINE GET_TIMDAT
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for get_subtim -- used to get time subset
 INTERFACE
  SUBROUTINE GET_SUBTIM(TB0,TB1,TIM,BEX,NTIM,TIME_START,TIME_COUNT,ISIM0,ISIM1)
   USE nrtype                                                   ! variable types etc.
   REAL(DP),DIMENSION(:),POINTER               :: TB0           ! start of time interval (in data file)
   REAL(DP),DIMENSION(:),POINTER               :: TB1           ! end of time interval (in data file)
   REAL(DP),DIMENSION(:),POINTER               :: TIM           ! time stamp (in data file)
   LOGICAL(LGT),INTENT(IN)                     :: BEX           ! .TRUE. if time bounds exist
   INTEGER(I4B),INTENT(IN)                     :: NTIM          ! number of time intervals in data
   INTEGER(I4B),INTENT(OUT)                    :: TIME_START    ! start index for time
   INTEGER(I4B),INTENT(OUT)                    :: TIME_COUNT    ! count for time
   INTEGER(I4B),INTENT(OUT)                    :: ISIM0         ! first index of simulation array
   INTEGER(I4B),INTENT(OUT)                    :: ISIM1         ! last index of simulation array
  END SUBROUTINE GET_SUBTIM
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for get_subset -- used to get data subset
 INTERFACE
  SUBROUTINE GET_SUBSET(TSERIES_FILE,CSYSTEM,VAR_ID,TDIM_IX,ZDIM_IX,YDIM_IX,XDIM_IX,EDIM_IX,&
                        TIME_START,TIME_COUNT,SPL1_START,SPL1_COUNT,SPL2_START,SPL2_COUNT,&
                        VERT_START,VERT_COUNT,IENS_START,IENS_COUNT,&
                        RAW_DAT,ENS_DAT,UNITSTR)
   USE nrtype                                                   ! variable types etc.
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   CHARACTER(LEN=12),INTENT(IN)                :: CSYSTEM       ! Coordinate system of raw data
   INTEGER(I4B),INTENT(IN)                     :: VAR_ID        ! NetCDF ID for data variable
   INTEGER(I4B),INTENT(IN)                     :: TDIM_IX       ! NetCDF dimension index for time variable
   INTEGER(I4B),INTENT(IN)                     :: ZDIM_IX       ! NetCDF dimension index for (height, depth)
   INTEGER(I4B),INTENT(IN)                     :: YDIM_IX       ! NetCDF dimension index for latitude
   INTEGER(I4B),INTENT(IN)                     :: XDIM_IX       ! NetCDF dimension index for longitude
   INTEGER(I4B),INTENT(IN)                     :: EDIM_IX       ! NetCDF dimension index for ensemble
   INTEGER(I4B),INTENT(IN)                     :: TIME_START    ! start index for time
   INTEGER(I4B),INTENT(IN)                     :: TIME_COUNT    ! count for time
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL1_START    ! start index for 1st spatial dimension
   INTEGER(I4B),INTENT(IN)                     :: SPL1_COUNT    ! count for 1st spatial dimension
   INTEGER(I4B),DIMENSION(:),POINTER           :: SPL2_START    ! start index for 2nd spatial dimension
   INTEGER(I4B),INTENT(IN)                     :: SPL2_COUNT    ! count for 2nd spatial dimension
   INTEGER(I4B),INTENT(IN)                     :: VERT_START    ! start index for vertical dimension
   INTEGER(I4B),INTENT(IN)                     :: VERT_COUNT    ! count for vertical dimension
   INTEGER(I4B),INTENT(IN)                     :: IENS_START    ! start index for ensemble dimension
   INTEGER(I4B),INTENT(IN)                     :: IENS_COUNT    ! count for ensemble dimension
   REAL(DP),DIMENSION(:,:,:),POINTER           :: RAW_DAT       ! raw data array (NSPL1,NSPL2,NTIME)
   REAL(DP),DIMENSION(:,:,:,:),POINTER         :: ENS_DAT       ! ensemble data array (NENSM,NSPL1,NSPL2,NTIME)
   CHARACTER(LEN=50),INTENT(OUT)               :: UNITSTR       ! units of data variable
  END SUBROUTINE GET_SUBSET
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for grid2basin -- used to get relationship between grids and basins
 INTERFACE
  SUBROUTINE GRID2BASIN(TSERIES_FILE,STANDARDNAME,CELL_METHODS,NSPL1,NSPL2,NRCH,DTYPE,GEN_2DG,GEN_SRF)
   USE nrtype                                                   ! variable types (DP, I4B, etc.)
   USE dat_2dgrid                                               ! generic 2-d data structures
   CHARACTER(LEN=120),INTENT(IN)               :: TSERIES_FILE  ! name of time series file
   CHARACTER(LEN=120),INTENT(IN)               :: STANDARDNAME  ! standard_name of data variable
   CHARACTER(LEN=120),INTENT(IN)               :: CELL_METHODS  ! cell methods for data variable
   INTEGER(I4B), INTENT(IN)                    :: NSPL1         ! number of x points
   INTEGER(I4B), INTENT(IN)                    :: NSPL2         ! number of y points
   INTEGER(I4B), INTENT(IN)                    :: NRCH          ! number of sub-basins
   INTEGER(I4B), INTENT(IN)                    :: DTYPE         ! input data (1-7: p,t,rh,srad,lrad,wind,pressuse)
   TYPE(GENDAT),POINTER                        :: GEN_2DG       ! data structure for generic 2-d grid
   TYPE(GENDAT),POINTER,OPTIONAL               :: GEN_SRF       ! generic 2-d grid -- mean surface
  END SUBROUTINE GRID2BASIN
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for grid2basin -- used to get relationship between grids and basins
 INTERFACE
  SUBROUTINE BASLINKAGE(NSPL1,NSPL2,NRCH,DTYPE,GEN_SRF)
   USE nrtype
   USE dat_2dgrid                                               ! generic 2-d data structures
   INTEGER(I4B), INTENT(IN)                    :: NSPL1         ! number of x points
   INTEGER(I4B), INTENT(IN)                    :: NSPL2         ! number of y points
   INTEGER(I4B), INTENT(IN)                    :: NRCH          ! number of sub-basins
   INTEGER(I4B), INTENT(IN)                    :: DTYPE         ! input data (1-7: p,t,rh,srad,lrad,wind,pressure)
   TYPE(GENDAT),POINTER,OPTIONAL               :: GEN_SRF       ! generic 2-d grid -- mean surface
  END SUBROUTINE BASLINKAGE
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for satfrc_bas -- used to get saturated areas
 INTERFACE
  SUBROUTINE SATFRC_BAS(IENS,IBAS,NBINS,UAREA,IAREA,IATNB,SAREA)
   USE nrtype
   INTEGER(I4B), INTENT(IN)                  :: IENS    ! Ensemble memeber
   INTEGER(I4B), INTENT(IN)                  :: IBAS    ! Catchment ID
   INTEGER(I4B), INTENT(IN)                  :: NBINS   ! Number of a/tan(b) classes
   REAL(DP), INTENT(OUT)                     :: UAREA   ! catchment area UNINFLUENCED
   REAL(DP), DIMENSION(NBINS),INTENT(OUT)    :: IAREA   ! catchment area INFLUENCED
   REAL(DP), DIMENSION(NBINS),INTENT(OUT)    :: IATNB   ! representative a/tan(b) values 
   REAL(DP), INTENT(OUT)                     :: SAREA   ! catchment area SATURATED
  END SUBROUTINE SATFRC_BAS
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for getusq_rch -- used to retrieve routed flow particles from u/s reaches and
 !                             non-routed flow particles from the current reach
 INTERFACE
  SUBROUTINE GETUSQ_RCH(IENS,JRCH,Q_JRCH,TENTRY,T_EXIT,RSTEP)
  USE nrtype
  INTEGER(I4B), INTENT(IN)                   :: IENS    ! ensemble member
  INTEGER(I4B), INTENT(IN)                   :: JRCH    ! reach to process
  REAL(DP), DIMENSION(:), POINTER            :: Q_JRCH  ! merged (non-routed) flow in JRCH
  REAL(DP), DIMENSION(:), POINTER            :: TENTRY  ! time flow particles entered JRCH
  REAL(DP), DIMENSION(:), POINTER            :: T_EXIT  ! time flow expected to exit JRCH
  INTEGER(I4B), INTENT(IN), OPTIONAL         :: RSTEP   ! retrospective time step offset
  END SUBROUTINE GETUSQ_RCH
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for qexmul_rch -- used to extract flow from multiple u/s reaches
 INTERFACE
  SUBROUTINE QEXMUL_RCH(IENS,JRCH,ND,QD,TD,RSTEP)
  USE nrtype
  INTEGER(I4B), INTENT(IN)                   :: IENS    ! ensemble member
  INTEGER(I4B), INTENT(IN)                   :: JRCH    ! reach to process
  INTEGER(I4B), INTENT(OUT)                  :: ND      ! number of routed particles
  REAL(DP), DIMENSION(:), POINTER            :: QD      ! flow particles just enetered JRCH
  REAL(DP), DIMENSION(:), POINTER            :: TD      ! time flow particles entered JRCH
  INTEGER(I4B), INTENT(IN), OPTIONAL         :: RSTEP   ! retrospective time step offset
  END SUBROUTINE QEXMUL_RCH
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for extract_from_reach -- used to extract flow from multiple u/s reaches
 INTERFACE
   SUBROUTINE EXTRACT_FROM_RCH(IENS,JRCH,NR,Q_JRCH,T_EXIT,T_END,TNEW)
   USE nrtype
   INTEGER(I4B), INTENT(IN)                    :: IENS   ! ensemble member
   INTEGER(I4B), INTENT(IN)                    :: JRCH   ! reach to process
   INTEGER(I4B), INTENT(IN)                    :: NR     ! number of routed particles
   REAL(DP),DIMENSION(:),POINTER               :: Q_JRCH ! flow in downstream reach JRCH
   REAL(DP),DIMENSION(:),POINTER               :: T_EXIT ! time particle expected exit JRCH
   REAL(DP)                                    :: T_END  ! end of time step
   REAL(DP),DIMENSION(2)                       :: TNEW   ! start/end of time step
   END SUBROUTINE EXTRACT_FROM_RCH
  END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for remove_rch -- used to extract flow from multiple u/s reaches
 INTERFACE
  SUBROUTINE REMOVE_RCH(Q_JRCH,TENTRY,T_EXIT)
  USE nrtype
  REAL(DP), DIMENSION(:), POINTER            :: Q_JRCH  ! merged (non-routed) flow in JRCH
  REAL(DP), DIMENSION(:), POINTER            :: TENTRY  ! time flow particles entered JRCH
  REAL(DP), DIMENSION(:), POINTER            :: T_EXIT  ! time flow particles exited JRCH
  END SUBROUTINE REMOVE_RCH
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for kinwav_rch -- used to propogate flow particles thru single river segment
 INTERFACE
  SUBROUTINE KINWAV_RCH(JRCH,T_START,T_END,Q_JRCH,TENTRY,&    ! input
                                           FROUTE,T_EXIT,NQ2) ! output
  USE nrtype
  INTEGER(I4B), INTENT(IN)                    :: JRCH   ! Reach to process
  REAL(DP), INTENT(IN)                        :: T_START! start of the time step
  REAL(DP), INTENT(IN)                        :: T_END  ! end of the time step
  REAL(DP), DIMENSION(:), INTENT(INOUT)       :: Q_JRCH ! flow to be routed
  REAL(DP), DIMENSION(:), INTENT(INOUT)       :: TENTRY ! time to be routed
  LOGICAL(LGT), DIMENSION(:), INTENT(INOUT)   :: FROUTE ! routing flag, T=routed
  REAL(DP), DIMENSION(:), INTENT(INOUT)       :: T_EXIT ! time pts expected exit segment
  INTEGER(I4B), INTENT(OUT)                   :: NQ2    ! # particles (<= input b/c merge)
  END SUBROUTINE KINWAV_RCH
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for interp_rch -- used to compute time step average flow
 INTERFACE
  SUBROUTINE INTERP_RCH(TOLD,QOLD,TNEW,QNEW,IERR)
  USE nrtype
  REAL(DP), DIMENSION(:), INTENT(IN)          :: TOLD   ! input time array
  REAL(DP), DIMENSION(:), INTENT(IN)          :: QOLD   ! input flow array
  REAL(DP), DIMENSION(:), INTENT(IN)          :: TNEW   ! desired output times
  REAL(DP), DIMENSION(:), INTENT(OUT)         :: QNEW   ! flow averaged for desired times
  INTEGER(I4B), INTENT(OUT)                   :: IERR   ! error, 1= bad bounds
  END SUBROUTINE INTERP_RCH
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for linearcorr -- used to compute the linear correlation coefficient
 INTERFACE
  SUBROUTINE LINEARCORR(X,Y,R)
  USE nrtype
  REAL(DP), INTENT(OUT)                       :: R      ! correlation coefficient
  REAL(DP), DIMENSION(:), INTENT(IN)          :: X,Y    ! input vectors
  END SUBROUTINE LINEARCORR
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for sqrtmatrix -- used to compute the square root of a matrix
 INTERFACE
  SUBROUTINE SQRTMATRIX(A,B,N)
  USE nrtype
  REAL(DP), DIMENSION(:,:), INTENT(IN)        :: A      ! input matrix
  REAL(DP), DIMENSION(:,:), INTENT(OUT)       :: B      ! square root of A (A=B*B)
  INTEGER(I4B), INTENT(IN)                    :: N      ! number of elements of A
  END SUBROUTINE SQRTMATRIX
 END INTERFACE
 ! ---------------------------------------------------------------------------------------
 ! interface for inv_matrix -- used to invert a matrix
 INTERFACE
  SUBROUTINE INV_MATRIX(A,B,N)
  USE nrtype
  REAL(DP), DIMENSION(:,:), INTENT(IN)        :: A      ! input matrix
  REAL(DP), DIMENSION(:,:), INTENT(OUT)       :: B      ! inverse of A (B = A**-1)
  INTEGER(I4B), INTENT(IN)                    :: N      ! number of elements of A
  END SUBROUTINE INV_MATRIX
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for raindisagg
 INTERFACE
  SUBROUTINE RAINDISAGG(NUM_DIS,RAIN,NLEVELS,TIMINT)
  USE nrtype
  INTEGER(I4B), INTENT(IN)                    :: NUM_DIS       ! Number of disaggregated elements
  REAL(DP), DIMENSION(NUM_DIS),INTENT(INOUT)  :: RAIN          ! Rainfall array
  INTEGER(I4B), INTENT(IN)                    :: NLEVELS       ! Number of the levels in the cascade
  REAL(DP), INTENT(IN)                        :: TIMINT        ! time interval of data (seconds)
  END SUBROUTINE RAINDISAGG
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for gridsubset
 INTERFACE ! required because of the pointer arguments in the subroutine
  SUBROUTINE GRIDSUBSET(LON2D_A,LAT2D_A,LON2D_B,LAT2D_B,THRSHLD,OFFSET,METHOD,MASK,MISSING,GRIDLIST)
    USE nrtype                                                 ! variable types (DP, I4B, etc.)
    USE grid_lists                                             ! structure for grid selection lists
    REAL(DP),DIMENSION(:,:),POINTER           :: LON2D_A       ! 2d array of lon for basic grid
    REAL(DP),DIMENSION(:,:),POINTER           :: LAT2D_A       ! 2d array of lat for basic grid
    REAL(DP),DIMENSION(:,:),POINTER           :: LON2D_B       ! 2d array of lon for input grid
    REAL(DP),DIMENSION(:,:),POINTER           :: LAT2D_B       ! 2d array of lat for input grid
    REAL(DP),INTENT(IN)                       :: THRSHLD       ! threshold distance
    REAL(DP),INTENT(IN)                       :: OFFSET        ! offset from initial search area
    INTEGER(I4B),INTENT(IN)                   :: METHOD        ! method used to select input grid points
    LOGICAL(LGT),DIMENSION(:,:),POINTER       :: MASK          ! input data grid points to consider
    LOGICAL(LGT),INTENT(IN)                   :: MISSING       ! if true missing points allowed
    TYPE(GRDARRAY),DIMENSION(:,:),POINTER     :: GRIDLIST      ! structure for output
  END SUBROUTINE GRIDSUBSET
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for model1step
 INTERFACE ! required because of the optional argument
  SUBROUTINE MODEL1STEP(NENS,NRCH,MRCH,NHYD,NLAK,NOBL,IOFFSET,ISTEP,IFORCE0,IFORCE1,RSTEP)
    USE nrtype
    INTEGER(I4B), INTENT(IN)               :: NENS        ! number of ensemble members
    INTEGER(I4B), INTENT(IN)               :: NRCH        ! number of basins
    INTEGER(I4B), INTENT(IN)               :: MRCH        ! number of basins selected for output
    INTEGER(I4B), INTENT(IN)               :: NHYD        ! number of water level recorders
    INTEGER(I4B), INTENT(IN)               :: NOBL        ! number of lake level observation locations
    INTEGER(I4B), INTENT(IN)               :: NLAK        ! number of lakes
    INTEGER(I4B), INTENT(IN)               :: IOFFSET     ! timestep offset from start of record
    INTEGER(I4B), INTENT(IN)               :: ISTEP       ! loop through time steps
    INTEGER(I4B), INTENT(IN)               :: IFORCE0     ! starting index for data
    INTEGER(I4B), INTENT(IN)               :: IFORCE1     ! end index for data
    INTEGER(I4B), INTENT(IN), OPTIONAL     :: RSTEP       ! retrospective time step offset
  END SUBROUTINE MODEL1STEP
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for q_rch_lake
 INTERFACE ! required because of the optional argument
  SUBROUTINE Q_RCH_LAKE(IENS,JRCH,NOBL,NHYD,ISTEP,RSTEP)
    USE nrtype
    INTEGER(I4B), INTENT(IN)               :: IENS        ! loop through ensemble members
    INTEGER(I4B), INTENT(IN)               :: JRCH        ! loop through the stream segments
    INTEGER(I4B), INTENT(IN)               :: NOBL        ! number of lake level observation locations
    INTEGER(I4B), INTENT(IN)               :: NHYD        ! number of streamq observation locations
    INTEGER(I4B), INTENT(IN)               :: ISTEP       ! loop through time steps
    INTEGER(I4B), INTENT(IN), OPTIONAL     :: RSTEP       ! retrospective time step offset
  END SUBROUTINE Q_RCH_LAKE
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for qroute_rch
 INTERFACE ! required because of the optional argument
  SUBROUTINE QROUTE_RCH(IENS,JRCH,RSTEP)
    USE nrtype
    INTEGER(I4B), INTENT(IN)                    :: IENS   ! ensemble member
    INTEGER(I4B), INTENT(IN)                    :: JRCH   ! reach to process
    INTEGER(I4B), INTENT(IN), OPTIONAL          :: RSTEP  ! retrospective time step offset
  END SUBROUTINE QROUTE_RCH
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface to write_ts_sp
 INTERFACE ! required because of the optional argument
  SUBROUTINE WRITE_TS_SP(ISTEP,NENS,MRCH,NRCH,NHYD,NPCP,NLAK,NOBL,SECTION,RSTEP)
    USE nrtype
    ! input variables
    INTEGER(I4B),INTENT(IN)                        :: ISTEP  ! loop through time steps
    INTEGER(I4B),INTENT(IN)                        :: NENS   ! number of ensemble members
    INTEGER(I4B),INTENT(IN)                        :: MRCH   ! number of basins selected for output
    INTEGER(I4B),INTENT(IN)                        :: NRCH   ! number of basins
    INTEGER(I4B),INTENT(IN)                        :: NHYD   ! number of water level recorders
    INTEGER(I4B),INTENT(IN)                        :: NPCP   ! number of precipitation stations
    INTEGER(I4B),INTENT(IN)                        :: NLAK   ! number of lakes
    INTEGER(I4B),INTENT(IN)                        :: NOBL   ! number of lake level observation locations
    INTEGER(I4B),INTENT(IN),OPTIONAL               :: SECTION  ! section of data to output
    INTEGER(I4B),INTENT(IN),OPTIONAL               :: RSTEP    ! retrospective time step offset
  END SUBROUTINE WRITE_TS_SP
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface to err_calc
 INTERFACE ! required because of the optional argument
  SUBROUTINE ERR_CALC(ISTEP,IOFFSET,NENS,MRCH,NRCH,NHYD,JRUN,RSTEP)
    USE nrtype
    INTEGER(I4B),INTENT(IN)                        :: ISTEP  ! loop through time steps
    INTEGER(I4B),INTENT(IN)                        :: IOFFSET ! timestep offset from start of record
    INTEGER(I4B),INTENT(IN)                        :: NENS   ! number of ensemble members
    INTEGER(I4B),INTENT(IN)                        :: MRCH   ! number of basins selected for output
    INTEGER(I4B),INTENT(IN)                        :: NRCH   ! number of basins
    INTEGER(I4B),INTENT(IN)                        :: NHYD   ! number of water level recorders
    INTEGER(I4B),INTENT(IN)                        :: JRUN   ! index of the model error parameter set used
    INTEGER(I4B),INTENT(IN),OPTIONAL               :: RSTEP    ! retrospective time step offset
  END SUBROUTINE ERR_CALC
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for sp_out_put
 INTERFACE ! required because of the optional argument
  SUBROUTINE SP_OUT_PUT(ISTEP,ITIM,NENS,NRCH,NHYD,SECTION,RSTEP)
    USE nrtype
    INTEGER(I4B), INTENT(IN)               :: ISTEP         ! Model time step
    INTEGER(I4B), INTENT(IN)               :: ITIM          ! Index of spatial data on snow covered area
    INTEGER(I4B), INTENT(IN)               :: NENS          ! Number of ensemble members
    INTEGER(I4B), INTENT(IN)               :: NRCH          ! Number of basins
    INTEGER(I4B), INTENT(IN)               :: NHYD          ! Number of gages
    INTEGER(I4B), INTENT(IN)               :: SECTION       ! Section of variables to write out
    INTEGER(I4B), INTENT(IN), OPTIONAL     :: RSTEP         ! retrospective time step
  END SUBROUTINE SP_OUT_PUT
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for ts_out_put
 INTERFACE ! required because of the optional argument
  SUBROUTINE TS_OUT_PUT(ISTEP,NENS,MRCH,NRCH,NHYD,NPCP,NLAK,NOBL,SECTION,RSTEP)
    USE nrtype
    INTEGER(I4B), INTENT(IN)               :: ISTEP         ! Model time step
    INTEGER(I4B), INTENT(IN)               :: NENS          ! Number of ensemble members
    INTEGER(I4B), INTENT(IN)               :: MRCH          ! Number of selected basins
    INTEGER(I4B), INTENT(IN)               :: NRCH          ! Number of basins
    INTEGER(I4B), INTENT(IN)               :: NHYD          ! Number of water level recorders
    INTEGER(I4B), INTENT(IN)               :: NOBL          ! Number of lake level observation locations
    INTEGER(I4B), INTENT(IN)               :: NPCP          ! Number of precipitation stations
    INTEGER(I4B), INTENT(IN)               :: NLAK          ! Number of lakes
    INTEGER(I4B), INTENT(IN)               :: SECTION       ! Section of variables to write out
    INTEGER(I4B), INTENT(IN), OPTIONAL     :: RSTEP         ! retrospective time step
  END SUBROUTINE TS_OUT_PUT
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for prberr_put
 INTERFACE ! required because of the optional argument
  SUBROUTINE PRBERR_PUT(IOFFSET,ISTEP,JRUN,NHYD,RSTEP)
    USE nrtype
    INTEGER(I4B), INTENT(IN)               :: IOFFSET       ! Offset for data write statement
    INTEGER(I4B), INTENT(IN)               :: ISTEP         ! Model time step
    INTEGER(I4B), INTENT(IN)               :: JRUN          ! Index of model error parameter set
    INTEGER(I4B), INTENT(IN)               :: NHYD          ! Number of gages
    INTEGER(I4B), INTENT(IN), OPTIONAL     :: RSTEP         ! Retrospective time step
  END SUBROUTINE PRBERR_PUT
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for rstart_get
 INTERFACE ! required because of the optional argument
  SUBROUTINE RSTART_GET(NENS,NRCH,NLAK,DTIME,LEXIST,FEXIST)
    USE nrtype
    INTEGER(I4B), INTENT(IN)               :: NENS      ! Number of ensemble members
    INTEGER(I4B), INTENT(IN)               :: NRCH      ! Number of basins
    INTEGER(I4B), INTENT(IN)               :: NLAK      ! Number of lakes
    REAL(DP), INTENT(IN)                   :: DTIME     ! Time difference between restart file and model reference time
    LOGICAL(LGT), INTENT(OUT)              :: LEXIST    ! .TRUE. if re-start file exists
    LOGICAL(LGT), INTENT(IN), OPTIONAL     :: FEXIST    ! if .TRUE., only check for file existance
  END SUBROUTINE RSTART_GET
 END INTERFACE
 ! ----------------------------------------------------------------------------------------
 ! interface for moddat_get
 INTERFACE ! required because of the optional argument
  SUBROUTINE MODDAT_GET(NENS,NRCH,DTIME,LEXIST,FEXIST)
    USE nrtype
    INTEGER(I4B), INTENT(IN)               :: NENS      ! Number of ensemble members
    INTEGER(I4B), INTENT(IN)               :: NRCH      ! Number of basins
    REAL(DP), INTENT(IN)                   :: DTIME     ! Time difference between data_save file and model reference time
    LOGICAL(LGT), INTENT(OUT)              :: LEXIST    ! .TRUE. if re-start file exists
    LOGICAL(LGT), INTENT(IN), OPTIONAL     :: FEXIST    ! if .TRUE., only check for file existance
  END SUBROUTINE MODDAT_GET
 END INTERFACE
! ----------------------------------------------------------------------------------------
END MODULE interblock
