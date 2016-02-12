Module interblock
  Implicit None
 ! ---------------------------------------------------------------------------------------
 ! DEFINES EXPLICIT INTERFACES BETWEEN SUB-PROGRAMS
 ! ---------------------------------------------------------------------------------------
 ! interface for get_2dgrid -- used to get a generic 2-d grid
  Interface
    Subroutine GET_2DGRID (TSERIES_FILE, STANDARDNAME, CELL_METHODS, GEN_2DG, NENSM, NSPL1, NSPL2, NTIME, IDAT0, IDAT1, &
   & ISIM0, ISIM1, UNITSTR, VEXIST, TEXIST, EEXIST)
      Use nrtype ! variable types etc.
      Use dat_2dgrid ! generic 2-d grid
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Character (Len=120), Intent (In) :: STANDARDNAME ! standard_name of data variable
      Character (Len=120), Intent (In) :: CELL_METHODS ! cell methods for data variable
      Type (GENDAT), Pointer :: GEN_2DG ! 2-d data structure
      Integer (I4B), Intent (In) :: NENSM ! # ensemble members
      Integer (I4B), Intent (Out) :: NSPL1 ! # points 1st spatial dimension
      Integer (I4B), Intent (Out) :: NSPL2 ! # points 2nd spatial dimension
      Integer (I4B), Intent (Out) :: NTIME ! number of time steps
      Integer (I4B), Intent (Out) :: IDAT0 ! first index in data array (time vector)
      Integer (I4B), Intent (Out) :: IDAT1 ! last index in data array (time vector)
      Integer (I4B), Intent (Out) :: ISIM0 ! first index of simulation array
      Integer (I4B), Intent (Out) :: ISIM1 ! last index of simulation array
      Character (Len=50), Intent (Out) :: UNITSTR ! units of data variable
      Logical (LGT), Intent (Out) :: VEXIST ! flag if the variable exists
      Logical (LGT), Intent (Out) :: TEXIST ! flag if variable depends on time
      Logical (LGT), Intent (Out) :: EEXIST ! flag if ensembles exist
    End Subroutine GET_2DGRID
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_station -- used to get lat-lon from station data
  Interface
    Subroutine XY_STATION (TSERIES_FILE, YVAR_ID, XVAR_ID, XY_LAT, XY_LON, SPL1_START, SPL1_COUNT, SPL2_START, &
   & SPL2_COUNT, VAR_ID)
      Use nrtype ! variable types etc.
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Integer (I4B), Intent (In) :: YVAR_ID ! NetCDF ID for latitude
      Integer (I4B), Intent (In) :: XVAR_ID ! NetCDF ID for longitude
      Integer (I4B), Intent (In) :: VAR_ID ! NetCDF ID for variable
      Real (DP), Dimension (:, :), Pointer :: XY_LAT ! 2-dimensional grid of latitudes
      Real (DP), Dimension (:, :), Pointer :: XY_LON ! 2-dimensional grid of longitudes
      Integer (I4B), Dimension (:), Pointer :: SPL1_START ! start index for 1st spatial dimension
      Integer (I4B), Intent (Out) :: SPL1_COUNT ! count for 1st spatial dimension
      Integer (I4B), Dimension (:), Pointer :: SPL2_START ! start index for 2nd spatial dimension
      Integer (I4B), Intent (Out) :: SPL2_COUNT ! count for 2nd spatial dimension
    End Subroutine XY_STATION
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_lat_lon -- used to get lat-lon from a regular grid
  Interface
    Subroutine XY_LAT_LON (TSERIES_FILE, YVAR_ID, XVAR_ID, XY_LAT, XY_LON, SPL1_START, SPL1_COUNT, SPL2_START, &
   & SPL2_COUNT)
      Use nrtype ! variable types etc.
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Integer (I4B), Intent (In) :: YVAR_ID ! NetCDF ID for latitude
      Integer (I4B), Intent (In) :: XVAR_ID ! NetCDF ID for longitude
      Real (DP), Dimension (:, :), Pointer :: XY_LAT ! 2-dimensional grid of latitudes
      Real (DP), Dimension (:, :), Pointer :: XY_LON ! 2-dimensional grid of longitudes
      Integer (I4B), Dimension (:), Pointer :: SPL1_START ! start index for 1st spatial dimension
      Integer (I4B), Intent (Out) :: SPL1_COUNT ! count for 1st spatial dimension
      Integer (I4B), Dimension (:), Pointer :: SPL2_START ! start index for 2nd spatial dimension
      Integer (I4B), Intent (Out) :: SPL2_COUNT ! count for 2nd spatial dimension
    End Subroutine XY_LAT_LON
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_rotated -- used to get lat-lon from a rotated grid
  Interface
    Subroutine XY_ROTATED (TSERIES_FILE, YVAR_ID, XVAR_ID, XY_LAT, XY_LON, SPL1_START, SPL1_COUNT, SPL2_START, &
   & SPL2_COUNT)
      Use nrtype ! variable types etc.
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Integer (I4B), Intent (In) :: YVAR_ID ! NetCDF ID for latitude
      Integer (I4B), Intent (In) :: XVAR_ID ! NetCDF ID for longitude
      Real (DP), Dimension (:, :), Pointer :: XY_LAT ! 2-dimensional grid of latitudes
      Real (DP), Dimension (:, :), Pointer :: XY_LON ! 2-dimensional grid of longitudes
      Integer (I4B), Dimension (:), Pointer :: SPL1_START ! start index for 1st spatial dimension
      Integer (I4B), Intent (Out) :: SPL1_COUNT ! count for 1st spatial dimension
      Integer (I4B), Dimension (:), Pointer :: SPL2_START ! start index for 2nd spatial dimension
      Integer (I4B), Intent (Out) :: SPL2_COUNT ! count for 2nd spatial dimension
    End Subroutine XY_ROTATED
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_elevtn -- used to get elevation data
  Interface
    Subroutine GET_ELEVTN (TSERIES_FILE, CSYSTEM, XY_ELV, SPL1_START, SPL1_COUNT, SPL2_START, SPL2_COUNT)
      Use nrtype ! variable types etc.
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Character (Len=12), Intent (In) :: CSYSTEM ! Coordinate system of raw data
      Integer (I4B), Dimension (:), Pointer :: SPL1_START ! start index for 1st spatial dimension
      Integer (I4B), Intent (In) :: SPL1_COUNT ! count for 1st spatial dimension
      Integer (I4B), Dimension (:), Pointer :: SPL2_START ! start index for 2nd spatial dimension
      Integer (I4B), Intent (In) :: SPL2_COUNT ! count for 2nd spatial dimension
      Real (DP), Dimension (:, :), Pointer :: XY_ELV ! 2-dimensional grid of elevation
    End Subroutine GET_ELEVTN
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_timdat -- used to get full time arrays
  Interface
    Subroutine GET_TIMDAT (TSERIES_FILE, TVAR_ID, TB0, TB1, TIM, BEX, NTIM_DATA)
      Use nrtype ! variable types etc.
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Integer (I4B), Intent (In) :: TVAR_ID ! NetCDF ID for time
      Real (DP), Dimension (:), Pointer :: TB0 ! start of time interval (in data file)
      Real (DP), Dimension (:), Pointer :: TB1 ! end of time interval (in data file)
      Real (DP), Dimension (:), Pointer :: TIM ! time stamp (in data file)
      Logical (LGT), Intent (Out) :: BEX ! .TRUE. if time bounds exist
      Integer (I4B), Intent (Out) :: NTIM_DATA ! number of time intervals in the data file
    End Subroutine GET_TIMDAT
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_subtim -- used to get time subset
  Interface
    Subroutine GET_SUBTIM (TB0, TB1, TIM, BEX, NTIM, TIME_START, TIME_COUNT, ISIM0, ISIM1)
      Use nrtype ! variable types etc.
      Real (DP), Dimension (:), Pointer :: TB0 ! start of time interval (in data file)
      Real (DP), Dimension (:), Pointer :: TB1 ! end of time interval (in data file)
      Real (DP), Dimension (:), Pointer :: TIM ! time stamp (in data file)
      Logical (LGT), Intent (In) :: BEX ! .TRUE. if time bounds exist
      Integer (I4B), Intent (In) :: NTIM ! number of time intervals in data
      Integer (I4B), Intent (Out) :: TIME_START ! start index for time
      Integer (I4B), Intent (Out) :: TIME_COUNT ! count for time
      Integer (I4B), Intent (Out) :: ISIM0 ! first index of simulation array
      Integer (I4B), Intent (Out) :: ISIM1 ! last index of simulation array
    End Subroutine GET_SUBTIM
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_subset -- used to get data subset
  Interface
    Subroutine GET_SUBSET (TSERIES_FILE, CSYSTEM, VAR_ID, TDIM_IX, ZDIM_IX, YDIM_IX, XDIM_IX, EDIM_IX, TIME_START, &
   & TIME_COUNT, SPL1_START, SPL1_COUNT, SPL2_START, SPL2_COUNT, VERT_START, VERT_COUNT, IENS_START, IENS_COUNT, &
   & RAW_DAT, ENS_DAT, UNITSTR)
      Use nrtype ! variable types etc.
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Character (Len=12), Intent (In) :: CSYSTEM ! Coordinate system of raw data
      Integer (I4B), Intent (In) :: VAR_ID ! NetCDF ID for data variable
      Integer (I4B), Intent (In) :: TDIM_IX ! NetCDF dimension index for time variable
      Integer (I4B), Intent (In) :: ZDIM_IX ! NetCDF dimension index for (height, depth)
      Integer (I4B), Intent (In) :: YDIM_IX ! NetCDF dimension index for latitude
      Integer (I4B), Intent (In) :: XDIM_IX ! NetCDF dimension index for longitude
      Integer (I4B), Intent (In) :: EDIM_IX ! NetCDF dimension index for ensemble
      Integer (I4B), Intent (In) :: TIME_START ! start index for time
      Integer (I4B), Intent (In) :: TIME_COUNT ! count for time
      Integer (I4B), Dimension (:), Pointer :: SPL1_START ! start index for 1st spatial dimension
      Integer (I4B), Intent (In) :: SPL1_COUNT ! count for 1st spatial dimension
      Integer (I4B), Dimension (:), Pointer :: SPL2_START ! start index for 2nd spatial dimension
      Integer (I4B), Intent (In) :: SPL2_COUNT ! count for 2nd spatial dimension
      Integer (I4B), Intent (In) :: VERT_START ! start index for vertical dimension
      Integer (I4B), Intent (In) :: VERT_COUNT ! count for vertical dimension
      Integer (I4B), Intent (In) :: IENS_START ! start index for ensemble dimension
      Integer (I4B), Intent (In) :: IENS_COUNT ! count for ensemble dimension
      Real (DP), Dimension (:, :, :), Pointer :: RAW_DAT ! raw data array (NSPL1,NSPL2,NTIME)
      Real (DP), Dimension (:, :, :, :), Pointer :: ENS_DAT ! ensemble data array (NENSM,NSPL1,NSPL2,NTIME)
      Character (Len=50), Intent (Out) :: UNITSTR ! units of data variable
    End Subroutine GET_SUBSET
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for grid2basin -- used to get relationship between grids and basins
  Interface
    Subroutine GRID2BASIN (TSERIES_FILE, STANDARDNAME, CELL_METHODS, NSPL1, NSPL2, NRCH, DTYPE, GEN_2DG, GEN_SRF)
      Use nrtype ! variable types (DP, I4B, etc.)
      Use dat_2dgrid ! generic 2-d data structures
      Character (Len=120), Intent (In) :: TSERIES_FILE ! name of time series file
      Character (Len=120), Intent (In) :: STANDARDNAME ! standard_name of data variable
      Character (Len=120), Intent (In) :: CELL_METHODS ! cell methods for data variable
      Integer (I4B), Intent (In) :: NSPL1 ! number of x points
      Integer (I4B), Intent (In) :: NSPL2 ! number of y points
      Integer (I4B), Intent (In) :: NRCH ! number of sub-basins
      Integer (I4B), Intent (In) :: DTYPE ! input data (1-7: p,t,rh,srad,lrad,wind,pressuse)
      Type (GENDAT), Pointer :: GEN_2DG ! data structure for generic 2-d grid
      Type (GENDAT), Pointer, Optional :: GEN_SRF ! generic 2-d grid -- mean surface
    End Subroutine GRID2BASIN
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for grid2basin -- used to get relationship between grids and basins
  Interface
    Subroutine BASLINKAGE (NSPL1, NSPL2, NRCH, DTYPE, GEN_SRF)
      Use nrtype
      Use dat_2dgrid ! generic 2-d data structures
      Integer (I4B), Intent (In) :: NSPL1 ! number of x points
      Integer (I4B), Intent (In) :: NSPL2 ! number of y points
      Integer (I4B), Intent (In) :: NRCH ! number of sub-basins
      Integer (I4B), Intent (In) :: DTYPE ! input data (1-7: p,t,rh,srad,lrad,wind,pressure)
      Type (GENDAT), Pointer, Optional :: GEN_SRF ! generic 2-d grid -- mean surface
    End Subroutine BASLINKAGE
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for satfrc_bas -- used to get saturated areas
  Interface
    Subroutine SATFRC_BAS (IENS, IBAS, NBINS, UAREA, IAREA, IATNB, SAREA)
      Use nrtype
      Integer (I4B), Intent (In) :: IENS ! Ensemble memeber
      Integer (I4B), Intent (In) :: IBAS ! Catchment ID
      Integer (I4B), Intent (In) :: NBINS ! Number of a/tan(b) classes
      Real (DP), Intent (Out) :: UAREA ! catchment area UNINFLUENCED
      Real (DP), Dimension (NBINS), Intent (Out) :: IAREA ! catchment area INFLUENCED
      Real (DP), Dimension (NBINS), Intent (Out) :: IATNB ! representative a/tan(b) values
      Real (DP), Intent (Out) :: SAREA ! catchment area SATURATED
    End Subroutine SATFRC_BAS
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for getusq_rch -- used to retrieve routed flow particles from u/s reaches and
 !                             non-routed flow particles from the current reach
  Interface
    Subroutine GETUSQ_RCH (IENS, JRCH, Q_JRCH, TENTRY, T_EXIT, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: IENS ! ensemble member
      Integer (I4B), Intent (In) :: JRCH ! reach to process
      Real (DP), Dimension (:), Pointer :: Q_JRCH ! merged (non-routed) flow in JRCH
      Real (DP), Dimension (:), Pointer :: TENTRY ! time flow particles entered JRCH
      Real (DP), Dimension (:), Pointer :: T_EXIT ! time flow expected to exit JRCH
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step offset
    End Subroutine GETUSQ_RCH
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for qexmul_rch -- used to extract flow from multiple u/s reaches
  Interface
    Subroutine QEXMUL_RCH (IENS, JRCH, ND, QD, TD, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: IENS ! ensemble member
      Integer (I4B), Intent (In) :: JRCH ! reach to process
      Integer (I4B), Intent (Out) :: ND ! number of routed particles
      Real (DP), Dimension (:), Pointer :: QD ! flow particles just enetered JRCH
      Real (DP), Dimension (:), Pointer :: TD ! time flow particles entered JRCH
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step offset
    End Subroutine QEXMUL_RCH
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for extract_from_reach -- used to extract flow from multiple u/s reaches
  Interface
    Subroutine EXTRACT_FROM_RCH (IENS, JRCH, NR, Q_JRCH, T_EXIT, T_END, TNEW)
      Use nrtype
      Integer (I4B), Intent (In) :: IENS ! ensemble member
      Integer (I4B), Intent (In) :: JRCH ! reach to process
      Integer (I4B), Intent (In) :: NR ! number of routed particles
      Real (DP), Dimension (:), Pointer :: Q_JRCH ! flow in downstream reach JRCH
      Real (DP), Dimension (:), Pointer :: T_EXIT ! time particle expected exit JRCH
      Real (DP) :: T_END ! end of time step
      Real (DP), Dimension (2) :: TNEW ! start/end of time step
    End Subroutine EXTRACT_FROM_RCH
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for remove_rch -- used to extract flow from multiple u/s reaches
  Interface
    Subroutine REMOVE_RCH (Q_JRCH, TENTRY, T_EXIT)
      Use nrtype
      Real (DP), Dimension (:), Pointer :: Q_JRCH ! merged (non-routed) flow in JRCH
      Real (DP), Dimension (:), Pointer :: TENTRY ! time flow particles entered JRCH
      Real (DP), Dimension (:), Pointer :: T_EXIT ! time flow particles exited JRCH
    End Subroutine REMOVE_RCH
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for kinwav_rch -- used to propogate flow particles thru single river segment
  Interface
    Subroutine KINWAV_RCH (JRCH, T_START, T_END, Q_JRCH, TENTRY, FROUTE, T_EXIT, NQ2)! output
      Use nrtype
      Integer (I4B), Intent (In) :: JRCH ! Reach to process
      Real (DP), Intent (In) :: T_START ! start of the time step
      Real (DP), Intent (In) :: T_END ! end of the time step
      Real (DP), Dimension (:), Intent (Inout) :: Q_JRCH ! flow to be routed
      Real (DP), Dimension (:), Intent (Inout) :: TENTRY ! time to be routed
      Logical (LGT), Dimension (:), Intent (Inout) :: FROUTE ! routing flag, T=routed
      Real (DP), Dimension (:), Intent (Inout) :: T_EXIT ! time pts expected exit segment
      Integer (I4B), Intent (Out) :: NQ2 ! # particles (<= input b/c merge)
    End Subroutine KINWAV_RCH
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for interp_rch -- used to compute time step average flow
  Interface
    Subroutine INTERP_RCH (TOLD, QOLD, TNEW, QNEW, IERR)
      Use nrtype
      Real (DP), Dimension (:), Intent (In) :: TOLD ! input time array
      Real (DP), Dimension (:), Intent (In) :: QOLD ! input flow array
      Real (DP), Dimension (:), Intent (In) :: TNEW ! desired output times
      Real (DP), Dimension (:), Intent (Out) :: QNEW ! flow averaged for desired times
      Integer (I4B), Intent (Out) :: IERR ! error, 1= bad bounds
    End Subroutine INTERP_RCH
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for linearcorr -- used to compute the linear correlation coefficient
  Interface
    Subroutine LINEARCORR (X, Y, R)
      Use nrtype
      Real (DP), Intent (Out) :: R ! correlation coefficient
      Real (DP), Dimension (:), Intent (In) :: X, Y ! input vectors
    End Subroutine LINEARCORR
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for sqrtmatrix -- used to compute the square root of a matrix
  Interface
    Subroutine SQRTMATRIX (A, B, N)
      Use nrtype
      Real (DP), Dimension (:, :), Intent (In) :: A ! input matrix
      Real (DP), Dimension (:, :), Intent (Out) :: B ! square root of A (A=B*B)
      Integer (I4B), Intent (In) :: N ! number of elements of A
    End Subroutine SQRTMATRIX
  End Interface
 ! ---------------------------------------------------------------------------------------
 ! interface for inv_matrix -- used to invert a matrix
  Interface
    Subroutine INV_MATRIX (A, B, N)
      Use nrtype
      Real (DP), Dimension (:, :), Intent (In) :: A ! input matrix
      Real (DP), Dimension (:, :), Intent (Out) :: B ! inverse of A (B = A**-1)
      Integer (I4B), Intent (In) :: N ! number of elements of A
    End Subroutine INV_MATRIX
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for raindisagg
  Interface
    Subroutine RAINDISAGG (NUM_DIS, RAIN, NLEVELS, TIMINT)
      Use nrtype
      Integer (I4B), Intent (In) :: NUM_DIS ! Number of disaggregated elements
      Real (DP), Dimension (NUM_DIS), Intent (Inout) :: RAIN ! Rainfall array
      Integer (I4B), Intent (In) :: NLEVELS ! Number of the levels in the cascade
      Real (DP), Intent (In) :: TIMINT ! time interval of data (seconds)
    End Subroutine RAINDISAGG
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for gridsubset
  Interface ! required because of the pointer arguments in the subroutine
    Subroutine GRIDSUBSET (LON2D_A, LAT2D_A, LON2D_B, LAT2D_B, THRSHLD, OFFSET, METHOD, MASK, MISSING, GRIDLIST)
      Use nrtype ! variable types (DP, I4B, etc.)
      Use grid_lists ! structure for grid selection lists
      Real (DP), Dimension (:, :), Pointer :: LON2D_A ! 2d array of lon for basic grid
      Real (DP), Dimension (:, :), Pointer :: LAT2D_A ! 2d array of lat for basic grid
      Real (DP), Dimension (:, :), Pointer :: LON2D_B ! 2d array of lon for input grid
      Real (DP), Dimension (:, :), Pointer :: LAT2D_B ! 2d array of lat for input grid
      Real (DP), Intent (In) :: THRSHLD ! threshold distance
      Real (DP), Intent (In) :: OFFSET ! offset from initial search area
      Integer (I4B), Intent (In) :: METHOD ! method used to select input grid points
      Logical (LGT), Dimension (:, :), Pointer :: MASK ! input data grid points to consider
      Logical (LGT), Intent (In) :: MISSING ! if true missing points allowed
      Type (GRDARRAY), Dimension (:, :), Pointer :: GRIDLIST ! structure for output
    End Subroutine GRIDSUBSET
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for model1step
  Interface ! required because of the optional argument
    Subroutine MODEL1STEP (NENS, NRCH, MRCH, NHYD, NLAK, NOBL, IOFFSET, ISTEP, IFORCE0, IFORCE1, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: NENS ! number of ensemble members
      Integer (I4B), Intent (In) :: NRCH ! number of basins
      Integer (I4B), Intent (In) :: MRCH ! number of basins selected for output
      Integer (I4B), Intent (In) :: NHYD ! number of water level recorders
      Integer (I4B), Intent (In) :: NOBL ! number of lake level observation locations
      Integer (I4B), Intent (In) :: NLAK ! number of lakes
      Integer (I4B), Intent (In) :: IOFFSET ! timestep offset from start of record
      Integer (I4B), Intent (In) :: ISTEP ! loop through time steps
      Integer (I4B), Intent (In) :: IFORCE0 ! starting index for data
      Integer (I4B), Intent (In) :: IFORCE1 ! end index for data
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step offset
    End Subroutine MODEL1STEP
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for q_rch_lake
  Interface ! required because of the optional argument
    Subroutine Q_RCH_LAKE (IENS, JRCH, NOBL, NHYD, ISTEP, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: IENS ! loop through ensemble members
      Integer (I4B), Intent (In) :: JRCH ! loop through the stream segments
      Integer (I4B), Intent (In) :: NOBL ! number of lake level observation locations
      Integer (I4B), Intent (In) :: NHYD ! number of streamq observation locations
      Integer (I4B), Intent (In) :: ISTEP ! loop through time steps
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step offset
    End Subroutine Q_RCH_LAKE
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for qroute_rch
  Interface ! required because of the optional argument
    Subroutine QROUTE_RCH (IENS, JRCH, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: IENS ! ensemble member
      Integer (I4B), Intent (In) :: JRCH ! reach to process
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step offset
    End Subroutine QROUTE_RCH
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface to write_ts_sp
  Interface ! required because of the optional argument
    Subroutine WRITE_TS_SP (ISTEP, NENS, MRCH, NRCH, NHYD, NPCP, NLAK, NOBL, SECTION, RSTEP)
      Use nrtype
    ! input variables
      Integer (I4B), Intent (In) :: ISTEP ! loop through time steps
      Integer (I4B), Intent (In) :: NENS ! number of ensemble members
      Integer (I4B), Intent (In) :: MRCH ! number of basins selected for output
      Integer (I4B), Intent (In) :: NRCH ! number of basins
      Integer (I4B), Intent (In) :: NHYD ! number of water level recorders
      Integer (I4B), Intent (In) :: NPCP ! number of precipitation stations
      Integer (I4B), Intent (In) :: NLAK ! number of lakes
      Integer (I4B), Intent (In) :: NOBL ! number of lake level observation locations
      Integer (I4B), Intent (In), Optional :: SECTION ! section of data to output
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step offset
    End Subroutine WRITE_TS_SP
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface to err_calc
  Interface ! required because of the optional argument
    Subroutine ERR_CALC (ISTEP, IOFFSET, NENS, MRCH, NRCH, NHYD, JRUN, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: ISTEP ! loop through time steps
      Integer (I4B), Intent (In) :: IOFFSET ! timestep offset from start of record
      Integer (I4B), Intent (In) :: NENS ! number of ensemble members
      Integer (I4B), Intent (In) :: MRCH ! number of basins selected for output
      Integer (I4B), Intent (In) :: NRCH ! number of basins
      Integer (I4B), Intent (In) :: NHYD ! number of water level recorders
      Integer (I4B), Intent (In) :: JRUN ! index of the model error parameter set used
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step offset
    End Subroutine ERR_CALC
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for sp_out_put
  Interface ! required because of the optional argument
    Subroutine SP_OUT_PUT (ISTEP, ITIM, NENS, NRCH, NHYD, SECTION, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: ISTEP ! Model time step
      Integer (I4B), Intent (In) :: ITIM ! Index of spatial data on snow covered area
      Integer (I4B), Intent (In) :: NENS ! Number of ensemble members
      Integer (I4B), Intent (In) :: NRCH ! Number of basins
      Integer (I4B), Intent (In) :: NHYD ! Number of gages
      Integer (I4B), Intent (In) :: SECTION ! Section of variables to write out
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step
    End Subroutine SP_OUT_PUT
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for ts_out_put
  Interface ! required because of the optional argument
    Subroutine TS_OUT_PUT (ISTEP, NENS, MRCH, NRCH, NHYD, NPCP, NLAK, NOBL, SECTION, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: ISTEP ! Model time step
      Integer (I4B), Intent (In) :: NENS ! Number of ensemble members
      Integer (I4B), Intent (In) :: MRCH ! Number of selected basins
      Integer (I4B), Intent (In) :: NRCH ! Number of basins
      Integer (I4B), Intent (In) :: NHYD ! Number of water level recorders
      Integer (I4B), Intent (In) :: NOBL ! Number of lake level observation locations
      Integer (I4B), Intent (In) :: NPCP ! Number of precipitation stations
      Integer (I4B), Intent (In) :: NLAK ! Number of lakes
      Integer (I4B), Intent (In) :: SECTION ! Section of variables to write out
      Integer (I4B), Intent (In), Optional :: RSTEP ! retrospective time step
    End Subroutine TS_OUT_PUT
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for prberr_put
  Interface ! required because of the optional argument
    Subroutine PRBERR_PUT (IOFFSET, ISTEP, JRUN, NHYD, RSTEP)
      Use nrtype
      Integer (I4B), Intent (In) :: IOFFSET ! Offset for data write statement
      Integer (I4B), Intent (In) :: ISTEP ! Model time step
      Integer (I4B), Intent (In) :: JRUN ! Index of model error parameter set
      Integer (I4B), Intent (In) :: NHYD ! Number of gages
      Integer (I4B), Intent (In), Optional :: RSTEP ! Retrospective time step
    End Subroutine PRBERR_PUT
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for rstart_get
  Interface ! required because of the optional argument
    Subroutine RSTART_GET (NENS, NRCH, NLAK, DTIME, LEXIST, FEXIST)
      Use nrtype
      Integer (I4B), Intent (In) :: NENS ! Number of ensemble members
      Integer (I4B), Intent (In) :: NRCH ! Number of basins
      Integer (I4B), Intent (In) :: NLAK ! Number of lakes
      Real (DP), Intent (In) :: DTIME ! Time difference between restart file and model reference time
      Logical (LGT), Intent (Out) :: LEXIST ! .TRUE. if re-start file exists
      Logical (LGT), Intent (In), Optional :: FEXIST ! if .TRUE., only check for file existance
    End Subroutine RSTART_GET
  End Interface
 ! ----------------------------------------------------------------------------------------
 ! interface for moddat_get
  Interface ! required because of the optional argument
    Subroutine MODDAT_GET (NENS, NRCH, DTIME, LEXIST, FEXIST)
      Use nrtype
      Integer (I4B), Intent (In) :: NENS ! Number of ensemble members
      Integer (I4B), Intent (In) :: NRCH ! Number of basins
      Real (DP), Intent (In) :: DTIME ! Time difference between data_save file and model reference time
      Logical (LGT), Intent (Out) :: LEXIST ! .TRUE. if re-start file exists
      Logical (LGT), Intent (In), Optional :: FEXIST ! if .TRUE., only check for file existance
    End Subroutine MODDAT_GET
  End Interface
! ----------------------------------------------------------------------------------------
End Module interblock
