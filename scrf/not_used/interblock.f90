module interblock
  implicit none
 ! ---------------------------------------------------------------------------------------
 ! DEFINES EXPLICIT INTERFACES BETWEEN SUB-PROGRAMS
 ! ---------------------------------------------------------------------------------------
 ! interface for get_2dgrid -- used to get a generic 2-d grid
  interface
    subroutine get_2dgrid (tseries_file, standardname, cell_methods, gen_2dg, nensm, nspl1, nspl2, &
   & ntime, idat0, idat1, isim0, isim1, unitstr, vexist, texist, eexist)
      use nrtype ! variable types etc.
      use dat_2dgrid ! generic 2-d grid
      character (len=120), intent (in) :: tseries_file ! name of time series file
      character (len=120), intent (in) :: standardname ! standard_name of data variable
      character (len=120), intent (in) :: cell_methods ! cell methods for data variable
      type (gendat), pointer :: gen_2dg ! 2-d data structure
      integer (i4b), intent (in) :: nensm ! # ensemble members
      integer (i4b), intent (out) :: nspl1 ! # points 1st spatial dimension
      integer (i4b), intent (out) :: nspl2 ! # points 2nd spatial dimension
      integer (i4b), intent (out) :: ntime ! number of time steps
      integer (i4b), intent (out) :: idat0 ! first index in data array (time vector)
      integer (i4b), intent (out) :: idat1 ! last index in data array (time vector)
      integer (i4b), intent (out) :: isim0 ! first index of simulation array
      integer (i4b), intent (out) :: isim1 ! last index of simulation array
      character (len=50), intent (out) :: unitstr ! units of data variable
      logical (lgt), intent (out) :: vexist ! flag if the variable exists
      logical (lgt), intent (out) :: texist ! flag if variable depends on time
      logical (lgt), intent (out) :: eexist ! flag if ensembles exist
    end subroutine get_2dgrid
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_station -- used to get lat-lon from station data
  interface
    subroutine xy_station (tseries_file, yvar_id, xvar_id, xy_lat, xy_lon, spl1_start, spl1_count, &
   & spl2_start, spl2_count, var_id)
      use nrtype ! variable types etc.
      character (len=120), intent (in) :: tseries_file ! name of time series file
      integer (i4b), intent (in) :: yvar_id ! NetCDF ID for latitude
      integer (i4b), intent (in) :: xvar_id ! NetCDF ID for longitude
      integer (i4b), intent (in) :: var_id ! NetCDF ID for variable
      real (dp), dimension (:, :), pointer :: xy_lat ! 2-dimensional grid of latitudes
      real (dp), dimension (:, :), pointer :: xy_lon ! 2-dimensional grid of longitudes
      integer (i4b), dimension (:), pointer :: spl1_start ! start index for 1st spatial dimension
      integer (i4b), intent (out) :: spl1_count ! count for 1st spatial dimension
      integer (i4b), dimension (:), pointer :: spl2_start ! start index for 2nd spatial dimension
      integer (i4b), intent (out) :: spl2_count ! count for 2nd spatial dimension
    end subroutine xy_station
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_lat_lon -- used to get lat-lon from a regular grid
  interface
    subroutine xy_lat_lon (tseries_file, yvar_id, xvar_id, xy_lat, xy_lon, spl1_start, spl1_count, &
   & spl2_start, spl2_count)
      use nrtype ! variable types etc.
      character (len=120), intent (in) :: tseries_file ! name of time series file
      integer (i4b), intent (in) :: yvar_id ! NetCDF ID for latitude
      integer (i4b), intent (in) :: xvar_id ! NetCDF ID for longitude
      real (dp), dimension (:, :), pointer :: xy_lat ! 2-dimensional grid of latitudes
      real (dp), dimension (:, :), pointer :: xy_lon ! 2-dimensional grid of longitudes
      integer (i4b), dimension (:), pointer :: spl1_start ! start index for 1st spatial dimension
      integer (i4b), intent (out) :: spl1_count ! count for 1st spatial dimension
      integer (i4b), dimension (:), pointer :: spl2_start ! start index for 2nd spatial dimension
      integer (i4b), intent (out) :: spl2_count ! count for 2nd spatial dimension
    end subroutine xy_lat_lon
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for xy_rotated -- used to get lat-lon from a rotated grid
  interface
    subroutine xy_rotated (tseries_file, yvar_id, xvar_id, xy_lat, xy_lon, spl1_start, spl1_count, &
   & spl2_start, spl2_count)
      use nrtype ! variable types etc.
      character (len=120), intent (in) :: tseries_file ! name of time series file
      integer (i4b), intent (in) :: yvar_id ! NetCDF ID for latitude
      integer (i4b), intent (in) :: xvar_id ! NetCDF ID for longitude
      real (dp), dimension (:, :), pointer :: xy_lat ! 2-dimensional grid of latitudes
      real (dp), dimension (:, :), pointer :: xy_lon ! 2-dimensional grid of longitudes
      integer (i4b), dimension (:), pointer :: spl1_start ! start index for 1st spatial dimension
      integer (i4b), intent (out) :: spl1_count ! count for 1st spatial dimension
      integer (i4b), dimension (:), pointer :: spl2_start ! start index for 2nd spatial dimension
      integer (i4b), intent (out) :: spl2_count ! count for 2nd spatial dimension
    end subroutine xy_rotated
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_elevtn -- used to get elevation data
  interface
    subroutine get_elevtn (tseries_file, csystem, xy_elv, spl1_start, spl1_count, spl2_start, &
   & spl2_count)
      use nrtype ! variable types etc.
      character (len=120), intent (in) :: tseries_file ! name of time series file
      character (len=12), intent (in) :: csystem ! Coordinate system of raw data
      integer (i4b), dimension (:), pointer :: spl1_start ! start index for 1st spatial dimension
      integer (i4b), intent (in) :: spl1_count ! count for 1st spatial dimension
      integer (i4b), dimension (:), pointer :: spl2_start ! start index for 2nd spatial dimension
      integer (i4b), intent (in) :: spl2_count ! count for 2nd spatial dimension
      real (dp), dimension (:, :), pointer :: xy_elv ! 2-dimensional grid of elevation
    end subroutine get_elevtn
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_timdat -- used to get full time arrays
  interface
    subroutine get_timdat (tseries_file, tvar_id, tb0, tb1, tim, bex, ntim_data)
      use nrtype ! variable types etc.
      character (len=120), intent (in) :: tseries_file ! name of time series file
      integer (i4b), intent (in) :: tvar_id ! NetCDF ID for time
      real (dp), dimension (:), pointer :: tb0 ! start of time interval (in data file)
      real (dp), dimension (:), pointer :: tb1 ! end of time interval (in data file)
      real (dp), dimension (:), pointer :: tim ! time stamp (in data file)
      logical (lgt), intent (out) :: bex ! .TRUE. if time bounds exist
      integer (i4b), intent (out) :: ntim_data ! number of time intervals in the data file
    end subroutine get_timdat
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_subtim -- used to get time subset
  interface
    subroutine get_subtim (tb0, tb1, tim, bex, ntim, time_start, time_count, isim0, isim1)
      use nrtype ! variable types etc.
      real (dp), dimension (:), pointer :: tb0 ! start of time interval (in data file)
      real (dp), dimension (:), pointer :: tb1 ! end of time interval (in data file)
      real (dp), dimension (:), pointer :: tim ! time stamp (in data file)
      logical (lgt), intent (in) :: bex ! .TRUE. if time bounds exist
      integer (i4b), intent (in) :: ntim ! number of time intervals in data
      integer (i4b), intent (out) :: time_start ! start index for time
      integer (i4b), intent (out) :: time_count ! count for time
      integer (i4b), intent (out) :: isim0 ! first index of simulation array
      integer (i4b), intent (out) :: isim1 ! last index of simulation array
    end subroutine get_subtim
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for get_subset -- used to get data subset
  interface
    subroutine get_subset (tseries_file, csystem, var_id, tdim_ix, zdim_ix, ydim_ix, xdim_ix, &
   & edim_ix, time_start, time_count, spl1_start, spl1_count, spl2_start, spl2_count, vert_start, &
   & vert_count, iens_start, iens_count, raw_dat, ens_dat, unitstr)
      use nrtype ! variable types etc.
      character (len=120), intent (in) :: tseries_file ! name of time series file
      character (len=12), intent (in) :: csystem ! Coordinate system of raw data
      integer (i4b), intent (in) :: var_id ! NetCDF ID for data variable
      integer (i4b), intent (in) :: tdim_ix ! NetCDF dimension index for time variable
      integer (i4b), intent (in) :: zdim_ix ! NetCDF dimension index for (height, depth)
      integer (i4b), intent (in) :: ydim_ix ! NetCDF dimension index for latitude
      integer (i4b), intent (in) :: xdim_ix ! NetCDF dimension index for longitude
      integer (i4b), intent (in) :: edim_ix ! NetCDF dimension index for ensemble
      integer (i4b), intent (in) :: time_start ! start index for time
      integer (i4b), intent (in) :: time_count ! count for time
      integer (i4b), dimension (:), pointer :: spl1_start ! start index for 1st spatial dimension
      integer (i4b), intent (in) :: spl1_count ! count for 1st spatial dimension
      integer (i4b), dimension (:), pointer :: spl2_start ! start index for 2nd spatial dimension
      integer (i4b), intent (in) :: spl2_count ! count for 2nd spatial dimension
      integer (i4b), intent (in) :: vert_start ! start index for vertical dimension
      integer (i4b), intent (in) :: vert_count ! count for vertical dimension
      integer (i4b), intent (in) :: iens_start ! start index for ensemble dimension
      integer (i4b), intent (in) :: iens_count ! count for ensemble dimension
      real (dp), dimension (:, :, :), pointer :: raw_dat ! raw data array (NSPL1,NSPL2,NTIME)
      real (dp), dimension (:, :, :, :), pointer :: ens_dat ! ensemble data array (NENSM,NSPL1,NSPL2,NTIME)
      character (len=50), intent (out) :: unitstr ! units of data variable
    end subroutine get_subset
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for grid2basin -- used to get relationship between grids and basins
  interface
    subroutine grid2basin (tseries_file, standardname, cell_methods, nspl1, nspl2, nrch, dtype, &
   & gen_2dg, gen_srf)
      use nrtype ! variable types (DP, I4B, etc.)
      use dat_2dgrid ! generic 2-d data structures
      character (len=120), intent (in) :: tseries_file ! name of time series file
      character (len=120), intent (in) :: standardname ! standard_name of data variable
      character (len=120), intent (in) :: cell_methods ! cell methods for data variable
      integer (i4b), intent (in) :: nspl1 ! number of x points
      integer (i4b), intent (in) :: nspl2 ! number of y points
      integer (i4b), intent (in) :: nrch ! number of sub-basins
      integer (i4b), intent (in) :: dtype ! input data (1-7: p,t,rh,srad,lrad,wind,pressuse)
      type (gendat), pointer :: gen_2dg ! data structure for generic 2-d grid
      type (gendat), pointer, optional :: gen_srf ! generic 2-d grid -- mean surface
    end subroutine grid2basin
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for grid2basin -- used to get relationship between grids and basins
  interface
    subroutine baslinkage (nspl1, nspl2, nrch, dtype, gen_srf)
      use nrtype
      use dat_2dgrid ! generic 2-d data structures
      integer (i4b), intent (in) :: nspl1 ! number of x points
      integer (i4b), intent (in) :: nspl2 ! number of y points
      integer (i4b), intent (in) :: nrch ! number of sub-basins
      integer (i4b), intent (in) :: dtype ! input data (1-7: p,t,rh,srad,lrad,wind,pressure)
      type (gendat), pointer, optional :: gen_srf ! generic 2-d grid -- mean surface
    end subroutine baslinkage
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for satfrc_bas -- used to get saturated areas
  interface
    subroutine satfrc_bas (iens, ibas, nbins, uarea, iarea, iatnb, sarea)
      use nrtype
      integer (i4b), intent (in) :: iens ! Ensemble memeber
      integer (i4b), intent (in) :: ibas ! Catchment ID
      integer (i4b), intent (in) :: nbins ! Number of a/tan(b) classes
      real (dp), intent (out) :: uarea ! catchment area UNINFLUENCED
      real (dp), dimension (nbins), intent (out) :: iarea ! catchment area INFLUENCED
      real (dp), dimension (nbins), intent (out) :: iatnb ! representative a/tan(b) values
      real (dp), intent (out) :: sarea ! catchment area SATURATED
    end subroutine satfrc_bas
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for getusq_rch -- used to retrieve routed flow particles from u/s reaches and
 !                             non-routed flow particles from the current reach
  interface
    subroutine getusq_rch (iens, jrch, q_jrch, tentry, t_exit, rstep)
      use nrtype
      integer (i4b), intent (in) :: iens ! ensemble member
      integer (i4b), intent (in) :: jrch ! reach to process
      real (dp), dimension (:), pointer :: q_jrch ! merged (non-routed) flow in JRCH
      real (dp), dimension (:), pointer :: tentry ! time flow particles entered JRCH
      real (dp), dimension (:), pointer :: t_exit ! time flow expected to exit JRCH
      integer (i4b), intent (in), optional :: rstep ! retrospective time step offset
    end subroutine getusq_rch
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for qexmul_rch -- used to extract flow from multiple u/s reaches
  interface
    subroutine qexmul_rch (iens, jrch, nd, qd, td, rstep)
      use nrtype
      integer (i4b), intent (in) :: iens ! ensemble member
      integer (i4b), intent (in) :: jrch ! reach to process
      integer (i4b), intent (out) :: nd ! number of routed particles
      real (dp), dimension (:), pointer :: qd ! flow particles just enetered JRCH
      real (dp), dimension (:), pointer :: td ! time flow particles entered JRCH
      integer (i4b), intent (in), optional :: rstep ! retrospective time step offset
    end subroutine qexmul_rch
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for extract_from_reach -- used to extract flow from multiple u/s reaches
  interface
    subroutine extract_from_rch (iens, jrch, nr, q_jrch, t_exit, t_end, tnew)
      use nrtype
      integer (i4b), intent (in) :: iens ! ensemble member
      integer (i4b), intent (in) :: jrch ! reach to process
      integer (i4b), intent (in) :: nr ! number of routed particles
      real (dp), dimension (:), pointer :: q_jrch ! flow in downstream reach JRCH
      real (dp), dimension (:), pointer :: t_exit ! time particle expected exit JRCH
      real (dp) :: t_end ! end of time step
      real (dp), dimension (2) :: tnew ! start/end of time step
    end subroutine extract_from_rch
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for remove_rch -- used to extract flow from multiple u/s reaches
  interface
    subroutine remove_rch (q_jrch, tentry, t_exit)
      use nrtype
      real (dp), dimension (:), pointer :: q_jrch ! merged (non-routed) flow in JRCH
      real (dp), dimension (:), pointer :: tentry ! time flow particles entered JRCH
      real (dp), dimension (:), pointer :: t_exit ! time flow particles exited JRCH
    end subroutine remove_rch
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for kinwav_rch -- used to propogate flow particles thru single river segment
  interface
    subroutine kinwav_rch (jrch, t_start, t_end, q_jrch, tentry, froute, t_exit, nq2)! output
      use nrtype
      integer (i4b), intent (in) :: jrch ! Reach to process
      real (dp), intent (in) :: t_start ! start of the time step
      real (dp), intent (in) :: t_end ! end of the time step
      real (dp), dimension (:), intent (inout) :: q_jrch ! flow to be routed
      real (dp), dimension (:), intent (inout) :: tentry ! time to be routed
      logical (lgt), dimension (:), intent (inout) :: froute ! routing flag, T=routed
      real (dp), dimension (:), intent (inout) :: t_exit ! time pts expected exit segment
      integer (i4b), intent (out) :: nq2 ! # particles (<= input b/c merge)
    end subroutine kinwav_rch
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for interp_rch -- used to compute time step average flow
  interface
    subroutine interp_rch (told, qold, tnew, qnew, ierr)
      use nrtype
      real (dp), dimension (:), intent (in) :: told ! input time array
      real (dp), dimension (:), intent (in) :: qold ! input flow array
      real (dp), dimension (:), intent (in) :: tnew ! desired output times
      real (dp), dimension (:), intent (out) :: qnew ! flow averaged for desired times
      integer (i4b), intent (out) :: ierr ! error, 1= bad bounds
    end subroutine interp_rch
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for linearcorr -- used to compute the linear correlation coefficient
  interface
    subroutine linearcorr (x, y, r)
      use nrtype
      real (dp), intent (out) :: r ! correlation coefficient
      real (dp), dimension (:), intent (in) :: x, y ! input vectors
    end subroutine linearcorr
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for sqrtmatrix -- used to compute the square root of a matrix
  interface
    subroutine sqrtmatrix (a, b, n)
      use nrtype
      real (dp), dimension (:, :), intent (in) :: a ! input matrix
      real (dp), dimension (:, :), intent (out) :: b ! square root of A (A=B*B)
      integer (i4b), intent (in) :: n ! number of elements of A
    end subroutine sqrtmatrix
  end interface
 ! ---------------------------------------------------------------------------------------
 ! interface for inv_matrix -- used to invert a matrix
  interface
    subroutine inv_matrix (a, b, n)
      use nrtype
      real (dp), dimension (:, :), intent (in) :: a ! input matrix
      real (dp), dimension (:, :), intent (out) :: b ! inverse of A (B = A**-1)
      integer (i4b), intent (in) :: n ! number of elements of A
    end subroutine inv_matrix
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for raindisagg
  interface
    subroutine raindisagg (num_dis, rain, nlevels, timint)
      use nrtype
      integer (i4b), intent (in) :: num_dis ! Number of disaggregated elements
      real (dp), dimension (num_dis), intent (inout) :: rain ! Rainfall array
      integer (i4b), intent (in) :: nlevels ! Number of the levels in the cascade
      real (dp), intent (in) :: timint ! time interval of data (seconds)
    end subroutine raindisagg
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for gridsubset
  interface ! required because of the pointer arguments in the subroutine
    subroutine gridsubset (lon2d_a, lat2d_a, lon2d_b, lat2d_b, thrshld, offset, method, mask, &
   & missing, gridlist)
      use nrtype ! variable types (DP, I4B, etc.)
      use grid_lists ! structure for grid selection lists
      real (dp), dimension (:, :), pointer :: lon2d_a ! 2d array of lon for basic grid
      real (dp), dimension (:, :), pointer :: lat2d_a ! 2d array of lat for basic grid
      real (dp), dimension (:, :), pointer :: lon2d_b ! 2d array of lon for input grid
      real (dp), dimension (:, :), pointer :: lat2d_b ! 2d array of lat for input grid
      real (dp), intent (in) :: thrshld ! threshold distance
      real (dp), intent (in) :: offset ! offset from initial search area
      integer (i4b), intent (in) :: method ! method used to select input grid points
      logical (lgt), dimension (:, :), pointer :: mask ! input data grid points to consider
      logical (lgt), intent (in) :: missing ! if true missing points allowed
      type (grdarray), dimension (:, :), pointer :: gridlist ! structure for output
    end subroutine gridsubset
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for model1step
  interface ! required because of the optional argument
    subroutine model1step (nens, nrch, mrch, nhyd, nlak, nobl, ioffset, istep, iforce0, iforce1, &
   & rstep)
      use nrtype
      integer (i4b), intent (in) :: nens ! number of ensemble members
      integer (i4b), intent (in) :: nrch ! number of basins
      integer (i4b), intent (in) :: mrch ! number of basins selected for output
      integer (i4b), intent (in) :: nhyd ! number of water level recorders
      integer (i4b), intent (in) :: nobl ! number of lake level observation locations
      integer (i4b), intent (in) :: nlak ! number of lakes
      integer (i4b), intent (in) :: ioffset ! timestep offset from start of record
      integer (i4b), intent (in) :: istep ! loop through time steps
      integer (i4b), intent (in) :: iforce0 ! starting index for data
      integer (i4b), intent (in) :: iforce1 ! end index for data
      integer (i4b), intent (in), optional :: rstep ! retrospective time step offset
    end subroutine model1step
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for q_rch_lake
  interface ! required because of the optional argument
    subroutine q_rch_lake (iens, jrch, nobl, nhyd, istep, rstep)
      use nrtype
      integer (i4b), intent (in) :: iens ! loop through ensemble members
      integer (i4b), intent (in) :: jrch ! loop through the stream segments
      integer (i4b), intent (in) :: nobl ! number of lake level observation locations
      integer (i4b), intent (in) :: nhyd ! number of streamq observation locations
      integer (i4b), intent (in) :: istep ! loop through time steps
      integer (i4b), intent (in), optional :: rstep ! retrospective time step offset
    end subroutine q_rch_lake
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for qroute_rch
  interface ! required because of the optional argument
    subroutine qroute_rch (iens, jrch, rstep)
      use nrtype
      integer (i4b), intent (in) :: iens ! ensemble member
      integer (i4b), intent (in) :: jrch ! reach to process
      integer (i4b), intent (in), optional :: rstep ! retrospective time step offset
    end subroutine qroute_rch
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface to write_ts_sp
  interface ! required because of the optional argument
    subroutine write_ts_sp (istep, nens, mrch, nrch, nhyd, npcp, nlak, nobl, section, rstep)
      use nrtype
    ! input variables
      integer (i4b), intent (in) :: istep ! loop through time steps
      integer (i4b), intent (in) :: nens ! number of ensemble members
      integer (i4b), intent (in) :: mrch ! number of basins selected for output
      integer (i4b), intent (in) :: nrch ! number of basins
      integer (i4b), intent (in) :: nhyd ! number of water level recorders
      integer (i4b), intent (in) :: npcp ! number of precipitation stations
      integer (i4b), intent (in) :: nlak ! number of lakes
      integer (i4b), intent (in) :: nobl ! number of lake level observation locations
      integer (i4b), intent (in), optional :: section ! section of data to output
      integer (i4b), intent (in), optional :: rstep ! retrospective time step offset
    end subroutine write_ts_sp
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface to err_calc
  interface ! required because of the optional argument
    subroutine err_calc (istep, ioffset, nens, mrch, nrch, nhyd, jrun, rstep)
      use nrtype
      integer (i4b), intent (in) :: istep ! loop through time steps
      integer (i4b), intent (in) :: ioffset ! timestep offset from start of record
      integer (i4b), intent (in) :: nens ! number of ensemble members
      integer (i4b), intent (in) :: mrch ! number of basins selected for output
      integer (i4b), intent (in) :: nrch ! number of basins
      integer (i4b), intent (in) :: nhyd ! number of water level recorders
      integer (i4b), intent (in) :: jrun ! index of the model error parameter set used
      integer (i4b), intent (in), optional :: rstep ! retrospective time step offset
    end subroutine err_calc
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for sp_out_put
  interface ! required because of the optional argument
    subroutine sp_out_put (istep, itim, nens, nrch, nhyd, section, rstep)
      use nrtype
      integer (i4b), intent (in) :: istep ! Model time step
      integer (i4b), intent (in) :: itim ! Index of spatial data on snow covered area
      integer (i4b), intent (in) :: nens ! Number of ensemble members
      integer (i4b), intent (in) :: nrch ! Number of basins
      integer (i4b), intent (in) :: nhyd ! Number of gages
      integer (i4b), intent (in) :: section ! Section of variables to write out
      integer (i4b), intent (in), optional :: rstep ! retrospective time step
    end subroutine sp_out_put
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for ts_out_put
  interface ! required because of the optional argument
    subroutine ts_out_put (istep, nens, mrch, nrch, nhyd, npcp, nlak, nobl, section, rstep)
      use nrtype
      integer (i4b), intent (in) :: istep ! Model time step
      integer (i4b), intent (in) :: nens ! Number of ensemble members
      integer (i4b), intent (in) :: mrch ! Number of selected basins
      integer (i4b), intent (in) :: nrch ! Number of basins
      integer (i4b), intent (in) :: nhyd ! Number of water level recorders
      integer (i4b), intent (in) :: nobl ! Number of lake level observation locations
      integer (i4b), intent (in) :: npcp ! Number of precipitation stations
      integer (i4b), intent (in) :: nlak ! Number of lakes
      integer (i4b), intent (in) :: section ! Section of variables to write out
      integer (i4b), intent (in), optional :: rstep ! retrospective time step
    end subroutine ts_out_put
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for prberr_put
  interface ! required because of the optional argument
    subroutine prberr_put (ioffset, istep, jrun, nhyd, rstep)
      use nrtype
      integer (i4b), intent (in) :: ioffset ! Offset for data write statement
      integer (i4b), intent (in) :: istep ! Model time step
      integer (i4b), intent (in) :: jrun ! Index of model error parameter set
      integer (i4b), intent (in) :: nhyd ! Number of gages
      integer (i4b), intent (in), optional :: rstep ! Retrospective time step
    end subroutine prberr_put
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for rstart_get
  interface ! required because of the optional argument
    subroutine rstart_get (nens, nrch, nlak, dtime, lexist, fexist)
      use nrtype
      integer (i4b), intent (in) :: nens ! Number of ensemble members
      integer (i4b), intent (in) :: nrch ! Number of basins
      integer (i4b), intent (in) :: nlak ! Number of lakes
      real (dp), intent (in) :: dtime ! Time difference between restart file and model reference time
      logical (lgt), intent (out) :: lexist ! .TRUE. if re-start file exists
      logical (lgt), intent (in), optional :: fexist ! if .TRUE., only check for file existance
    end subroutine rstart_get
  end interface
 ! ----------------------------------------------------------------------------------------
 ! interface for moddat_get
  interface ! required because of the optional argument
    subroutine moddat_get (nens, nrch, dtime, lexist, fexist)
      use nrtype
      integer (i4b), intent (in) :: nens ! Number of ensemble members
      integer (i4b), intent (in) :: nrch ! Number of basins
      real (dp), intent (in) :: dtime ! Time difference between data_save file and model reference time
      logical (lgt), intent (out) :: lexist ! .TRUE. if re-start file exists
      logical (lgt), intent (in), optional :: fexist ! if .TRUE., only check for file existance
    end subroutine moddat_get
  end interface
! ----------------------------------------------------------------------------------------
end module interblock
