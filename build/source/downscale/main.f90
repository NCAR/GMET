! Main driver / program file for the ensemble regression component of the 
! Gridded Meteorological Ensemble Tool (GMET)

program gmet
  use type
  use string_mod
  use utim
  implicit none
 
  interface
    subroutine get_time_list (startdate, enddate, t)
      use utim
      use type
      character (len=100), intent (in) :: startdate, enddate
      real (dp), allocatable, intent (out) :: t (:)
    end subroutine get_time_list
 
    subroutine read_refcst (startdate, enddate, file_var, perturbation, var_name, forecast, v, x, &
   & y, t, error)
      use type
      character (len=100), intent (in) :: startdate, enddate, file_var, perturbation
      character (len=*), intent (in) :: var_name
      integer (i4b), intent (in) :: forecast
      real (dp), allocatable, intent (out) :: v (:, :), x (:), y (:)
      real (dp), allocatable, intent (out) :: t (:)
      integer, intent (out) :: error
    end subroutine read_refcst
 
    subroutine read_station_list (file_name, id, name, lat, lon, elev, sslp_n, sslp_e, n_stations, &
   & error, vars) !AWW
      use type
      character (len=500), intent (in) :: file_name
      character (len=100), allocatable, intent (out) :: id (:), name (:)
      character (len=2), allocatable, intent (out) :: vars (:) !AWW holds PT identifiers
      real (dp), allocatable, intent (out) :: lat (:), lon (:), elev (:), sslp_n (:), sslp_e (:)
      integer (i4b), intent (out) :: n_stations
      integer, intent (out) :: error
    end subroutine read_station_list
 
    subroutine read_grid_list (file_name, lats, lons, elevs, slp_n, slp_e, nx, ny, error)
      use type
      character (len=500), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lats (:), lons (:), elevs (:), slp_n (:), slp_e (:)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_grid_list
 
    subroutine read_domain_grid (grid_list, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
      use netcdf
      use type
 
      character (len=500), intent (in) :: grid_list
      real (dp), allocatable, intent (out) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), &
     & grad_e (:, :), mask (:, :)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_domain_grid
 
    subroutine estimate_coefficients (d, nvars, lats, lons, times, st_rec, end_rec, stnid, stnlat, &
   & stnlon, stnvar, directory, c, poc, error)
      use type
      real (dp), intent (in) :: d (:, :, :), lats (:), lons (:)
      real (dp), intent (in) :: times (:)
      integer (i4b), intent (in) :: nvars
      integer (i4b), intent (in) :: st_rec, end_rec
      character (len=100), intent (in) :: stnid (:)
      real (dp), intent (in) :: stnlat (:), stnlon (:)
      character (len=100), intent (in) :: stnvar
      character (len=500), intent (in) :: directory !AWW added
      real (dp), allocatable, intent (out) :: c (:, :, :), poc (:, :, :)
      integer, intent (out) :: error
    end subroutine estimate_coefficients
 
    subroutine save_coefficients (n_vars, var_names, coefs, startdate, enddate, times, site_list, &
   & stnid, stnlat, stnlon, stnelev, forecast, file, error)
      use netcdf
      use type
      character (len=100), intent (in) :: var_names (:), stnid (:)
      integer, intent (in) :: n_vars, forecast
      character (len=100), intent (in) :: startdate, enddate
      character (len=500), intent (in) :: file, site_list
      real (dp), intent (in) :: stnlat (:), stnlon (:), stnelev (:)
      real (dp), intent (in) :: coefs (:, :, :)
      real (dp), intent (in) :: times (:)
      integer, intent (out) :: error
    end subroutine save_coefficients
 
    subroutine estimate_forcing_regression (nPredict, gen_sta_weights, sta_weight_name, nwp_input_list, &
   & n_nwp, nwp_vars, nwp_prcp_var, x, z, ngrid, maxdistance, times, st_rec, end_rec, &
   & stnid, stnvar, directory, kfold_trials, kfold_hold, pcp, pop, pcperr, obs_max_pcp, tmean, tmean_err, trange, &
   & trange_err, mean_autocorr, mean_tp_corr, error, pcp_2, &
   & pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2, use_stn_weights)
      use type
      character (len=500), intent(in)   :: gen_sta_weights        ! station weight generation flag
      character (len = 500), intent(in) :: sta_weight_name        ! station weight file name
      character (len=500), intent(in)   :: use_stn_weights        ! station weight usage option
      character (len = 2000),intent(in) :: nwp_input_list         ! name of file containint list of NWP input files
      character (len=100), intent(in)   :: nwp_vars(:)            !list of nwp predictor variables
      character (len=100), intent(in)   :: nwp_prcp_var           !name of nwp precipitation predictor
      integer(i4b), intent(inout)          :: nPredict            ! number of total predictors
      integer(i4b), intent(in)             :: n_nwp               ! number of NWP predictors
      real (dp), intent (inout) :: x (:, :), z (:, :)
      real (dp), intent (in) :: maxdistance
      integer (i4b), intent (in) :: ngrid
      real (dp), intent (in) :: times (:)
      integer (i4b), intent (in) :: st_rec, end_rec
      character (len=100), intent (in) :: stnid (:)
      character (len=100), intent (in) :: stnvar
      character (len=500), intent (in) :: directory
      integer(I4B), intent(in)          :: kfold_trials        !number of kfold xval trials
      integer(I4B), intent(in)          :: kfold_hold          !number of stations to withhold from regression
      real (sp), allocatable, intent (out) :: pcp (:, :), pop (:, :), pcperr (:, :)
      real (sp), allocatable, intent (out) :: pcp_2 (:, :), pop_2 (:, :), pcperr_2 (:, :)
      real (sp), allocatable, intent (out) :: tmean (:, :), tmean_err (:, :)        ! OLS tmean estimate and error
      real (sp), allocatable, intent (out) :: trange (:, :), trange_err (:, :)      ! OLS trange estimate and error
      real (sp), allocatable, intent (out) :: tmean_2 (:, :), tmean_err_2 (:, :)    ! OLS tmean estimate and error
      real (sp), allocatable, intent (out) :: trange_2 (:, :), trange_err_2 (:, :)  ! OLS trange estimate and error
      integer, intent (out) :: error
      real (dp), intent (out) :: mean_autocorr (:)  ! mean autocorrelation from all stations over entire time period
      real (dp), intent (out) :: mean_tp_corr (:)   ! mean correlation for mean temp and precip
      real (dp), intent (out) :: obs_max_pcp (:, :) ! max of normalized time step precip
    end subroutine estimate_forcing_regression
 
    ! modified AJN Sept 2013
    ! subroutine save_precip(pcp, pop, pcperror, tmean, tmean_err, trange, trange_err, &
    ! AWW Feb2016, renamed to be more descriptive.  Note did NOT call scrf/save_precip.f90 subroutine
    subroutine save_forcing_regression (pcp, pop, pcperror, obs_max_pcp, tmean, tmean_err, trange, trange_err, &
   & nx, ny, grdlat, grdlon, grdelev, times, mean_autocorr, mean_tp_corr, &
   & file, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, &
   & trange_err_2)
      use netcdf
      use type
      real (sp), intent (in) :: pcp (:, :), pop (:, :), pcperror (:, :)
      real (sp), intent (in) :: tmean (:, :), tmean_err (:, :), trange (:, :), trange_err (:, :)
      real (sp), intent (in) :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
      real (sp), intent (in) :: tmean_2 (:, :), tmean_err_2 (:, :), trange_2 (:, :), trange_err_2 &
     & (:, :)
      integer (i4b), intent (in) :: nx, ny
      real (dp), intent (in) :: grdlat (:), grdlon (:), grdelev (:)
      real (dp), intent (in) :: times (:)
      real (dp), intent (in) :: mean_autocorr (:), mean_tp_corr (:)
      !real (dp), intent (in) :: y_mean (:, :), y_std (:, :), y_std_all (:, :)
      !real (dp), intent (in) :: y_min (:, :), y_max (:, :)
      real (dp), intent (in) :: obs_max_pcp (:, :)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_forcing_regression
 
  end interface
  ! === end of interface, start the program ====

  ! Local variables
 
  character (len=100) :: config_file
  integer, parameter  :: nconfigs = 27
  character (len=500) :: config_names (nconfigs)
  character (len=500) :: config_values (nconfigs)
  character (len=500) :: site_list, output_file, output_file2, grid_list
  character (len=500) :: directory

  character (len=100) :: startdate, enddate ! desired output dates, YYYYMMDD, read from config file
  character (len=20)  :: stn_startdate, stn_enddate ! input station period dates (YYYYMMDD) from config file
  character (len=100) :: perturbation, station_var, site_var, site_var_t
  character (len=100), allocatable :: file_var (:), var_name (:)
  character (len=100), allocatable :: nwp_vars(:)    !list of nwp predictor variables
  character (len=100)         :: nwp_prcp_var = ''   !name of NWP precipitation predictor variable
  character (len=100), allocatable :: stnid (:), stnname (:)
  character (len=2), allocatable :: vars (:) !AWW-feb2016 for station P/T indicators

  character (len = 500) :: gen_sta_weights     ! flag for generating station weight file
  character (len = 500) :: sta_weight_name     ! name of station weight file
  character (len = 500) :: use_stn_weights     ! option for using station weights
 
  character (len=2000) :: arg !command line arg for configuration file
  character (len=2000) :: output_file_tmp !temporary output file name
  character (len=2000) :: sys_str !string for system commands
  character (len=2000) :: nwp_input_list !name of file containint list of NWP input files

  integer :: i, error, n_vars, nfile_var, nvar_name, forecast, mode
  integer :: nstations, lenfile
 
  real (dp), allocatable :: y (:, :, :), vals (:, :), lats (:), lons (:)
  real (dp), allocatable :: coefs (:, :, :), prob_coefs (:, :, :)
  real (dp), allocatable :: stnlat (:), stnlon (:), stnelev (:), stn_slp_n (:), stn_slp_e (:)
  real (dp), allocatable :: times (:)
  !modified AJN Sept 2013
  real (dp), allocatable :: mean_autocorr (:)!mean auto correlation for all stations over entire time period
  real (dp), allocatable :: mean_tp_corr (:)!mean correlation between precip and trange (31-day moving avg anomaly)
  real (dp), allocatable :: obs_max_pcp (:, :)
 
  !added for grid netcdf read  Oct 2015 AJN
  real (dp), allocatable :: lat (:, :)
  real (dp), allocatable :: lon (:, :)
  real (dp), allocatable :: elev (:, :)
  real (dp), allocatable :: grad_n (:, :)
  real (dp), allocatable :: grad_e (:, :)
  real (dp), allocatable :: mask (:, :)
 
  integer(I4B)           :: kfold_trials        !number of kfold xval trials
  integer(I4B)           :: kfold_hold          !number of stations to withhold from regression

  real (dp), allocatable :: x (:, :), z (:, :)
  integer :: ngrid
 
  integer (i4b) :: nx, ny, ntimes
  integer (i4b) :: nPredict           ! total number of predictors for spatial regression
  integer (i4b) :: n_nwp           ! number of NWP predictors for spatial regression
  integer (i4b) :: st_stndata_utime, end_stndata_utime, st_rec, end_rec ! AWW
 
  real (dp) :: maxdistance
  real (dp), allocatable :: grdlat (:), grdlon (:), grdelev (:), grd_slp_n (:), grd_slp_e (:), &
                            & mask_1d (:)
  real (sp), allocatable :: pcp (:, :), pop (:, :), pcperror (:, :)
  !modified AJN Sept 2013
  real (sp), allocatable :: tmean (:, :), tmean_err (:, :)
  real (sp), allocatable :: trange (:, :), trange_err (:, :)
 
  real (sp), allocatable :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
  !modified AJN Sept 2013
  real (sp), allocatable :: tmean_2 (:, :), tmean_err_2 (:, :)
  real (sp), allocatable :: trange_2 (:, :), trange_err_2 (:, :)
 
  ! ========== code starts below ==============================
 
  ! get config_file filename from command line
  i = 0
  do
    call get_command_argument (i, arg)
    if (i .eq. 1) config_file = arg
    if (len_trim(arg) == 0) exit
    i = i + 1
  end do

  ! initialize 
  error    = 0
  n_vars   = 0
  forecast = -1
 
  ! get configuration information, make assignments
  call read_config (config_file, nconfigs, config_names, config_values)
  call value (config_values(1), mode, error)
  if (error /= 0) then
    print *, "ERROR: Failed to read mode from config file."
    stop
  end if
 
  ! store rest of config file values (ordered)
  startdate       = config_values(2)
  enddate         = config_values(3)
  site_list       = config_values(4)
  site_var        = config_values(5)
  station_var     = config_values(6)
  output_file     = config_values(12)
  site_var_t      = config_values(15)
  directory       = config_values(16)
  stn_startdate   = config_values(17)
  stn_enddate     = config_values(18)
  gen_sta_weights = config_values(19)
  sta_weight_name = config_values(20)
  use_stn_weights = config_values(27)

  call value(config_values(21),nPredict,error)
  if (error /= 0) then
    print *, "ERROR: Failed to read npredict from config file."
    stop
  end if

  !check to see if output file path is valid
  !create the output file and see if an error occurs
  output_file_tmp = trim(output_file) // ".txt"
  open(unit=34,file=trim(output_file_tmp),form='unformatted',iostat=error)
    
  ! check to see if it was created
  if(error /= 0) then
    print *, "Error: Output path is not valid." 
    print *, trim(output_file_tmp), " cannot be created in output directory"
    stop
  else
    sys_str = "rm " // trim(output_file_tmp)
    call system(sys_str)
  end if
 
  ! print *,trim(site_var_t),' ',trim(site_var)
 
  ! check that dates were entered
  if (len(trim(startdate)) == 0 .or. len(trim(enddate)) == 0 .or. len(trim(site_list)) == 0 .or. &
    & len(trim(output_file)) == 0) then
    print *, "ERROR: Failed to read in one more more required config parameters. &
      & START_DATE, END_DATE, SITE_LIST or OUTPUT_FILE"
    stop
  end if
 
  ! --- AWW:  calculate start and end utimes & records for requested station data read period ---
  call get_time_list (startdate, enddate, times)  ! makes unix-time list for desired records
  ntimes = size (times)
  print *, 'startdate=', startdate, 'enddate=', enddate, 'ntimes=', ntimes  ! YYYYMMDD dates
  print *, "---------"

  ! translate times to a start and end record for station files

  ! NOTE: the station file start & end dates are given in the config file so that this calculation 
  !   doesn't have to be done for every single station (ie reading st/end from the station file
  !   but it would be more flexible to do it for every station file, perhaps -- then they would
  !   not all have to be the same lengths when input

  ! --- station data start and end utimes are now derived from config file dates
  st_stndata_utime  = date_to_unix (stn_startdate) ! returns secs-since-1970 for st date of station files
  end_stndata_utime = date_to_unix (stn_enddate)  ! returns secs-since-1970 for end date of station files

  ! calculate start & end recs for processing station data files 
  st_rec  = floor((times(1)-st_stndata_utime)/86400) + 1
  end_rec = floor((times(ntimes)-st_stndata_utime)/86400) + 1
  print *, 'st_rec and end_rec days since 1970/1/1: ', st_rec, end_rec
  ! --- end AWW add ---
 
 
  ! === CHOOSE BETWEEN MODE 1 and MODE 2 ====
  !   1) ens. source is gridded variables --- NOT USED/Maintained -- superceded by GARD software
  !   2) ens. source is station data (and optionally gridded predictors)

  if (mode == 1) then
  
    print *, 'Mode 1 (gridded downscaling) is not currently being maintained in this code'
    print *, 'The functionality has been recreated and expanded in the NCAR GARD software'
    stop

    ! =================== Ensemble Source Is Gridded Model Variables =================

    perturbation = config_values (7)
    call value (config_values(8), forecast, error)
    if (error /= 0) forecast = - 1
    call value (config_values(9), n_vars, error)
    if (error /= 0) n_vars = 0
 
    if (len(trim(perturbation)) == 0 .or. n_vars == 0 .or. forecast ==-1) then
      print *, "ERROR: Failed to read in one or more required model config parameters. &
               &(PERTURBATION, NUMBER_VARS or FORECAST)"
      stop
    end if
 
    if (n_vars > 1) then
      allocate (file_var(n_vars))
      allocate (var_name(n_vars))
      call parse (config_values(10), ",", file_var, nfile_var)
      call parse (config_values(11), ",", var_name, nvar_name)
      if (nfile_var /= n_vars .or. nvar_name /= n_vars) then
        print *, "ERROR: Number of variables in config file does not match."
        stop
      end if
    else
      allocate (file_var(1))
      allocate (var_name(1))
      file_var (1) = config_values (10)
      var_name (1) = config_values (11)
      if (len(trim(file_var(1))) == 0 .or. len(trim(var_name(1))) == 0) then
        print *, "ERROR: Failed to read in one or more required model config parameters. &
                 &(FILE_VARIABLE or VARIABLE_NAME)"
        stop
      end if
 
    end if
     !  print *,'start date ',startdate
 
    do i = 1, n_vars, 1
      call read_refcst (startdate, enddate, file_var(i), perturbation, var_name(i), forecast, vals, &
     & lats, lons, times, error)
      if (error /= 0) then
        print*, 'ERROR: subroutine read_refcst() returned an error: ', error
        stop
      end if 
      if (i == 1) then
        allocate (y(n_vars, size(lons)*size(lats), size(times)))
      end if
      y (i, :, :) = vals (:, :)
      deallocate (vals)
    end do
 
    call read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnelev, stn_slp_n, &
   & stn_slp_e, nstations, error, vars)  ! AWW-feb2016 handled reading station variables
    if (error /= 0) then
      print*, "ERROR in reading station list: ", error
      stop
    end if
 
    call estimate_coefficients (y, n_vars, lats, lons, times, st_rec, end_rec, stnid, stnlat, &
   & stnlon, station_var, directory, coefs, prob_coefs, error)
    if (error /= 0) then
      print*, "ERROR calling estimate_coefficients() routine: ", error
      stop
    end if
 
    if (trim(station_var) .eq. "PRCP") then
      lenfile = len_trim (output_file)
      output_file2 (:) = " "
      output_file2 (1:5) = "prob_"
      output_file2 (6:lenfile+6) = output_file (1:lenfile)
      print *, trim (output_file2)
 
      ! store PoP
      call save_coefficients (n_vars, var_name, prob_coefs, startdate, enddate, times, site_list, &
     & stnid, stnlat, stnlon, stnelev, forecast, output_file2, error)
      if (error /= 0) then
        print*, "ERROR calling save_coefficients() subroutine, pop: ", error
        stop
      end if 
    end if
 
    ! store PCP  ( should this be inside the if block above?)
    call save_coefficients (n_vars, var_name, coefs, startdate, enddate, times, site_list, &
   & stnid, stnlat, stnlon, stnelev, forecast, output_file, error)
    if (error /= 0) then
      print*, "ERROR calling save_coefficients() subroutine, pcp: ", error
      stop
    end if
 
    ! AWW: note Mode 1 has not been fully coded to use st_rec & end_rec ... just passed now because
    !   an internal subroutine wants them and may use them in the future
    !   also Mode 1 has not be used much with this program, and may be removed as it has been 
    !   superceded by GARD
    !   lastly some of the subroutines (eg save_coefficients & estimate_coefficients are 
    !   now duplicated by the save... and estimate... routines called in mode2
 
  else if (mode == 2) then
 
    ! =================== Ensemble Forcing Generation =====================
 
    call value (config_values(14), maxdistance, error)  ! convert config str to number
    print*, "Max Distance = ", maxdistance
    if (error /= 0) then
      print*, "Max Distance not correctly read ... quitting", error
      stop
    end if
    maxdistance = maxdistance * 0.539957   ! convert from km to nautical miles
 
    n_nwp = nPredict-6
    allocate(nwp_vars(n_nwp))
    call parse (config_values(22), ",", nwp_vars, nfile_var)
    if(nfile_var /= n_nwp) then
      print *,"Failed to read in NWP variable list or nPredict does not equal 6 + number of NWP predictors:", nfile_var
      stop
    end if

    nwp_prcp_var   = config_values(23)
    nwp_input_list = config_values(24)

    call read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnelev, stn_slp_n, &
   & stn_slp_e, nstations, error, vars) ! AWW added vars
    if (error /= 0) then
      print *, "ERROR: Failed to read station list ... quitting", error
      stop
    end if

    ! read grid domain file 
    grid_list = config_values (13)
    if (len(trim(grid_list)) == 0) then
      print *, "ERROR: Failed to read GRID_LIST (domain grid) name"
      stop
    end if

    call read_domain_grid (grid_list, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
    if(error /= 0) then
      print *, "ERROR: Failed to read domain grid ... quitting", error
      stop
    end if

    !grab kfold cross validation config parameters
    call value(config_values(25), kfold_trials, error)  ! convert config str to number
    call value(config_values(26), kfold_hold,   error)  ! convert config str to number

    !check kfold values
    !limit number of trial as 10-30
    if(kfold_trials .gt. 30 .or. kfold_trials .lt. 2) then
      print *,'Error:  K-fold trials limited to range of 10-30 trials'
      stop
    end if
    !limit number of stations withheld to 1-10
    if(kfold_hold .gt. 10 .or. kfold_hold .lt. 1) then
      print *,'Error:  K-fold station withholding limited to range of 1-10 stations withheld'
      stop
    end if

    ! allocate 1-d grid variables
    print*, "allocating vector variables matching grids of nx= ",nx," by ny= ",ny
    allocate (grdlat(nx*ny))
    allocate (grdlon(nx*ny))
    allocate (grdelev(nx*ny))
    allocate (grd_slp_n(nx*ny))
    allocate (grd_slp_e(nx*ny))
    allocate (mask_1d(nx*ny))

    ! reshape the domain grids into vectors
    grdlat    = reshape (lat, (/ nx*ny /))
    grdlon    = reshape (lon, (/ nx*ny /))
    grdelev    = reshape (elev, (/ nx*ny /))
    grd_slp_n = reshape (grad_n, (/ nx*ny /))
    grd_slp_e = reshape (grad_e, (/ nx*ny /))
    mask_1d   = reshape (mask, (/ nx*ny /))
 
    ngrid = nx * ny

    allocate (x(nstations, nPredict))   ! x arrays for station variables
    allocate (z(ngrid, nPredict))

    !define first six predictors as geophysical predictors
    x(:, 1) = 1.0
    x(:, 2) = stnlat (:)
    x(:, 3) = stnlon (:)
    x(:, 4) = stnelev (:)
    x(:, 5) = stn_slp_n (:)
    x(:, 6) = stn_slp_e (:)
 
    z(:, 1) = 1.0               ! z arrays for grid variables
    z(:, 2) = grdlat (:)
    z(:, 3) = grdlon (:)
    z(:, 4) = grdelev (:)
    z(:, 5) = grd_slp_n (:)
    z(:, 6) = grd_slp_e (:)
 
    allocate (mean_autocorr(ntimes))
    allocate (mean_tp_corr(ntimes))
    allocate (obs_max_pcp(ngrid, ntimes))
 
    call estimate_forcing_regression (nPredict, gen_sta_weights, sta_weight_name, nwp_input_list, &
   & n_nwp, nwp_vars, nwp_prcp_var, x, z, ngrid, maxdistance, times, st_rec, end_rec, &
   & stnid, station_var, directory, kfold_trials, kfold_hold, pcp, pop, pcperror, obs_max_pcp, tmean, &
   & tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, &
   & error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, trange_err_2, use_stn_weights)

    if (error /= 0) then
      print *, "ERROR: subroutine estimate_forcing_regression() returned error", error
      stop
    end if
 
    print *, 'Creating output file'
 
    call save_forcing_regression (pcp, pop, pcperror, obs_max_pcp, tmean, tmean_err, trange, trange_err, nx, &
   & ny, grdlat, grdlon, grdelev, times, mean_autocorr, mean_tp_corr, &
   & output_file, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, trange_err_2)

    if (error /= 0) then
      print *, "ERROR: subroutine save_forcing_regression() returned error", error
      stop
    end if
    
    ! ===== end Mode 2:  ensemble regression
 
  else

    ! mode not recognized...stop
    print*, 'Mode given in config file = ',mode,' Not recognized (can be 1 or 2).  Quitting.' 
    stop

  end if
  
end program gmet


! =============================================== 
! ================= SUBROUTINES =================
! =============================================== 
 
subroutine get_time_list (startdate, enddate, times)
  ! makes a list of data times in secs since 1970-1-1
  ! corresponding to requested period
  use utim
  use type
  implicit none
 
  character (len=100), intent (in) :: startdate, enddate
  real (dp), allocatable, intent (out) :: times (:)
  integer (i4b) :: t, ntimes, sday, eday, error
  integer (i4b) :: sec, min, hour, day, month, year
  real (dp) :: utime
 
  call parse_date (startdate, year, month, day, hour, min, sec, error)
  sday = julian_date (day, month, year)
  call parse_date (enddate, year, month, day, hour, min, sec, error)
  eday = julian_date (day, month, year)
  ntimes = eday - sday + 1
  allocate (times(ntimes))
 
  utime = date_to_unix (startdate)  ! secs since 1970-1-1
  print *, 'times!', eday, sday, utime, date_to_unix (enddate), ntimes
  do t = 1, ntimes, 1
    if (utime > date_to_unix(enddate)) exit
    times (t) = utime
    utime = utime + 86400
  end do
 
  print *, 'time list:', times  !seconds since 1970-1-1
 
end subroutine get_time_list
