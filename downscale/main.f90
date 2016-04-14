! Main program file for the Gridded Meteorological Ensemble Tool (GMET)
program gmet
  use type
  use strings
  use utim ! AWW
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
 
    subroutine read_station_list (file_name, id, name, lat, lon, alt, sslp_n, sslp_e, n_stations, &
   & error, vars) !AWW
      use type
      character (len=500), intent (in) :: file_name
      character (len=100), allocatable, intent (out) :: id (:), name (:)
      character (len=2), allocatable, intent (out) :: vars (:) !AWW holds PT identifiers
      real (dp), allocatable, intent (out) :: lat (:), lon (:), alt (:), sslp_n (:), sslp_e (:)
      integer (i4b), intent (out) :: n_stations
      integer, intent (out) :: error
    end subroutine read_station_list
 
    subroutine read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
      use type
      character (len=500), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_grid_list
 
    subroutine read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
      use netcdf
      use type
 
      character (len=500), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), &
     & grad_e (:, :), mask (:, :)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_nc_grid
 
    ! subroutine estimate_coefficients(D, nvars, Lats, Lons, Times, stnid, stnlat, stnlon, &
    ! AWW modified feb-2016
    subroutine estimate_coefficients (d, nvars, lats, lons, times, st_rec, end_rec, stnid, stnlat, &
   & stnlon, stnalt, stnvar, site_var, site_list, directory, c, poc, error)
      use type
      real (dp), intent (in) :: d (:, :, :), lats (:), lons (:)
      real (dp), intent (in) :: times (:)
      integer (i4b), intent (in) :: st_rec, end_rec ! AWW
      integer (i4b), intent (in) :: nvars
      character (len=100), intent (in) :: stnid (:)
      real (dp), intent (in) :: stnlat (:), stnlon (:), stnalt (:)
      character (len=100), intent (in) :: stnvar
      character (len=100), intent (in) :: site_var
      character (len=500), intent (in) :: site_list
      character (len=500), intent (in) :: directory !AWW added
      real (dp), allocatable, intent (out) :: c (:, :, :), poc (:, :, :)
      integer, intent (out) :: error
    end subroutine estimate_coefficients
 
    subroutine save_coefficients (n_vars, var_names, coefs, startdate, enddate, times, site_list, &
   & station_var, stnid, stnlat, stnlon, stnalt, forecast, file, error)
      use netcdf
      use type
      character (len=100), intent (in) :: var_names (:), stnid (:)
      integer, intent (in) :: n_vars, forecast
      character (len=100), intent (in) :: startdate, enddate, station_var
      character (len=500), intent (in) :: file, site_list
      real (dp), intent (in) :: stnlat (:), stnlon (:), stnalt (:)
      real (dp), intent (in) :: coefs (:, :, :)
      real (dp), intent (in) :: times (:)
      integer, intent (out) :: error
    end subroutine save_coefficients
 
    !modified AJN Sept 2013
    ! subroutine estimate_precip(X, Z, nsta, ngrid, maxDistance, Times,  &
    ! subroutine estimate_precip(X, Z, nsta, ngrid, maxDistance, Times, st_rec, end_rec, &
    subroutine estimate_forcing_regression (x, z, nsta, ngrid, maxdistance, times, st_rec, end_rec, &
   & stnid, stnvar, site_var, site_var_t, site_list, directory, pcp, pop, pcperr, tmean, tmean_err, trange, &
   & trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, error, pcp_2, &
   & pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
      use type
      real (dp), intent (in) :: x (:, :), z (:, :)
      real (dp), intent (in) :: maxdistance
      integer (i4b), intent (in) :: nsta, ngrid
      real (dp), intent (in) :: times (:)
      integer (i4b) :: st_stndata_utime, end_stndata_utime, st_rec, end_rec ! AWW added
      character (len=100), intent (in) :: stnid (:)
      character (len=100), intent (in) :: stnvar, site_var, site_var_t
      character (len=500), intent (in) :: site_list
      character (len=500), intent (in) :: directory !AWW
      real (sp), allocatable, intent (out) :: pcp (:, :), pop (:, :), pcperr (:, :)
      real (sp), allocatable, intent (out) :: tmean (:, :), tmean_err (:, :)!OLS tmean estimate and error
      real (sp), allocatable, intent (out) :: trange (:, :), trange_err (:, :)!OLS trange estimate and error
      real (sp), allocatable, intent (out) :: pcp_2 (:, :), pop_2 (:, :), pcperr_2 (:, :)
      real (sp), allocatable, intent (out) :: tmean_2 (:, :), tmean_err_2 (:, :)!OLS tmean estimate and error
      real (sp), allocatable, intent (out) :: trange_2 (:, :), trange_err_2 (:, :)!OLS trange estimate and error
      integer, intent (out) :: error
      real (dp), intent (out) :: mean_autocorr (:)!mean autocorrelation from all stations over entire time period
      real (dp), intent (out) :: mean_tp_corr (:)!mean correlation for mean temp and precip
      real (dp), intent (out) :: y_mean (:, :), y_std (:, :), y_std_all (:, :)!std and mean of transformed time step precipitation
      real (dp), intent (out) :: y_min (:, :), y_max (:, :)!min,max of normalized time step precip
    end subroutine estimate_forcing_regression
    !end subroutine estimate_precip  ! AWW-Feb2016 renamed to be more descriptive
 
    !modified AJN Sept 2013
    !subroutine save_precip(pcp, pop, pcperror, tmean, tmean_err, trange, trange_err, &
    ! AWW Feb2016, renamed to be more descriptive.  Note did NOT call scrf/save_precip.f90 subroutine
    subroutine save_forcing_regression (pcp, pop, pcperror, tmean, tmean_err, trange, trange_err, &
   & nx, ny, grdlat, grdlon, grdalt, times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, &
   & y_min, y_max, file, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, &
   & trange_err_2)
      use netcdf
      use type
      real (sp), intent (in) :: pcp (:, :), pop (:, :), pcperror (:, :)
      real (sp), intent (in) :: tmean (:, :), tmean_err (:, :), trange (:, :), trange_err (:, :)
      real (sp), intent (in) :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
      real (sp), intent (in) :: tmean_2 (:, :), tmean_err_2 (:, :), trange_2 (:, :), trange_err_2 &
     & (:, :)
      integer (i4b), intent (in) :: nx, ny
      real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
      real (dp), intent (in) :: times (:)
      real (dp), intent (in) :: mean_autocorr (:), mean_tp_corr (:)
      real (dp), intent (in) :: y_mean (:, :), y_std (:, :), y_std_all (:, :)
      real (dp), intent (in) :: y_min (:, :), y_max (:, :)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_forcing_regression
    !end subroutine save_precip
 
  end interface
  ! === end of interfaces, start of program ====
 
  character (len=100) :: config_file
  !integer, parameter :: nconfigs = 15 !modified AJN Sept 2013
  integer, parameter :: nconfigs = 18 !modified AWW Feb 2016
  character (len=500) :: config_names (nconfigs)
  character (len=500) :: config_values (nconfigs)
  character (len=500) :: site_list, output_file, output_file2, grid_list
  character (len=500) :: directory ! AWW feb 2016

  character (len=100) :: startdate, enddate ! desired output dates, YYYYMMDD, read from config file
  character (len=20) :: stn_startdate, stn_enddate ! input station period dates (YYYYMMDD) from config file
  character (len=100) :: perturbation, station_var, site_var, site_var_t
  character (len=100), allocatable :: file_var (:), var_name (:)
  character (len=100), allocatable :: stnid (:), stnname (:)
  character (len=2), allocatable :: vars (:) !AWW-feb2016 for station P/T indicators
 
  character (len=2000) :: arg !command line arg for configuration file
 
  integer :: i, error, n_vars, nfile_var, nvar_name, forecast, mode
  integer :: nstations, lenfile
 
  real (dp), allocatable :: y (:, :, :), vals (:, :), lats (:), lons (:)
  real (dp), allocatable :: coefs (:, :, :), prob_coefs (:, :, :)
  real (dp), allocatable :: stnlat (:), stnlon (:), stnalt (:), stn_slp_n (:), stn_slp_e (:)
  real (dp), allocatable :: times (:)
  !modified AJN Sept 2013
  real (dp), allocatable :: mean_autocorr (:)!mean auto correlation for all stations over entire time period
  real (dp), allocatable :: mean_tp_corr (:)!mean correlation between precip and trange (31-day moving avg anomaly)
  real (dp), allocatable :: y_mean (:, :)
  real (dp), allocatable :: y_std (:, :)
  real (dp), allocatable :: y_std_all (:, :)
  real (dp), allocatable :: y_max (:, :)
  real (dp), allocatable :: y_min (:, :)
 
  !added for grid netcdf read  Oct 2015 AJN
  real (dp), allocatable :: lat (:, :)
  real (dp), allocatable :: lon (:, :)
  real (dp), allocatable :: elev (:, :)
  real (dp), allocatable :: grad_n (:, :)
  real (dp), allocatable :: grad_e (:, :)
  real (dp), allocatable :: mask (:, :)
 
 
  real (dp), allocatable :: x (:, :), z (:, :)
  integer :: ngrid
 
  integer (i4b) :: nx, ny, ntimes
 
  integer (i4b) :: st_stndata_utime, end_stndata_utime, st_rec, end_rec ! AWW
 
  real (dp) :: maxdistance
  real (dp), allocatable :: grdlat (:), grdlon (:), grdalt (:), grd_slp_n (:), grd_slp_e (:), &
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
 
  ! test configs
  ! mode 1
  !  config_file = "config_prcp.txt" 
  ! mode 2
  !  config_file = "config_pnw.txt"
  
  ! get config_file filename from command line
  i = 0
  do
    call get_command_argument (i, arg)
    if (i .eq. 1) config_file = arg
    if (len_trim(arg) == 0) exit
    i = i + 1
  end do
 
  config_names(1) = "MODE"
  config_names(2) = "START_DATE"
  config_names(3) = "END_DATE"
  config_names(4) = "SITE_LIST"
  config_names(5) = "SITE_VAR"
  config_names(6) = "STATION_VAR"
  config_names(7) = "PERTURBATION"
  config_names(8) = "FORECAST"
  config_names(9) = "NUMBER_VARS"
  config_names(10) = "FILE_VARIABLE"
  config_names(11) = "VARIABLE_NAME"
  config_names(12) = "OUTPUT_FILE"
  config_names(13) = "GRID_LIST"
  config_names(14) = "MAX_DISTANCE"
  config_names(15) = "SITE_VAR_T" !modified AJN Sept 2013
  config_names(16) = "DATA_DIRECTORY" !AWW-feb2016, free data dir from site list path
  config_names(17) = "STN_START_DATE" !AWW-apr2016, add station period limits
  config_names(18) = "STN_END_DATE"   !  to get them out of being hardwired
 
  error    = 0
  n_vars   = 0
  forecast = -1
 
  call read_config (config_file, nconfigs, config_names, config_values)
  call value (config_values(1), mode, error)
  if (error /= 0) then
    print *, "ERROR: Failed to read mode from config file."
    return
  end if
 
  startdate   = config_values(2)
  enddate     = config_values(3)
  site_list   = config_values(4)
  site_var    = config_values(5)
  site_var_t  = config_values(15)
  station_var = config_values(6)
  output_file = config_values(12)
  directory   = config_values(16)
  stn_startdate   = config_values(17)
  stn_enddate     = config_values(18)
 
  ! print *,trim(site_var_t),' ',trim(site_var)
 
  if (len(trim(startdate)) == 0 .or. len(trim(enddate)) == 0 .or. len(trim(site_list)) == 0 .or. &
 & len(trim(output_file)) == 0) then
    print *, "ERROR: Failed to read in one more more required config parameters. (START_DATE, END_D&
   &ATE, SITE_LIST or OUTPUT_FILE)"
    return
  end if
 
  ! --- AWW:  figure out start and end utimes & records for requested station data read period ---
  call get_time_list (startdate, enddate, times)  ! makes unix-time list for desired records
  ntimes = size (times)
  print *, 'startdate=', startdate, 'enddate=', enddate, 'ntimes=', ntimes  ! YYYYMMDD dates

  !translate times to a start and end record for station files

  ! the station file start & end dates are given in the config file so that this calculation 
  !   doesn't have to be done for every single station (ie reading st/end from the station file
  !   but it would be more flexible to do it for every station file, perhaps -- then they would
  !   not all have to be the same lengths

  ! retro dates
  !st_stndata_utime = date_to_unix ('19800101')  ! returns secs-since-1970 for st date of station files
                                                ! hardwired for testing
  !end_stndata_utime = date_to_unix ('20141231') ! returns secs-since-1970 for end date of station files
                                                ! hardwired for testing
  ! spinup dates
  !st_stndata_utime = date_to_unix ('20150101')  ! returns secs-since-1970 for st date of station files
                                                 ! hardwired for testing
  !end_stndata_utime = date_to_unix ('20160430') ! returns secs-since-1970 for end date of station files
                                                 ! hardwired for testing
  ! --- station data start and end dates
  st_stndata_utime = date_to_unix (stn_startdate) ! returns secs-since-1970 for st date of station files
  end_stndata_utime = date_to_unix (stn_enddate)  ! returns secs-since-1970 for end date of station files

  ! calculate start & end recs for processing station data files 
  st_rec = floor ((times(1)-st_stndata_utime)/86400) + 1
  end_rec = floor ((times(ntimes)-st_stndata_utime)/86400) + 1
 
  print *, 'st_rec and end_rec =', st_rec, end_rec
  ! --- end AWW add ---
 
 
  ! === CHOOSE BETWEEN MODE 1 and MODE 2 ====
  if (mode == 1) then

    ! =================== Ensemble Source Is Gridded Model Variables =================

    perturbation = config_values (7)
    call value (config_values(8), forecast, error)
    if (error /= 0) forecast = - 1
    call value (config_values(9), n_vars, error)
    if (error /= 0) n_vars = 0
 
    if (len(trim(perturbation)) == 0 .or. n_vars == 0 .or. forecast ==-1) then
      print *, "ERROR: Failed to read in one or more required model config parameters. (PERTURBATIO&
     &N, NUMBER_VARS or FORECAST)"
      return
    end if
 
    if (n_vars > 1) then
      allocate (file_var(n_vars))
      allocate (var_name(n_vars))
      call parse (config_values(10), ",", file_var, nfile_var)
      call parse (config_values(11), ",", var_name, nvar_name)
      if (nfile_var /= n_vars .or. nvar_name /= n_vars) then
        print *, "ERROR: Number of variables in config file does not match."
        return
      end if
    else
      allocate (file_var(1))
      allocate (var_name(1))
      file_var (1) = config_values (10)
      var_name (1) = config_values (11)
      if (len(trim(file_var(1))) == 0 .or. len(trim(var_name(1))) == 0) then
        print *, "ERROR: Failed to read in one or more required model config parameters. (FILE_VARI&
       &ABLE or VARIABLE_NAME)"
        return
      end if
 
    end if
     !  print *,'start date ',startdate
 
    do i = 1, n_vars, 1
      call read_refcst (startdate, enddate, file_var(i), perturbation, var_name(i), forecast, vals, &
     & lats, lons, times, error)
      if (error /= 0) return
      if (i == 1) then
        allocate (y(n_vars, size(lons)*size(lats), size(times)))
      end if
      y (i, :, :) = vals (:, :)
      deallocate (vals)
    end do
 
    call read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, &
   & stn_slp_e, nstations, error, vars)  ! AWW-feb2016 handled reading station variables
    !   & stn_slp_e, nstations, error)
    if (error /= 0) return
 
    !call estimate_coefficients(Y, n_vars, Lats, Lons, Times, stnid, stnlat, stnlon, &
    call estimate_coefficients (y, n_vars, lats, lons, times, st_rec, end_rec, stnid, stnlat, &
   & stnlon, stnalt, station_var, site_var, site_list, directory, coefs, prob_coefs, error)
    if (error /= 0) return
 
    if (trim(station_var) .eq. "PRCP") then
      lenfile = len_trim (output_file)
      output_file2 (:) = " "
      output_file2 (1:5) = "prob_"
      output_file2 (6:lenfile+6) = output_file (1:lenfile)
      print *, trim (output_file2)
 
      call save_coefficients (n_vars, var_name, prob_coefs, startdate, enddate, times, site_list, &
     & station_var, stnid, stnlat, stnlon, stnalt, forecast, output_file2, error)
      if (error /= 0) return
    end if
 
    call save_coefficients (n_vars, var_name, coefs, startdate, enddate, times, site_list, &
   & station_var, stnid, stnlat, stnlon, stnalt, forecast, output_file, error)
    if (error /= 0) return
 
    ! AWW: note Mode 1 has not been fully coded to use st_rec & end_rec...just passed now because
    !   an internal subroutine wants them and may use them in the future
 
  else if (mode == 2) then
 
    ! =================== ensemble forcing generation =====================
 
    grid_list = config_values (13)
    call value (config_values(14), maxdistance, error)
    maxdistance = maxdistance * 0.539957   ! AWW...why?
    if (error /= 0) maxdistance = - 1
 
    if (len(trim(grid_list)) == 0) then
      print *, "ERROR: Failed to read in one more more required model config parameters.&
     & (GRID_LIST)"
      return
    end if
  
    call read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, &
   & stn_slp_e, nstations, error, vars) !AWW added vars
    if (error /= 0) return
 
    ! call read_grid_list(grid_list, grdlat, grdlon, grdalt, grd_slp_n, grd_slp_e, nx, ny, error)
    ! if(error /= 0) return
 
    call read_nc_grid (grid_list, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
 
    print *, size (lat), size (elev), nx, ny, nx * ny
 
    !allocate 1-d grid variables
    allocate (grdlat(nx*ny))
    allocate (grdlon(nx*ny))
    allocate (grdalt(nx*ny))
    allocate (grd_slp_n(nx*ny))
    allocate (grd_slp_e(nx*ny))
    allocate (mask_1d(nx*ny))
 
    print *, size (grdlat)
 
    grdlat = reshape (lat, (/ nx*ny /))
    grdlon = reshape (lon, (/ nx*ny /))
    grdalt = reshape (elev, (/ nx*ny /))
    grd_slp_n = reshape (grad_n, (/ nx*ny /))
    grd_slp_e = reshape (grad_e, (/ nx*ny /))
    mask_1d = reshape (mask, (/ nx*ny /))
 
    print *, size (lat), size (elev), nx, ny
 
    ngrid = nx * ny
    allocate (x(nstations, 6))   ! x arrays for station variables
    allocate (z(ngrid, 6))
    x(:, 1) = 1.0
    x(:, 2) = stnlat (:)
    x(:, 3) = stnlon (:)
    x(:, 4) = stnalt (:)
    x(:, 5) = stn_slp_n (:)
    x(:, 6) = stn_slp_e (:)
 
    z(:, 1) = 1.0               ! z arrays for grid variables
    z(:, 2) = grdlat (:)
    z(:, 3) = grdlon (:)
    z(:, 4) = grdalt (:)
    z(:, 5) = grd_slp_n (:)
    z(:, 6) = grd_slp_e (:)
 
    print *, size (lat), size (elev), nx, ny
 
    !       call get_time_list(startdate, enddate, Times)
    !	ntimes = size(Times)
    !        print *,'startdate=',startdate,'enddate=',enddate,'ntimes=',ntimes
  
    ! -- AWW:  translate times to a start and end record for station files
    !        st_stndata_utime = date_to_unix('19800101')   ! returns secs-since-1970 for enddate of stn files
                                                       ! hardwired for testing
    !        end_stndata_utime = date_to_unix('20141231')   ! returns secs-since-1970 for enddate of stn files
                                                       ! hardwired for testing
    !        st_rec  = FLOOR((Times(1) - st_stndata_utime)/86400) + 1
    !        end_rec = FLOOR((Times(ntimes) - st_stndata_utime)/86400) + 1
    !        print*, 'st_rec and end_rec =', st_rec, end_rec
    ! -- AWW:  end addition
 
    ! modified AJN Sept 2013
    allocate (mean_autocorr(ntimes))
    allocate (mean_tp_corr(ntimes))
    allocate (y_mean(ngrid, ntimes))
    allocate (y_std(ngrid, ntimes))
    allocate (y_std_all(ngrid, ntimes))
    allocate (y_max(ngrid, ntimes))
    allocate (y_min(ngrid, ntimes))
 
    ! call estimate_precip(X, Z, nstations, ngrid, maxDistance, Times, st_rec, end_rec, &
    call estimate_forcing_regression (x, z, nstations, ngrid, maxdistance, times, st_rec, &
   & end_rec, stnid, station_var, site_var, site_var_t, site_list, directory, pcp, pop, pcperror, tmean, &
   & tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, &
   & y_max, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
    if (error /= 0) return
 
    print *, 'Creating output file'
 
    ! call save_precip(pcp, pop, pcperror, tmean, tmean_err, trange, trange_err, &
    call save_forcing_regression (pcp, pop, pcperror, tmean, tmean_err, trange, trange_err, nx, &
   & ny, grdlat, grdlon, grdalt, times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, &
   & y_min, y_max, output_file, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, &
   & trange_err_2)
    if (error /= 0) return
 
  else

    ! mode not recognized...stop
    print*, 'Mode given in config file = ',mode,' Not recognized (can be 1 or 2).  Quitting.' 
    stop

  end if  ! end if mode=2
  
end program gmet
 
! ================= SUBROUTINES =================
 
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
