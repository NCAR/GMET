 
PROGRAM precip
  USE type
  USE strings
  IMPLICIT NONE
 
  INTERFACE
   SUBROUTINE get_time_list (startdate, enddate, T)
    USE utim
    USE type
    CHARACTER (LEN=100), INTENT (IN) :: startdate, enddate
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: T (:)
   END SUBROUTINE get_time_list
 
   SUBROUTINE read_refcst (startdate, enddate, file_var, perturbation, var_name, forecast, V, X, Y, T, error)
    USE type
    CHARACTER (LEN=100), INTENT (IN) :: startdate, enddate, file_var, perturbation
    CHARACTER (LEN=*), INTENT (IN) :: var_name
    INTEGER (I4B), INTENT (IN) :: forecast
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: V (:, :), X (:), Y (:)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: T (:)
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_refcst
 
   SUBROUTINE read_station_list (file_name, id, name, lat, lon, alt, sslp_n, sslp_e, n_stations, error)
    USE type
    CHARACTER (LEN=500), INTENT (IN) :: file_name
    CHARACTER (LEN=100), ALLOCATABLE, INTENT (OUT) :: id (:), name (:)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: lat (:), lon (:), alt (:), sslp_n (:), sslp_e (:)
    INTEGER (I4B), INTENT (OUT) :: n_stations
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_station_list
 
   SUBROUTINE read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
    USE type
    CHARACTER (LEN=500), INTENT (IN) :: file_name
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
    INTEGER (I4B), INTENT (OUT) :: nx, ny
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_grid_list
 
   SUBROUTINE read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
    USE netcdf
    USE type
 
    CHARACTER (LEN=500), INTENT (IN) :: file_name
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), grad_e (:, :), mask (:, &
   & :)
    INTEGER (I4B), INTENT (OUT) :: nx, ny
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_nc_grid
 
 
   SUBROUTINE estimate_coefficients (D, nvars, lats, lons, Times, stnid, stnlat, stnlon, stnalt, stnvar, site_var, &
  & site_list, C, POC, error)
    USE type
    REAL (DP), INTENT (IN) :: D (:, :, :), lats (:), lons (:)
    REAL (DP), INTENT (IN) :: Times (:)
    INTEGER (I4B), INTENT (IN) :: nvars
    CHARACTER (LEN=100), INTENT (IN) :: stnid (:)
    REAL (DP), INTENT (IN) :: stnlat (:), stnlon (:), stnalt (:)
    CHARACTER (LEN=100), INTENT (IN) :: stnvar
    CHARACTER (LEN=100), INTENT (IN) :: site_var
    CHARACTER (LEN=500), INTENT (IN) :: site_list
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: C (:, :, :), POC (:, :, :)
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE estimate_coefficients
 
   SUBROUTINE save_coefficients (n_vars, var_names, coefs, startdate, enddate, Times, site_list, station_var, stnid, &
  & stnlat, stnlon, stnalt, forecast, file, error)
    USE netcdf
    USE type
    CHARACTER (LEN=100), INTENT (IN) :: var_names (:), stnid (:)
    INTEGER, INTENT (IN) :: n_vars, forecast
    CHARACTER (LEN=100), INTENT (IN) :: startdate, enddate, station_var
    CHARACTER (LEN=500), INTENT (IN) :: file, site_list
    REAL (DP), INTENT (IN) :: stnlat (:), stnlon (:), stnalt (:)
    REAL (DP), INTENT (IN) :: coefs (:, :, :)
    REAL (DP), INTENT (IN) :: Times (:)
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE save_coefficients
 
   SUBROUTINE estimate_precip (X, Z, nsta, ngrid, maxDistance, Times, stnid, stnvar, site_var, site_var_t, site_list, &
  & PCP, POP, PCPERR, tmean, tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, &
  & y_min, y_max, error, pcp_2, pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
    USE type
    REAL (DP), INTENT (IN) :: X (:, :), Z (:, :)
    REAL (DP), INTENT (IN) :: maxDistance
    INTEGER (I4B), INTENT (IN) :: nsta, ngrid
    REAL (DP), INTENT (IN) :: Times (:)
    CHARACTER (LEN=100), INTENT (IN) :: stnid (:)
    CHARACTER (LEN=100), INTENT (IN) :: stnvar, site_var, site_var_t
    CHARACTER (LEN=500), INTENT (IN) :: site_list
    REAL (SP), ALLOCATABLE, INTENT (OUT) :: PCP (:, :), POP (:, :), PCPERR (:, :)
    REAL (SP), ALLOCATABLE, INTENT (OUT) :: tmean (:, :), tmean_err (:, :)!OLS tmean estimate and error
    REAL (SP), ALLOCATABLE, INTENT (OUT) :: trange (:, :), trange_err (:, :)!OLS trange estimate and error
 
    REAL (SP), ALLOCATABLE, INTENT (OUT) :: pcp_2 (:, :), pop_2 (:, :), pcperr_2 (:, :)
    REAL (SP), ALLOCATABLE, INTENT (OUT) :: tmean_2 (:, :), tmean_err_2 (:, :)!OLS tmean estimate and error
    REAL (SP), ALLOCATABLE, INTENT (OUT) :: trange_2 (:, :), trange_err_2 (:, :)!OLS trange estimate and error
 
 
    INTEGER, INTENT (OUT) :: error
    REAL (DP), INTENT (OUT) :: mean_autocorr (:)!mean autocorrelation from all stations over entire time period
    REAL (DP), INTENT (OUT) :: mean_tp_corr (:)!mean correlation for mean temp and precip
 
    REAL (DP), INTENT (OUT) :: y_mean (:, :), y_std (:, :), y_std_all (:, :)!std and mean of transformed time step precipitation
    REAL (DP), INTENT (OUT) :: y_min (:, :), y_max (:, :)!min,max of normalized time step precip
   END SUBROUTINE estimate_precip
 
   SUBROUTINE save_precip (PCP, POP, pcperror, tmean, tmean_err, trange, trange_err, nx, ny, grdlat, grdlon, grdalt, &
  & Times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, file, error, pcp_2, pop_2, pcperror_2, &
  & tmean_2, tmean_err_2, trange_2, trange_err_2)
    USE netcdf
    USE type
 
    REAL (SP), INTENT (IN) :: PCP (:, :), POP (:, :), pcperror (:, :)
    REAL (SP), INTENT (IN) :: tmean (:, :), tmean_err (:, :), trange (:, :), trange_err (:, :)
 
    REAL (SP), INTENT (IN) :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
    REAL (SP), INTENT (IN) :: tmean_2 (:, :), tmean_err_2 (:, :), trange_2 (:, :), trange_err_2 (:, :)
 
    INTEGER (I4B), INTENT (IN) :: nx, ny
    REAL (DP), INTENT (IN) :: grdlat (:), grdlon (:), grdalt (:)
    REAL (DP), INTENT (IN) :: Times (:)
    REAL (DP), INTENT (IN) :: mean_autocorr (:), mean_tp_corr (:)
 
    REAL (DP), INTENT (IN) :: y_mean (:, :), y_std (:, :), y_std_all (:, :)
    REAL (DP), INTENT (IN) :: y_min (:, :), y_max (:, :)
    CHARACTER (LEN=500), INTENT (IN) :: file
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE save_precip
 
  END INTERFACE
 
  CHARACTER (LEN=100) :: config_file
  INTEGER, PARAMETER :: nconfigs = 15
  CHARACTER (LEN=500) :: config_names (nconfigs)
  CHARACTER (LEN=500) :: config_values (nconfigs)
  CHARACTER (LEN=500) :: site_list, output_file, output_file2, grid_list
  CHARACTER (LEN=100) :: startdate, enddate, perturbation, station_var, site_var, site_var_t
  CHARACTER (LEN=100), ALLOCATABLE :: file_var (:), var_name (:)
  CHARACTER (LEN=100), ALLOCATABLE :: stnid (:), stnname (:)
 
  CHARACTER (LEN=2000) :: arg !command line arg for configuration file
 
  INTEGER :: i, error, n_vars, nfile_var, nvar_name, forecast, mode
  INTEGER :: nstations, lenfile
 
  REAL (DP), ALLOCATABLE :: Y (:, :, :), Vals (:, :), lats (:), lons (:)
  REAL (DP), ALLOCATABLE :: coefs (:, :, :), prob_coefs (:, :, :)
  REAL (DP), ALLOCATABLE :: stnlat (:), stnlon (:), stnalt (:), stn_slp_n (:), stn_slp_e (:)
  REAL (DP), ALLOCATABLE :: Times (:)
 
  REAL (DP), ALLOCATABLE :: mean_autocorr (:)!mean auto correlation for all stations over entire time period
  REAL (DP), ALLOCATABLE :: mean_tp_corr (:)!mean correlation between precip and trange (31-day moving avg anomaly)
  REAL (DP), ALLOCATABLE :: y_mean (:, :)
  REAL (DP), ALLOCATABLE :: y_std (:, :)
  REAL (DP), ALLOCATABLE :: y_std_all (:, :)
  REAL (DP), ALLOCATABLE :: y_max (:, :)
  REAL (DP), ALLOCATABLE :: y_min (:, :)
 
  REAL (DP), ALLOCATABLE :: lat (:, :)
  REAL (DP), ALLOCATABLE :: lon (:, :)
  REAL (DP), ALLOCATABLE :: elev (:, :)
  REAL (DP), ALLOCATABLE :: grad_n (:, :)
  REAL (DP), ALLOCATABLE :: grad_e (:, :)
  REAL (DP), ALLOCATABLE :: mask (:, :)
 
 
  REAL (DP), ALLOCATABLE :: X (:, :), Z (:, :)
  INTEGER :: ngrid
 
  INTEGER (I4B) :: nx, ny, ntimes
  REAL (DP) :: maxDistance
  REAL (DP), ALLOCATABLE :: grdlat (:), grdlon (:), grdalt (:), grd_slp_n (:), grd_slp_e (:), mask_1d (:)
  REAL (SP), ALLOCATABLE :: PCP (:, :), POP (:, :), pcperror (:, :)
  REAL (SP), ALLOCATABLE :: tmean (:, :), tmean_err (:, :)
  REAL (SP), ALLOCATABLE :: trange (:, :), trange_err (:, :)
 
  REAL (SP), ALLOCATABLE :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
  REAL (SP), ALLOCATABLE :: tmean_2 (:, :), tmean_err_2 (:, :)
  REAL (SP), ALLOCATABLE :: trange_2 (:, :), trange_err_2 (:, :)
 
!code starts below
 
!mode 2
  config_file = "config_pnw.txt"
!mode 1
!  config_file = "config_prcp.txt"
 
 
!get config_file filename from command line
  i = 0
  DO
   CALL get_command_argument (i, arg)
   IF (i .EQ. 1) config_file = arg
   IF (LEN_TRIM(arg) == 0) EXIT
   i = i + 1
  END DO
 
  config_names (1) = "MODE"
  config_names (2) = "START_DATE"
  config_names (3) = "END_DATE"
  config_names (4) = "SITE_LIST"
  config_names (5) = "SITE_VAR"
  config_names (6) = "STATION_VAR"
  config_names (7) = "PERTURBATION"
  config_names (8) = "FORECAST"
  config_names (9) = "NUMBER_VARS"
  config_names (10) = "FILE_VARIABLE"
  config_names (11) = "VARIABLE_NAME"
  config_names (12) = "OUTPUT_FILE"
  config_names (13) = "GRID_LIST"
  config_names (14) = "MAX_DISTANCE"
  config_names (15) = "SITE_VAR_T"
 
  error = 0
  n_vars = 0
  forecast = - 1
 
  CALL read_config (config_file, nconfigs, config_names, config_values)
  CALL value (config_values(1), mode, error)
  IF (error /= 0) THEN
   PRINT *, "ERROR: Failed to read mode from config file."
   RETURN
  END IF
 
  startdate = config_values (2)
  enddate = config_values (3)
  site_list = config_values (4)
  site_var = config_values (5)
  site_var_t = config_values (15)
  station_var = config_values (6)
  output_file = config_values (12)
 
 
  IF (len(trim(startdate)) == 0 .OR. len(trim(enddate)) == 0 .OR. len(trim(site_list)) == 0 .OR. len(trim(output_file)) &
 & == 0) THEN
   PRINT *, "ERROR: Failed to read in one more more required config parameters. (START_DATE, END_DATE, SITE_LIST or OUT&
  &PUT_FILE)"
   RETURN
  END IF
 
  IF (mode == 1) THEN
     ! Model Variables
   perturbation = config_values (7)
   CALL value (config_values(8), forecast, error)
   IF (error /= 0) forecast = - 1
   CALL value (config_values(9), n_vars, error)
   IF (error /= 0) n_vars = 0
 
   IF (len(trim(perturbation)) == 0 .OR. n_vars == 0 .OR. forecast ==-1) THEN
    PRINT *, "ERROR: Failed to read in one or more required model config parameters. (PERTURBATION, NUMBER_VARS or FORE&
   &CAST)"
    RETURN
   END IF
 
   IF (n_vars > 1) THEN
    ALLOCATE (file_var(n_vars))
    ALLOCATE (var_name(n_vars))
    CALL parse (config_values(10), ",", file_var, nfile_var)
    CALL parse (config_values(11), ",", var_name, nvar_name)
    IF (nfile_var /= n_vars .OR. nvar_name /= n_vars) THEN
     PRINT *, "ERROR: Number of variables in config file does not match."
     RETURN
    END IF
   ELSE
    ALLOCATE (file_var(1))
    ALLOCATE (var_name(1))
    file_var (1) = config_values (10)
    var_name (1) = config_values (11)
    IF (len(trim(file_var(1))) == 0 .OR. len(trim(var_name(1))) == 0) THEN
     PRINT *, "ERROR: Failed to read in one or more required model config parameters. (FILE_VARIABLE or VARIABLE_NAME)"
     RETURN
    END IF
 
   END IF
 
   DO i = 1, n_vars, 1
    CALL read_refcst (startdate, enddate, file_var(i), perturbation, var_name(i), forecast, Vals, lats, lons, Times, &
   & error)
    IF (error /= 0) RETURN
    IF (i == 1) THEN
     ALLOCATE (Y(n_vars, size(lons)*size(lats), size(Times)))
    END IF
    Y (i, :, :) = Vals (:, :)
    DEALLOCATE (Vals)
   END DO
 
   CALL read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, stn_slp_e, nstations, error)
   IF (error /= 0) RETURN
 
   CALL estimate_coefficients (Y, n_vars, lats, lons, Times, stnid, stnlat, stnlon, stnalt, station_var, site_var, &
  & site_list, coefs, prob_coefs, error)
   IF (error /= 0) RETURN
 
 
   IF (trim(station_var) .EQ. "PRCP") THEN
    lenfile = LEN_TRIM (output_file)
    output_file2 (:) = " "
    output_file2 (1:5) = "prob_"
    output_file2 (6:lenfile+6) = output_file (1:lenfile)
 
    CALL save_coefficients (n_vars, var_name, prob_coefs, startdate, enddate, Times, site_list, station_var, stnid, &
   & stnlat, stnlon, stnalt, forecast, output_file2, error)
    IF (error /= 0) RETURN
   END IF
 
   CALL save_coefficients (n_vars, var_name, coefs, startdate, enddate, Times, site_list, station_var, stnid, stnlat, &
  & stnlon, stnalt, forecast, output_file, error)
   IF (error /= 0) RETURN
 
  ELSE
   IF (mode == 2) THEN
 
    grid_list = config_values (13)
    CALL value (config_values(14), maxDistance, error)
    maxDistance = maxDistance * 0.539957
    IF (error /= 0) maxDistance = - 1
 
    IF (len(trim(grid_list)) == 0) THEN
     PRINT *, "ERROR: Failed to read in one more more required model config parameters. (GRID_LIST)"
     RETURN
    END IF
 
 
    CALL read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, stn_slp_e, nstations, error)
    IF (error /= 0) RETURN
 
    CALL read_nc_grid (grid_list, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
 
      !allocate 1-d grid variables
    ALLOCATE (grdlat(nx*ny))
    ALLOCATE (grdlon(nx*ny))
    ALLOCATE (grdalt(nx*ny))
    ALLOCATE (grd_slp_n(nx*ny))
    ALLOCATE (grd_slp_e(nx*ny))
    ALLOCATE (mask_1d(nx*ny))
 
 
    grdlat = reshape (lat, (/ nx*ny /))
    grdlon = reshape (lon, (/ nx*ny /))
    grdalt = reshape (elev, (/ nx*ny /))
    grd_slp_n = reshape (grad_n, (/ nx*ny /))
    grd_slp_e = reshape (grad_e, (/ nx*ny /))
    mask_1d = reshape (mask, (/ nx*ny /))
 
    ngrid = nx * ny
    ALLOCATE (X(nstations, 6))
    ALLOCATE (Z(ngrid, 6))
    X (:, 1) = 1.0
    X (:, 2) = stnlat (:)
    X (:, 3) = stnlon (:)
    X (:, 4) = stnalt (:)
    X (:, 5) = stn_slp_n (:)
    X (:, 6) = stn_slp_e (:)
 
    Z (:, 1) = 1.0
    Z (:, 2) = grdlat (:)
    Z (:, 3) = grdlon (:)
    Z (:, 4) = grdalt (:)
    Z (:, 5) = grd_slp_n (:)
    Z (:, 6) = grd_slp_e (:)
 
 
    CALL get_time_list (startdate, enddate, Times)
 
    ntimes = size (Times)
 
    ALLOCATE (mean_autocorr(ntimes))
    ALLOCATE (mean_tp_corr(ntimes))
    ALLOCATE (y_mean(ngrid, ntimes))
    ALLOCATE (y_std(ngrid, ntimes))
    ALLOCATE (y_std_all(ngrid, ntimes))
    ALLOCATE (y_max(ngrid, ntimes))
    ALLOCATE (y_min(ngrid, ntimes))
 
    CALL estimate_precip (X, Z, nstations, ngrid, maxDistance, Times, stnid, station_var, site_var, site_var_t, &
   & site_list, PCP, POP, pcperror, tmean, tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, &
   & y_std_all, y_min, y_max, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
    IF (error /= 0) RETURN
 
 
    PRINT *, 'Creating output file'
 
    CALL save_precip (PCP, POP, pcperror, tmean, tmean_err, trange, trange_err, nx, ny, grdlat, grdlon, grdalt, Times, &
   & mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, output_file, error, pcp_2, pop_2, pcperror_2, &
   & tmean_2, tmean_err_2, trange_2, trange_err_2)
    IF (error /= 0) RETURN
 
 
   END IF
  END IF
 
 
END PROGRAM precip
 
 
 
SUBROUTINE get_time_list (startdate, enddate, Times)
  USE utim
  USE type
  IMPLICIT NONE
 
  CHARACTER (LEN=100), INTENT (IN) :: startdate, enddate
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: Times (:)
 
  INTEGER (I4B) :: T, ntimes, sday, eday, error
  INTEGER (I4B) :: sec, Min, hour, day, month, year
  REAL (DP) :: utime
 
  CALL parse_date (startdate, year, month, day, hour, Min, sec, error)
  sday = julian_date (day, month, year)
  CALL parse_date (enddate, year, month, day, hour, Min, sec, error)
  eday = julian_date (day, month, year)
  ntimes = eday - sday + 1
  ALLOCATE (Times(ntimes))
 
 
  utime = date_to_unix (startdate)
  DO T = 1, ntimes, 1
   IF (utime > date_to_unix(enddate)) EXIT
   Times (T) = utime
   utime = utime + 86400
  END DO
 
END SUBROUTINE get_time_list
