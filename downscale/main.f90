 
Program precip
  Use type
  Use strings
  Implicit None
 
  Interface
    Subroutine get_time_list (startdate, enddate, T)
      Use utim
      Use type
      Character (Len=100), Intent (In) :: startdate, enddate
      Real (DP), Allocatable, Intent (Out) :: T (:)
    End Subroutine get_time_list
 
    Subroutine read_refcst (startdate, enddate, file_var, perturbation, var_name, forecast, V, X, Y, T, error)
      Use type
      Character (Len=100), Intent (In) :: startdate, enddate, file_var, perturbation
      Character (Len=*), Intent (In) :: var_name
      Integer (I4B), Intent (In) :: forecast
      Real (DP), Allocatable, Intent (Out) :: V (:, :), X (:), Y (:)
      Real (DP), Allocatable, Intent (Out) :: T (:)
      Integer, Intent (Out) :: error
    End Subroutine read_refcst
 
    Subroutine read_station_list (file_name, id, name, lat, lon, alt, sslp_n, sslp_e, n_stations, error)
      Use type
      Character (Len=500), Intent (In) :: file_name
      Character (Len=100), Allocatable, Intent (Out) :: id (:), name (:)
      Real (DP), Allocatable, Intent (Out) :: lat (:), lon (:), alt (:), sslp_n (:), sslp_e (:)
      Integer (I4B), Intent (Out) :: n_stations
      Integer, Intent (Out) :: error
    End Subroutine read_station_list
 
    Subroutine read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
      Use type
      Character (Len=500), Intent (In) :: file_name
      Real (DP), Allocatable, Intent (Out) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
      Integer (I4B), Intent (Out) :: nx, ny
      Integer, Intent (Out) :: error
    End Subroutine read_grid_list
 
    Subroutine read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
      Use netcdf
      Use type
 
      Character (Len=500), Intent (In) :: file_name
      Real (DP), Allocatable, Intent (Out) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), grad_e (:, :), mask &
     & (:, :)
      Integer (I4B), Intent (Out) :: nx, ny
      Integer, Intent (Out) :: error
    End Subroutine read_nc_grid
 
 
    Subroutine estimate_coefficients (D, nvars, lats, lons, Times, stnid, stnlat, stnlon, stnalt, stnvar, site_var, &
   & site_list, C, POC, error)
      Use type
      Real (DP), Intent (In) :: D (:, :, :), lats (:), lons (:)
      Real (DP), Intent (In) :: Times (:)
      Integer (I4B), Intent (In) :: nvars
      Character (Len=100), Intent (In) :: stnid (:)
      Real (DP), Intent (In) :: stnlat (:), stnlon (:), stnalt (:)
      Character (Len=100), Intent (In) :: stnvar
      Character (Len=100), Intent (In) :: site_var
      Character (Len=500), Intent (In) :: site_list
      Real (DP), Allocatable, Intent (Out) :: C (:, :, :), POC (:, :, :)
      Integer, Intent (Out) :: error
    End Subroutine estimate_coefficients
 
    Subroutine save_coefficients (n_vars, var_names, coefs, startdate, enddate, Times, site_list, station_var, stnid, &
   & stnlat, stnlon, stnalt, forecast, file, error)
      Use netcdf
      Use type
      Character (Len=100), Intent (In) :: var_names (:), stnid (:)
      Integer, Intent (In) :: n_vars, forecast
      Character (Len=100), Intent (In) :: startdate, enddate, station_var
      Character (Len=500), Intent (In) :: file, site_list
      Real (DP), Intent (In) :: stnlat (:), stnlon (:), stnalt (:)
      Real (DP), Intent (In) :: coefs (:, :, :)
      Real (DP), Intent (In) :: Times (:)
      Integer, Intent (Out) :: error
    End Subroutine save_coefficients
 
    Subroutine estimate_precip (X, Z, nsta, ngrid, maxDistance, Times, stnid, stnvar, site_var, site_var_t, site_list, &
   & PCP, POP, PCPERR, tmean, tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, &
   & y_min, y_max, error, pcp_2, pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
      Use type
      Real (DP), Intent (In) :: X (:, :), Z (:, :)
      Real (DP), Intent (In) :: maxDistance
      Integer (I4B), Intent (In) :: nsta, ngrid
      Real (DP), Intent (In) :: Times (:)
      Character (Len=100), Intent (In) :: stnid (:)
      Character (Len=100), Intent (In) :: stnvar, site_var, site_var_t
      Character (Len=500), Intent (In) :: site_list
      Real (SP), Allocatable, Intent (Out) :: PCP (:, :), POP (:, :), PCPERR (:, :)
      Real (SP), Allocatable, Intent (Out) :: tmean (:, :), tmean_err (:, :)!OLS tmean estimate and error
      Real (SP), Allocatable, Intent (Out) :: trange (:, :), trange_err (:, :)!OLS trange estimate and error
 
      Real (SP), Allocatable, Intent (Out) :: pcp_2 (:, :), pop_2 (:, :), pcperr_2 (:, :)
      Real (SP), Allocatable, Intent (Out) :: tmean_2 (:, :), tmean_err_2 (:, :)!OLS tmean estimate and error
      Real (SP), Allocatable, Intent (Out) :: trange_2 (:, :), trange_err_2 (:, :)!OLS trange estimate and error
 
 
      Integer, Intent (Out) :: error
      Real (DP), Intent (Out) :: mean_autocorr (:)!mean autocorrelation from all stations over entire time period
      Real (DP), Intent (Out) :: mean_tp_corr (:)!mean correlation for mean temp and precip
 
      Real (DP), Intent (Out) :: y_mean (:, :), y_std (:, :), y_std_all (:, :)!std and mean of transformed time step precipitation
      Real (DP), Intent (Out) :: y_min (:, :), y_max (:, :)!min,max of normalized time step precip
    End Subroutine estimate_precip
 
    Subroutine save_precip (PCP, POP, pcperror, tmean, tmean_err, trange, trange_err, nx, ny, grdlat, grdlon, grdalt, &
   & Times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, file, error, pcp_2, pop_2, pcperror_2, &
   & tmean_2, tmean_err_2, trange_2, trange_err_2)
      Use netcdf
      Use type
 
      Real (SP), Intent (In) :: PCP (:, :), POP (:, :), pcperror (:, :)
      Real (SP), Intent (In) :: tmean (:, :), tmean_err (:, :), trange (:, :), trange_err (:, :)
 
      Real (SP), Intent (In) :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
      Real (SP), Intent (In) :: tmean_2 (:, :), tmean_err_2 (:, :), trange_2 (:, :), trange_err_2 (:, :)
 
      Integer (I4B), Intent (In) :: nx, ny
      Real (DP), Intent (In) :: grdlat (:), grdlon (:), grdalt (:)
      Real (DP), Intent (In) :: Times (:)
      Real (DP), Intent (In) :: mean_autocorr (:), mean_tp_corr (:)
 
      Real (DP), Intent (In) :: y_mean (:, :), y_std (:, :), y_std_all (:, :)
      Real (DP), Intent (In) :: y_min (:, :), y_max (:, :)
      Character (Len=500), Intent (In) :: file
      Integer, Intent (Out) :: error
    End Subroutine save_precip
 
  End Interface
 
  Character (Len=100) :: config_file
  Integer, Parameter :: nconfigs = 15
  Character (Len=500) :: config_names (nconfigs)
  Character (Len=500) :: config_values (nconfigs)
  Character (Len=500) :: site_list, output_file, output_file2, grid_list
  Character (Len=100) :: startdate, enddate, perturbation, station_var, site_var, site_var_t
  Character (Len=100), Allocatable :: file_var (:), var_name (:)
  Character (Len=100), Allocatable :: stnid (:), stnname (:)
 
  Character (Len=2000) :: arg !command line arg for configuration file
 
  Integer :: i, error, n_vars, nfile_var, nvar_name, forecast, mode
  Integer :: nstations, lenfile
 
  Real (DP), Allocatable :: Y (:, :, :), Vals (:, :), lats (:), lons (:)
  Real (DP), Allocatable :: coefs (:, :, :), prob_coefs (:, :, :)
  Real (DP), Allocatable :: stnlat (:), stnlon (:), stnalt (:), stn_slp_n (:), stn_slp_e (:)
  Real (DP), Allocatable :: Times (:)
 
  Real (DP), Allocatable :: mean_autocorr (:)!mean auto correlation for all stations over entire time period
  Real (DP), Allocatable :: mean_tp_corr (:)!mean correlation between precip and trange (31-day moving avg anomaly)
  Real (DP), Allocatable :: y_mean (:, :)
  Real (DP), Allocatable :: y_std (:, :)
  Real (DP), Allocatable :: y_std_all (:, :)
  Real (DP), Allocatable :: y_max (:, :)
  Real (DP), Allocatable :: y_min (:, :)
 
  Real (DP), Allocatable :: lat (:, :)
  Real (DP), Allocatable :: lon (:, :)
  Real (DP), Allocatable :: elev (:, :)
  Real (DP), Allocatable :: grad_n (:, :)
  Real (DP), Allocatable :: grad_e (:, :)
  Real (DP), Allocatable :: mask (:, :)
 
 
  Real (DP), Allocatable :: X (:, :), Z (:, :)
  Integer :: ngrid
 
  Integer (I4B) :: nx, ny, ntimes
  Real (DP) :: maxDistance
  Real (DP), Allocatable :: grdlat (:), grdlon (:), grdalt (:), grd_slp_n (:), grd_slp_e (:), mask_1d (:)
  Real (SP), Allocatable :: PCP (:, :), POP (:, :), pcperror (:, :)
  Real (SP), Allocatable :: tmean (:, :), tmean_err (:, :)
  Real (SP), Allocatable :: trange (:, :), trange_err (:, :)
 
  Real (SP), Allocatable :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
  Real (SP), Allocatable :: tmean_2 (:, :), tmean_err_2 (:, :)
  Real (SP), Allocatable :: trange_2 (:, :), trange_err_2 (:, :)
 
!code starts below
 
!mode 2
  config_file = "config_pnw.txt"
!mode 1
!  config_file = "config_prcp.txt"
 
 
!get config_file filename from command line
  i = 0
  Do
    Call get_command_argument (i, arg)
    If (i .Eq. 1) config_file = arg
    If (LEN_TRIM(arg) == 0) Exit
    i = i + 1
  End Do
 
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
 
  Call read_config (config_file, nconfigs, config_names, config_values)
  Call value (config_values(1), mode, error)
  If (error /= 0) Then
    Print *, "ERROR: Failed to read mode from config file."
    Return
  End If
 
  startdate = config_values (2)
  enddate = config_values (3)
  site_list = config_values (4)
  site_var = config_values (5)
  site_var_t = config_values (15)
  station_var = config_values (6)
  output_file = config_values (12)
 
 
  If (len(trim(startdate)) == 0 .Or. len(trim(enddate)) == 0 .Or. len(trim(site_list)) == 0 .Or. len(trim(output_file)) &
 & == 0) Then
    Print *, "ERROR: Failed to read in one more more required config parameters. (START_DATE, END_DATE, SITE_LIST or OU&
   &TPUT_FILE)"
    Return
  End If
 
  If (mode == 1) Then
     ! Model Variables
    perturbation = config_values (7)
    Call value (config_values(8), forecast, error)
    If (error /= 0) forecast = - 1
    Call value (config_values(9), n_vars, error)
    If (error /= 0) n_vars = 0
 
    If (len(trim(perturbation)) == 0 .Or. n_vars == 0 .Or. forecast ==-1) Then
      Print *, "ERROR: Failed to read in one or more required model config parameters. (PERTURBATION, NUMBER_VARS or FO&
     &RECAST)"
      Return
    End If
 
    If (n_vars > 1) Then
      Allocate (file_var(n_vars))
      Allocate (var_name(n_vars))
      Call parse (config_values(10), ",", file_var, nfile_var)
      Call parse (config_values(11), ",", var_name, nvar_name)
      If (nfile_var /= n_vars .Or. nvar_name /= n_vars) Then
        Print *, "ERROR: Number of variables in config file does not match."
        Return
      End If
    Else
      Allocate (file_var(1))
      Allocate (var_name(1))
      file_var (1) = config_values (10)
      var_name (1) = config_values (11)
      If (len(trim(file_var(1))) == 0 .Or. len(trim(var_name(1))) == 0) Then
        Print *, "ERROR: Failed to read in one or more required model config parameters. (FILE_VARIABLE or VARIABLE_NAM&
       &E)"
        Return
      End If
 
    End If
 
    Do i = 1, n_vars, 1
      Call read_refcst (startdate, enddate, file_var(i), perturbation, var_name(i), forecast, Vals, lats, lons, Times, &
     & error)
      If (error /= 0) Return
      If (i == 1) Then
        Allocate (Y(n_vars, size(lons)*size(lats), size(Times)))
      End If
      Y (i, :, :) = Vals (:, :)
      Deallocate (Vals)
    End Do
 
    Call read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, stn_slp_e, nstations, error)
    If (error /= 0) Return
 
    Call estimate_coefficients (Y, n_vars, lats, lons, Times, stnid, stnlat, stnlon, stnalt, station_var, site_var, &
   & site_list, coefs, prob_coefs, error)
    If (error /= 0) Return
 
 
    If (trim(station_var) .Eq. "PRCP") Then
      lenfile = LEN_TRIM (output_file)
      output_file2 (:) = " "
      output_file2 (1:5) = "prob_"
      output_file2 (6:lenfile+6) = output_file (1:lenfile)
 
      Call save_coefficients (n_vars, var_name, prob_coefs, startdate, enddate, Times, site_list, station_var, stnid, &
     & stnlat, stnlon, stnalt, forecast, output_file2, error)
      If (error /= 0) Return
    End If
 
    Call save_coefficients (n_vars, var_name, coefs, startdate, enddate, Times, site_list, station_var, stnid, stnlat, &
   & stnlon, stnalt, forecast, output_file, error)
    If (error /= 0) Return
 
  Else
    If (mode == 2) Then
 
      grid_list = config_values (13)
      Call value (config_values(14), maxDistance, error)
      maxDistance = maxDistance * 0.539957
      If (error /= 0) maxDistance = - 1
 
      If (len(trim(grid_list)) == 0) Then
        Print *, "ERROR: Failed to read in one more more required model config parameters. (GRID_LIST)"
        Return
      End If
 
 
      Call read_station_list (site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, stn_slp_e, nstations, &
     & error)
      If (error /= 0) Return
 
      Call read_nc_grid (grid_list, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
 
      !allocate 1-d grid variables
      Allocate (grdlat(nx*ny))
      Allocate (grdlon(nx*ny))
      Allocate (grdalt(nx*ny))
      Allocate (grd_slp_n(nx*ny))
      Allocate (grd_slp_e(nx*ny))
      Allocate (mask_1d(nx*ny))
 
 
      grdlat = reshape (lat, (/ nx*ny /))
      grdlon = reshape (lon, (/ nx*ny /))
      grdalt = reshape (elev, (/ nx*ny /))
      grd_slp_n = reshape (grad_n, (/ nx*ny /))
      grd_slp_e = reshape (grad_e, (/ nx*ny /))
      mask_1d = reshape (mask, (/ nx*ny /))
 
      ngrid = nx * ny
      Allocate (X(nstations, 6))
      Allocate (Z(ngrid, 6))
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
 
 
      Call get_time_list (startdate, enddate, Times)
 
      ntimes = size (Times)
 
      Allocate (mean_autocorr(ntimes))
      Allocate (mean_tp_corr(ntimes))
      Allocate (y_mean(ngrid, ntimes))
      Allocate (y_std(ngrid, ntimes))
      Allocate (y_std_all(ngrid, ntimes))
      Allocate (y_max(ngrid, ntimes))
      Allocate (y_min(ngrid, ntimes))
 
      Call estimate_precip (X, Z, nstations, ngrid, maxDistance, Times, stnid, station_var, site_var, site_var_t, &
     & site_list, PCP, POP, pcperror, tmean, tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, &
     & y_std_all, y_min, y_max, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
      If (error /= 0) Return
 
 
      Print *, 'Creating output file'
 
      Call save_precip (PCP, POP, pcperror, tmean, tmean_err, trange, trange_err, nx, ny, grdlat, grdlon, grdalt, &
     & Times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, output_file, error, pcp_2, pop_2, &
     & pcperror_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
      If (error /= 0) Return
 
 
    End If
  End If
 
 
End Program precip
 
 
 
Subroutine get_time_list (startdate, enddate, Times)
  Use utim
  Use type
  Implicit None
 
  Character (Len=100), Intent (In) :: startdate, enddate
  Real (DP), Allocatable, Intent (Out) :: Times (:)
 
  Integer (I4B) :: T, ntimes, sday, eday, error
  Integer (I4B) :: sec, Min, hour, day, month, year
  Real (DP) :: utime
 
  Call parse_date (startdate, year, month, day, hour, Min, sec, error)
  sday = julian_date (day, month, year)
  Call parse_date (enddate, year, month, day, hour, Min, sec, error)
  eday = julian_date (day, month, year)
  ntimes = eday - sday + 1
  Allocate (Times(ntimes))
 
 
  utime = date_to_unix (startdate)
  Do T = 1, ntimes, 1
    If (utime > date_to_unix(enddate)) Exit
    Times (T) = utime
    utime = utime + 86400
  End Do
 
End Subroutine get_time_list
