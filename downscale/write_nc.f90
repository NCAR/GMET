 
Subroutine save_coefficients (n_vars, var_names, coefs, startdate, enddate, times, site_list, station_var, stnid, &
& stnlat, stnlon, stnalt, forecast, file, error)
  Use netcdf
  Use type
  Implicit None
 
  Character (Len=100), Intent (In) :: var_names (:), stnid (:)
  Integer, Intent (In) :: n_vars, forecast
  Character (Len=100), Intent (In) :: startdate, enddate, station_var
  Character (Len=500), Intent (In) :: file, site_list
  Real (DP), Intent (In) :: stnlat (:), stnlon (:), stnalt (:)
  Real (DP), Intent (In) :: coefs (:, :, :)
  Real (DP), Intent (In) :: times (:)
 
  Integer, Intent (Out) :: error
 
  ! Dimension names
  Character (Len=*), Parameter :: STN_NAME = "station"
  Character (Len=*), Parameter :: VAR_NAME = "variable"
  Character (Len=*), Parameter :: CHAR_NAME = "string"
  Character (Len=*), Parameter :: TIME_NAME = "time"
  Character (Len=*), Parameter :: REC_NAME = "run"
  ! Variable Names
  Character (Len=*), Parameter :: VARS_NAME = "variable_name"
  Character (Len=*), Parameter :: STN_FILE_NAME = "station_file_name"
  Character (Len=*), Parameter :: STN_ID_NAME = "station_id"
  Character (Len=*), Parameter :: STN_LAT_NAME = "station_latitude"
  Character (Len=*), Parameter :: STN_LON_NAME = "station_longitude"
  Character (Len=*), Parameter :: STN_ALT_NAME = "station_altitude"
  Character (Len=*), Parameter :: FORECAST_NAME = "forecast_hr"
  Character (Len=*), Parameter :: STARTDATE_NAME = "run_start_date"
  Character (Len=*), Parameter :: ENDDATE_NAME = "run_end_date"
  Character (Len=*), Parameter :: TIME_VAR_NAME = "time"
  Character (Len=*), Parameter :: COEFS_NAME = "coefficient"
  Character (Len=*), Parameter :: CONSTANT_NAME = "constant_term"
  ! Units
  Character (Len=*), Parameter :: UNITS = "units"
  Character (Len=*), Parameter :: COEFS_UNITS = "linear_equation_coefficient"
  Character (Len=*), Parameter :: LAT_UNITS = "degrees_north"
  Character (Len=*), Parameter :: LON_UNITS = "degrees_east"
  Character (Len=*), Parameter :: ALT_UNITS = "feet"
  Character (Len=*), Parameter :: FORECAST_UNITS = "hours"
  Character (Len=*), Parameter :: DATE_UNITS = "YYYYMMDD"
  Character (Len=*), Parameter :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  Character (Len=*), Parameter :: FILL = "_FillValue"
 
  Integer :: n_stns, n_chars, n_times
  Integer :: ncid, stn_dimid, vars_dimid, char_dimid, time_dimid, rec_dimid
  Integer :: coefs_varid, time_varid, vars_varid, forecast_varid, startdate_varid, enddate_varid
  Integer :: stn_id_varid, stn_lat_varid, stn_lon_varid, stn_alt_varid
  Integer :: charids (2)
  Integer :: dimids (4)
  Integer :: count4 (4), start4 (4), count3 (3), start3 (3), count2 (2), start2 (2), count1 (1), start1 (1)
  Integer :: nstns, nvars, nrecs, ntimes, rec, i
 
  Character (Len=100) :: startdates, enddates
  Integer, Allocatable :: forecasts (:)
 
  Integer :: forecast_arr (1)
 
  n_chars = 100
  n_stns = size (stnlat)
  n_times = size (times)
 
  error = nf90_open (file, nf90_write, ncid)
  If (error /= nf90_noerr) Then
    error = 0
     ! Create the file.
    Call check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
    If (error /= 0) Return
 
     ! Define the dimensions.
    Call check (nf90_def_dim(ncid, STN_NAME, n_stns, stn_dimid), "station dim def error", error)
    Call check (nf90_def_dim(ncid, VAR_NAME, n_vars+1, vars_dimid), "variable dim def error", error)
    Call check (nf90_def_dim(ncid, CHAR_NAME, n_chars, char_dimid), "char dim def error", error)
    Call check (nf90_def_dim(ncid, TIME_NAME, n_times, time_dimid), "time dim def error", error)
    Call check (nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid), "rec dim def error", error)
    If (error /= 0) Return
 
     ! Define the variables.
    charids = (/ char_dimid, vars_dimid /)
    Call check (nf90_def_var(ncid, VARS_NAME, NF90_CHAR, charids, vars_varid), "variable_name var def error", error)
    charids = (/ char_dimid, stn_dimid /)
    Call check (nf90_def_var(ncid, STN_ID_NAME, NF90_CHAR, charids, stn_id_varid), "station_id var def error", error)
    Call check (nf90_def_var(ncid, STN_LAT_NAME, NF90_DOUBLE, stn_dimid, stn_lat_varid), "station_latitude var def erro&
   &r", error)
    Call check (nf90_def_var(ncid, STN_LON_NAME, NF90_DOUBLE, stn_dimid, stn_lon_varid), "station_longitude var def err&
   &or", error)
    Call check (nf90_def_var(ncid, STN_ALT_NAME, NF90_DOUBLE, stn_dimid, stn_alt_varid), "station_altitude var def erro&
   &r", error)
    If (error /= 0) Return
 
    charids = (/ char_dimid, rec_dimid /)
    Call check (nf90_def_var(ncid, STARTDATE_NAME, NF90_CHAR, charids, startdate_varid), "start_date var def error", &
   & error)
    Call check (nf90_def_var(ncid, ENDDATE_NAME, NF90_CHAR, charids, enddate_varid), "end_date var def error", error)
    Call check (nf90_def_var(ncid, FORECAST_NAME, NF90_INT, rec_dimid, forecast_varid), "forecast var def error", &
   & error)
    If (error /= 0) Return
 
    Call check (nf90_def_var(ncid, TIME_VAR_NAME, NF90_DOUBLE, time_dimid, time_varid), "time var def error", error)
    dimids = (/ stn_dimid, time_dimid, vars_dimid, rec_dimid /)
    Call check (nf90_def_var(ncid, COEFS_NAME, NF90_DOUBLE, dimids, coefs_varid), "coefficient var def error", error)
    If (error /= 0) Return
 
     ! Add attributes.
    Call check (nf90_put_att(ncid, stn_lat_varid, UNITS, LAT_UNITS), "station_lat units attribute error", error)
    Call check (nf90_put_att(ncid, stn_lon_varid, UNITS, LON_UNITS), "station_lon units attribute error", error)
    Call check (nf90_put_att(ncid, stn_alt_varid, UNITS, ALT_UNITS), "station_alt units attribute error", error)
    Call check (nf90_put_att(ncid, startdate_varid, UNITS, DATE_UNITS), "start_date units attribute error", error)
    Call check (nf90_put_att(ncid, enddate_varid, UNITS, DATE_UNITS), "end_date units attribute error", error)
    Call check (nf90_put_att(ncid, forecast_varid, UNITS, FORECAST_UNITS), "forecase units attribute error", error)
    Call check (nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)
    Call check (nf90_put_att(ncid, coefs_varid, UNITS, COEFS_UNITS), "coefficient units attribute error", error)
    Call check (nf90_put_att(ncid, NF90_GLOBAL, STN_FILE_NAME, site_list), "station_file_name global attribute error", &
   & error)
 
     ! End define mode.
    Call check (nf90_enddef(ncid), "end define mode error", error)
    If (error /= 0) Return
 
    Call check (nf90_put_var(ncid, stn_lat_varid, stnlat), "put staion_lat error", error)
    Call check (nf90_put_var(ncid, stn_lon_varid, stnlon), "put staion_lon error", error)
    Call check (nf90_put_var(ncid, stn_alt_varid, stnalt), "put staion_alt error", error)
    Call check (nf90_put_var(ncid, time_varid, times), "put times error", error)
 
    count2 = (/ 1, 1 /)
    start2 = (/ 1, 1 /)
    Do i = 1, n_stns, 1
      count2 (1) = len (trim(stnid(i)))
      start2 (2) = i
      Call check (nf90_put_var(ncid, stn_id_varid, stnid(i), start=start2, count=count2), "put staion_id error", error)
      If (error /= 0) Return
    End Do
 
    count2 (1) = len (trim(CONSTANT_NAME))
    start2 (2) = 1
    Call check (nf90_put_var(ncid, vars_varid, CONSTANT_NAME, start=start2, count=count2), "put variable_name error", &
   & error)
    If (error /= 0) Return
    Do i = 1, n_vars, 1
      count2 (1) = len (trim(var_names(i)))
      start2 (2) = i + 1
      Call check (nf90_put_var(ncid, vars_varid, var_names(i), start=start2, count=count2), "put variable_name error", &
     & error)
      If (error /= 0) Return
    End Do
 
    rec = 1
    nrecs = 0
 
  Else
 
     ! File already exists, get dim and var ids
    Call check (nf90_inq_dimid(ncid, STN_NAME, stn_dimid), "station dim inq error", error)
    Call check (nf90_inq_dimid(ncid, VAR_NAME, vars_dimid), "variable dim inq error", error)
    Call check (nf90_inq_dimid(ncid, CHAR_NAME, char_dimid), "char dim inq error", error)
    Call check (nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
    Call check (nf90_inq_dimid(ncid, REC_NAME, rec_dimid), "run dim inq error", error)
    If (error /= 0) Return
 
    Call check (nf90_inq_varid(ncid, VARS_NAME, vars_varid), "variable_name var inq error", error)
    Call check (nf90_inq_varid(ncid, STN_ID_NAME, stn_id_varid), "station_id var inq error", error)
    Call check (nf90_inq_varid(ncid, STN_LAT_NAME, stn_lat_varid), "station_latitude var inq error", error)
    Call check (nf90_inq_varid(ncid, STN_LON_NAME, stn_lon_varid), "station_longitude var inq error", error)
    Call check (nf90_inq_varid(ncid, STN_ALT_NAME, stn_alt_varid), "station_altitude var inq error", error)
    Call check (nf90_inq_varid(ncid, STARTDATE_NAME, startdate_varid), "start_date var inq error", error)
    Call check (nf90_inq_varid(ncid, ENDDATE_NAME, enddate_varid), "end_date var inq error", error)
    Call check (nf90_inq_varid(ncid, FORECAST_NAME, forecast_varid), "forecast var inq error", error)
    Call check (nf90_inq_varid(ncid, TIME_VAR_NAME, time_varid), "time var inq error", error)
    Call check (nf90_inq_varid(ncid, COEFS_NAME, coefs_varid), "coefficient var inq error", error)
    If (error /= 0) Return
 
     ! Verify Dimensions match
    Call check (nf90_Inquire_Dimension(ncid, stn_dimid, len=nstns), "station dim len error", error)
    Call check (nf90_Inquire_Dimension(ncid, vars_dimid, len=nvars), "variable dim len error", error)
    Call check (nf90_Inquire_Dimension(ncid, time_dimid, len=ntimes), "time dim len error", error)
    Call check (nf90_Inquire_Dimension(ncid, rec_dimid, len=nrecs), "run dim len error", error)
    If (error /= 0) Return
 
    If (n_stns /= nstns .Or. nvars /= n_vars+1 .Or. ntimes /= n_times) Then
      Print *, "Error dimensions in output file do not match current run."
      error = 1
      Return
    End If
 
    Allocate (forecasts(nrecs))
 
    Call check (nf90_get_var(ncid, forecast_varid, forecasts), "error getting file forecasts", error)
    If (error /= 0) Return
 
    rec = nrecs + 1
    Do i = 1, nrecs, 1
      count2 = (/ 100, 1 /)
      start2 = (/ 1, i /)
      Call check (nf90_get_var(ncid, startdate_varid, startdates, start=start2, count=count2), "error getting file star&
     &tdates", error)
      Call check (nf90_get_var(ncid, enddate_varid, enddates, start=start2, count=count2), "error getting file enddates&
     &", error)
      If (forecasts(i) == forecast .And. startdates(1:8) == startdate(1:8) .And. enddates(1:8) == enddate(1:8)) Then
        Print *, "WARNING, overwriting data in output file, record ", i
        rec = i
      End If
    End Do
 
 
  End If
 
  count2 = (/ len (trim(startdate)), 1 /)
  start2 = (/ 1, rec /)
  Call check (nf90_put_var(ncid, startdate_varid, startdate, start=start2, count=count2), "put start_date error", &
 & error)
  If (error /= 0) Return
 
  count2 = (/ len (trim(enddate)), 1 /)
  start2 = (/ 1, rec /)
  Call check (nf90_put_var(ncid, enddate_varid, enddate, start=start2, count=count2), "put start_date error", error)
  If (error /= 0) Return
 
  count1 (1) = 1
  start1 (1) = rec
  forecast_arr = forecast
  Call check (nf90_put_var(ncid, forecast_varid, forecast_arr, start=start1, count=count1), "put forecast error", &
 & error)
  If (error /= 0) Return
 
  count4 = (/ n_stns, n_times, n_vars + 1, 1 /)
  start4 = (/ 1, 1, 1, rec /)
  forecast_arr = forecast
  Call check (nf90_put_var(ncid, coefs_varid, coefs, start=start4, count=count4), "put coefficients error", error)
  If (error /= 0) Return
 
  Call check (nf90_close(ncid), "closing file error", error)
 
Contains
  Subroutine check (status, info, error)
    Integer, Intent (In) :: status
    Character (Len=*), Intent (In) :: info
    Integer, Intent (Out) :: error
 
    If (status /= nf90_noerr) Then
      Print *, trim (info) // ": " // trim (nf90_strerror(status))
      error = 1
    End If
  End Subroutine check
End Subroutine save_coefficients
 
Subroutine save_precip (pcp, pop, pcperror, tmean, tmean_error, trange, trange_error, nx, ny, grdlat, grdlon, grdalt, &
& times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, file, error, pcp_2, pop_2, pcperror_2, &
& tmean_2, tmean_error_2, trange_2, trange_error_2)
  Use netcdf
  Use type
  Implicit None
 
  Real (SP), Intent (In) :: pcp (:, :), pop (:, :), pcperror (:, :)
  Real (SP), Intent (In) :: tmean (:, :), tmean_error (:, :), trange (:, :), trange_error (:, :)
 
  Real (SP), Intent (In) :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
  Real (SP), Intent (In) :: tmean_2 (:, :), tmean_error_2 (:, :), trange_2 (:, :), trange_error_2 (:, :)
 
  Integer (I4B), Intent (In) :: nx, ny
  Real (DP), Intent (In) :: grdlat (:), grdlon (:), grdalt (:)
  Real (DP), Intent (In) :: times (:)
  Real (DP), Intent (In) :: mean_autocorr (:), mean_tp_corr (:)
 
  Real (DP), Intent (In) :: y_mean (:, :), y_std (:, :), y_min (:, :), y_max (:, :), y_std_all (:, :)
  Character (Len=500), Intent (In) :: file
  Integer, Intent (Out) :: error
 
 
  ! Dimension names
  Character (Len=*), Parameter :: Y_NAME = "y"
  Character (Len=*), Parameter :: X_NAME = "x"
  Character (Len=*), Parameter :: TIME_NAME = "time"
 
  ! Variable Names
  Character (Len=*), Parameter :: LAT_NAME = "latitude"
  Character (Len=*), Parameter :: LON_NAME = "longitude"
  Character (Len=*), Parameter :: autoc_name = "auto_corr"
  Character (Len=*), Parameter :: tpc_name = "tp_corr"
  Character (Len=*), Parameter :: ALT_NAME = "altitude"
  Character (Len=*), Parameter :: PCP_NAME = "pcp"
  Character (Len=*), Parameter :: POP_NAME = "pop"
  Character (Len=*), Parameter :: PCP_ERROR_NAME = "pcp_error"
  Character (Len=*), Parameter :: tmean_name = "tmean"
  Character (Len=*), Parameter :: tmean_error_name = "tmean_error"
  Character (Len=*), Parameter :: trange_name = "trange"
  Character (Len=*), Parameter :: trange_error_name = "trange_error"
  Character (Len=*), Parameter :: y_mean_name = "ymean"
  Character (Len=*), Parameter :: y_std_name = "ystd"
  Character (Len=*), Parameter :: y_stdall_name = "ystd_all"
  Character (Len=*), Parameter :: y_min_name = "ymax"
  Character (Len=*), Parameter :: y_max_name = "ymin"
  Character (Len=*), Parameter :: PCP_NAME_2 = "pcp_2"
  Character (Len=*), Parameter :: POP_NAME_2 = "pop_2"
  Character (Len=*), Parameter :: PCP_ERROR_NAME_2 = "pcp_error_2"
  Character (Len=*), Parameter :: tmean_name_2 = "tmean_2"
  Character (Len=*), Parameter :: tmean_error_name_2 = "tmean_error_2"
  Character (Len=*), Parameter :: trange_name_2 = "trange_2"
  Character (Len=*), Parameter :: trange_error_name_2 = "trange_error_2"
 
 
  Character (Len=*), Parameter :: LONG_NAME = "long_name"
  Character (Len=*), Parameter :: PCP_LONG_NAME = "estimated precip in normal space"
  Character (Len=*), Parameter :: POP_LONG_NAME = "probability of precipitation occurrence"
  Character (Len=*), Parameter :: PCP_ERROR_LONG_NAME = "error in estimated precip"
  Character (Len=*), Parameter :: tmean_long_name = "estimated daily mean temperature"
  Character (Len=*), Parameter :: tmean_error_long_name = "error in estimated daily mean temp"
  Character (Len=*), Parameter :: trange_long_name = "estimated diurnal range"
  Character (Len=*), Parameter :: trange_error_long_name = "error in estimated diurnal range"
  Character (Len=*), Parameter :: autoc_long_name = "Lag-1 autocorrelation of temperature"
  Character (Len=*), Parameter :: tpc_long_name = "Correlation of diurnal range and precipitation"
  Character (Len=*), Parameter :: y_mean_long_name = "mean of transformed non-zero precip"
  Character (Len=*), Parameter :: y_std_long_name = "std. dev. of transformed non-zero precip"
  Character (Len=*), Parameter :: y_max_long_name = "max of normalized transformed non-zero precip"
  Character (Len=*), Parameter :: y_min_long_name = "min of normalized transformed non-zero precip"
  Character (Len=*), Parameter :: PCP_LONG_NAME_2 = "estimated precip in normal space (no slope)"
  Character (Len=*), Parameter :: POP_LONG_NAME_2 = "probability of precipitation occurrence (no slope)"
  Character (Len=*), Parameter :: PCP_ERROR_LONG_NAME_2 = "error in estimated precip (no slope)"
  Character (Len=*), Parameter :: tmean_long_name_2 = "estimated daily mean temperature (no slope)"
  Character (Len=*), Parameter :: tmean_error_long_name_2 = "error in estimated daily mean temp (no slope)"
  Character (Len=*), Parameter :: trange_long_name_2 = "estimated diurnal range (no slope)"
  Character (Len=*), Parameter :: trange_error_long_name_2 = "error in estimated diurnal range (no slope)"
 
 
  ! Units
  Character (Len=*), Parameter :: UNITS = "units"
  Character (Len=*), Parameter :: PCP_UNITS = ""
  Character (Len=*), Parameter :: POP_UNITS = ""
  Character (Len=*), Parameter :: autoc_units = ""
  Character (Len=*), Parameter :: tpc_units = ""
  Character (Len=*), Parameter :: PCP_ERROR_UNITS = ""
  Character (Len=*), Parameter :: tmean_units = "deg_C"
  Character (Len=*), Parameter :: trange_units = "deg_C"
  Character (Len=*), Parameter :: tmean_error_units = "deg_C"
  Character (Len=*), Parameter :: trange_error_units = "deg_C"
  Character (Len=*), Parameter :: y_mean_units = ""
  Character (Len=*), Parameter :: y_std_units = ""
  Character (Len=*), Parameter :: y_max_units = ""
  Character (Len=*), Parameter :: y_min_units = ""
 
  Character (Len=*), Parameter :: LAT_UNITS = "degrees_north"
  Character (Len=*), Parameter :: LON_UNITS = "degrees_east"
  Character (Len=*), Parameter :: ALT_UNITS = "meters"
  Character (Len=*), Parameter :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  Character (Len=*), Parameter :: FILL = "_FillValue"
 
  Real (DP), Allocatable :: file_times (:)
 
  Integer :: n_chars, n_times, inx, iny
  Integer :: ncid, x_dimid, y_dimid, time_dimid
  Integer :: lat_varid, lon_varid, autoc_varid, alt_varid, time_varid, pcp_varid, pop_varid, pcp_error_varid, tpc_varid
  Integer :: tmean_varid, tmean_error_varid, trange_varid, trange_error_varid
  Integer :: ymean_varid, ystd_varid, ymax_varid, ymin_varid, ystdall_varid
  Integer :: count1 (1), start1 (1), count2 (2), start2 (2), count3 (3), start3 (3), dimids2 (2), dimids3 (3)
  Integer :: trec, nrecs, file_nx, file_ny, file_ntimes, i
 
  Integer :: pcp_varid_2, pop_varid_2, pcp_error_varid_2
  Integer :: tmean_varid_2, tmean_error_varid_2, trange_varid_2, trange_error_varid_2
 
  trec = 0
  n_chars = 100
  n_times = size (times)
  inx = nx
  iny = ny
 
  If (size(grdlat) /= inx*iny) Then
    Print *, "Error "
  End If
 
  error = nf90_open (file, nf90_write, ncid)
  If (error /= nf90_noerr) Then
    error = 0
     ! Create the file.
    Call check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
    If (error /= 0) Return
 
     ! Define the dimensions.
    Call check (nf90_def_dim(ncid, Y_NAME, iny, y_dimid), "y dim def error", error)
    Call check (nf90_def_dim(ncid, X_NAME, inx, x_dimid), "x dim def error", error)
    Call check (nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), "time dim def error", error)
    If (error /= 0) Return
 
     ! Define the variables.
    dimids2 = (/ x_dimid, y_dimid /)
    Call check (nf90_def_var(ncid, LAT_NAME, NF90_DOUBLE, dimids2, lat_varid), "lat var def error", error)
    Call check (nf90_def_var(ncid, LON_NAME, NF90_DOUBLE, dimids2, lon_varid), "lon var def error", error)
    Call check (nf90_def_var(ncid, ALT_NAME, NF90_DOUBLE, dimids2, alt_varid), "lon var def error", error)
    Call check (nf90_def_var(ncid, TIME_NAME, NF90_DOUBLE, time_dimid, time_varid), "time var def error", error)
 
     !correlation variables
    Call check (nf90_def_var(ncid, autoc_name, NF90_DOUBLE, time_dimid, autoc_varid), "auto correlaion var def error", &
   & error)
    Call check (nf90_def_var(ncid, tpc_name, NF90_DOUBLE, time_dimid, tpc_varid), "tp correlaion var def error", error)
 
 
 
    If (error /= 0) Return
 
    dimids3 = (/ x_dimid, y_dimid, time_dimid /)
    Call check (nf90_def_var(ncid, PCP_NAME, NF90_DOUBLE, dimids3, pcp_varid), "pcp var def error", error)
    Call check (nf90_def_var(ncid, POP_NAME, NF90_DOUBLE, dimids3, pop_varid), "pop var def error", error)
    Call check (nf90_def_var(ncid, PCP_ERROR_NAME, NF90_DOUBLE, dimids3, pcp_error_varid), "pcp_error var def error", &
   & error)
    If (error /= 0) Return
 
    Call check (nf90_def_var(ncid, tmean_name, NF90_DOUBLE, dimids3, tmean_varid), "tmean var def error", error)
    Call check (nf90_def_var(ncid, tmean_error_name, NF90_DOUBLE, dimids3, tmean_error_varid), "tmean error var def err&
   &or", error)
    Call check (nf90_def_var(ncid, trange_name, NF90_DOUBLE, dimids3, trange_varid), "trange var def error", error)
    Call check (nf90_def_var(ncid, trange_error_name, NF90_DOUBLE, dimids3, trange_error_varid), "trange error var def &
   &error", error)
    If (error /= 0) Return
 
    Call check (nf90_def_var(ncid, PCP_NAME_2, NF90_DOUBLE, dimids3, pcp_varid_2), "pcp var def error", error)
    Call check (nf90_def_var(ncid, POP_NAME_2, NF90_DOUBLE, dimids3, pop_varid_2), "pop var def error", error)
    Call check (nf90_def_var(ncid, PCP_ERROR_NAME_2, NF90_DOUBLE, dimids3, pcp_error_varid_2), "pcp_error var def error&
   &", error)
    If (error /= 0) Return
 
    Call check (nf90_def_var(ncid, tmean_name_2, NF90_DOUBLE, dimids3, tmean_varid_2), "tmean var def error", error)
    Call check (nf90_def_var(ncid, tmean_error_name_2, NF90_DOUBLE, dimids3, tmean_error_varid_2), "tmean error var def&
   & error", error)
    Call check (nf90_def_var(ncid, trange_name_2, NF90_DOUBLE, dimids3, trange_varid_2), "trange var def error", error)
    Call check (nf90_def_var(ncid, trange_error_name_2, NF90_DOUBLE, dimids3, trange_error_varid_2), "trange error var &
   &def error", error)
    If (error /= 0) Return
 
 
     !transformed mean, std variables, min, max of normalized y
    Call check (nf90_def_var(ncid, y_mean_name, NF90_DOUBLE, dimids3, ymean_varid), "y_mean var def error", error)
    Call check (nf90_def_var(ncid, y_std_name, NF90_DOUBLE, dimids3, ystd_varid), "y_std var def error", error)
    Call check (nf90_def_var(ncid, y_max_name, NF90_DOUBLE, dimids3, ymax_varid), "y_max var def error", error)
    Call check (nf90_def_var(ncid, y_min_name, NF90_DOUBLE, dimids3, ymin_varid), "y_min var def error", error)
    Call check (nf90_def_var(ncid, y_stdall_name, NF90_DOUBLE, dimids3, ystdall_varid), "y_std_all var def error", &
   & error)
 
    ! Add attributes.
 
     !long names
    Call check (nf90_put_att(ncid, pcp_varid, LONG_NAME, PCP_LONG_NAME), "pcp long_name attribute error", error)
    Call check (nf90_put_att(ncid, pop_varid, LONG_NAME, POP_LONG_NAME), "pcp long_name attribute error", error)
    Call check (nf90_put_att(ncid, pcp_error_varid, LONG_NAME, PCP_ERROR_LONG_NAME), "pcp_error long_name attribute err&
   &or", error)
 
    Call check (nf90_put_att(ncid, tmean_varid, LONG_NAME, tmean_long_name), "tmean long_name attribute error", error)
    Call check (nf90_put_att(ncid, tmean_error_varid, LONG_NAME, tmean_error_long_name), "tmean long_name attribute err&
   &or", error)
    Call check (nf90_put_att(ncid, trange_varid, LONG_NAME, trange_long_name), "trange long_name attribute error", &
   & error)
    Call check (nf90_put_att(ncid, trange_error_varid, LONG_NAME, trange_error_long_name), "trange long_name attribute &
   &error", error)
 
    Call check (nf90_put_att(ncid, pcp_varid_2, LONG_NAME, PCP_LONG_NAME_2), "pcp long_name attribute error", error)
    Call check (nf90_put_att(ncid, pop_varid_2, LONG_NAME, POP_LONG_NAME_2), "pcp long_name attribute error", error)
    Call check (nf90_put_att(ncid, pcp_error_varid_2, LONG_NAME, PCP_ERROR_LONG_NAME_2), "pcp_error long_name attribute&
   & error", error)
 
    Call check (nf90_put_att(ncid, tmean_varid_2, LONG_NAME, tmean_long_name_2), "tmean long_name attribute error", &
   & error)
    Call check (nf90_put_att(ncid, tmean_error_varid_2, LONG_NAME, tmean_error_long_name_2), "tmean long_name attribute&
   & error", error)
    Call check (nf90_put_att(ncid, trange_varid_2, LONG_NAME, trange_long_name_2), "trange long_name attribute error", &
   & error)
    Call check (nf90_put_att(ncid, trange_error_varid_2, LONG_NAME, trange_error_long_name_2), "trange long_name attrib&
   &ute error", error)
 
 
     !correlation variables
    Call check (nf90_put_att(ncid, autoc_varid, LONG_NAME, autoc_long_name), "auto_corr long_name attribute error", &
   & error)
    Call check (nf90_put_att(ncid, tpc_varid, LONG_NAME, tpc_long_name), "tp_corr long_name attribute error", error)
 
     !transformed mean, std variables, min, max of normalized y
    Call check (nf90_put_att(ncid, ymean_varid, LONG_NAME, y_mean_long_name), "ymean long_name attribute error", error)
    Call check (nf90_put_att(ncid, ystd_varid, LONG_NAME, y_std_long_name), "ystd long_name attribute error", error)
    Call check (nf90_put_att(ncid, ystdall_varid, LONG_NAME, y_std_long_name), "ystd_all long_name attribute error", &
   & error)
    Call check (nf90_put_att(ncid, ymax_varid, LONG_NAME, y_max_long_name), "ymax long_name attribute error", error)
    Call check (nf90_put_att(ncid, ymin_varid, LONG_NAME, y_min_long_name), "ymin long_name attribute error", error)
     !units
    Call check (nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS), "lat units attribute error", error)
    Call check (nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS), "lon units attribute error", error)
    Call check (nf90_put_att(ncid, alt_varid, UNITS, ALT_UNITS), "alt units attribute error", error)
    Call check (nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)
 
 
    Call check (nf90_put_att(ncid, pcp_varid, UNITS, PCP_UNITS), "pcp units attribute error", error)
    Call check (nf90_put_att(ncid, pop_varid, UNITS, POP_UNITS), "pcp units attribute error", error)
    Call check (nf90_put_att(ncid, pcp_error_varid, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)
 
    Call check (nf90_put_att(ncid, tmean_varid, UNITS, tmean_units), "tmean units attribute error", error)
    Call check (nf90_put_att(ncid, tmean_error_varid, UNITS, tmean_error_units), "tmean_error units attribute error", &
   & error)
    Call check (nf90_put_att(ncid, trange_varid, UNITS, trange_units), "trange units attribute error", error)
    Call check (nf90_put_att(ncid, trange_error_varid, UNITS, trange_error_units), "trange_error units attribute error",&
   &  error)
 
    Call check (nf90_put_att(ncid, pcp_varid_2, UNITS, PCP_UNITS), "pcp units attribute error", error)
    Call check (nf90_put_att(ncid, pop_varid_2, UNITS, POP_UNITS), "pcp units attribute error", error)
    Call check (nf90_put_att(ncid, pcp_error_varid_2, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", &
   & error)
 
    Call check (nf90_put_att(ncid, tmean_varid_2, UNITS, tmean_units), "tmean units attribute error", error)
    Call check (nf90_put_att(ncid, tmean_error_varid_2, UNITS, tmean_error_units), "tmean_error units attribute error", &
   & error)
    Call check (nf90_put_att(ncid, trange_varid_2, UNITS, trange_units), "trange units attribute error", error)
    Call check (nf90_put_att(ncid, trange_error_varid_2, UNITS, trange_error_units), "trange_error units attribute erro&
   &r", error)
 
 
     !correlation variables
    Call check (nf90_put_att(ncid, autoc_varid, UNITS, autoc_units), "auto correlation units attribute error", error)
    Call check (nf90_put_att(ncid, tpc_varid, UNITS, tpc_units), "tp correlation units attribute error", error)
 
     !transformed mean,std variables, min, max of normalized y
    Call check (nf90_put_att(ncid, ymean_varid, UNITS, y_mean_units), "ymean units attribute error", error)
    Call check (nf90_put_att(ncid, ystd_varid, UNITS, y_std_units), "ystd units attribute error", error)
    Call check (nf90_put_att(ncid, ystdall_varid, UNITS, y_std_units), "ystd_all units attribute error", error)
    Call check (nf90_put_att(ncid, ymax_varid, UNITS, y_max_units), "ymax units attribute error", error)
    Call check (nf90_put_att(ncid, ymin_varid, UNITS, y_min_units), "ymin units attribute error", error)
 
    If (error /= 0) Return
 
 
 
 
     ! End define mode.
    Call check (nf90_enddef(ncid), "end define mode error", error)
    If (error /= 0) Return
 
    count2 = (/ inx, iny /)
    start2 = (/ 1, 1 /)
 
    Call check (nf90_put_var(ncid, lat_varid, grdlat, start=start2, count=count2), "put lat error", error)
    Call check (nf90_put_var(ncid, lon_varid, grdlon, start=start2, count=count2), "put lon error", error)
    Call check (nf90_put_var(ncid, alt_varid, grdalt, start=start2, count=count2), "put alt error", error)
 
    trec = 1
    nrecs = 0
 
 
  Else
 
     ! File already exists, get dim and var ids
    Call check (nf90_inq_dimid(ncid, X_NAME, x_dimid), "x dim inq error", error)
    Call check (nf90_inq_dimid(ncid, Y_NAME, y_dimid), "y dim inq error", error)
    Call check (nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
    If (error /= 0) Return
 
    Call check (nf90_inq_varid(ncid, LAT_NAME, lat_varid), "lat var inq error", error)
    Call check (nf90_inq_varid(ncid, LON_NAME, lon_varid), "lon var inq error", error)
    Call check (nf90_inq_varid(ncid, TIME_NAME, time_varid), "time var inq error", error)
    Call check (nf90_inq_varid(ncid, PCP_NAME, pcp_varid), "pcp var inq error", error)
    Call check (nf90_inq_varid(ncid, POP_NAME, pop_varid), "pop var inq error", error)
    Call check (nf90_inq_varid(ncid, PCP_ERROR_NAME, pcp_error_varid), "pcp_error var inq error", error)
    Call check (nf90_inq_varid(ncid, tmean_name, tmean_varid), "tmean var inq error", error)
    Call check (nf90_inq_varid(ncid, tmean_error_name, tmean_error_varid), "tmean error var inq error", error)
    Call check (nf90_inq_varid(ncid, trange_name, trange_varid), "trange var inq error", error)
    Call check (nf90_inq_varid(ncid, trange_error_name, trange_error_varid), "trange error var inq error", error)
 
    Call check (nf90_inq_varid(ncid, PCP_NAME_2, pcp_varid_2), "pcp var inq error", error)
    Call check (nf90_inq_varid(ncid, POP_NAME_2, pop_varid_2), "pop var inq error", error)
    Call check (nf90_inq_varid(ncid, PCP_ERROR_NAME_2, pcp_error_varid_2), "pcp_error var inq error", error)
    Call check (nf90_inq_varid(ncid, tmean_name_2, tmean_varid_2), "tmean var inq error", error)
    Call check (nf90_inq_varid(ncid, tmean_error_name_2, tmean_error_varid_2), "tmean error var inq error", error)
    Call check (nf90_inq_varid(ncid, trange_name_2, trange_varid_2), "trange var inq error", error)
    Call check (nf90_inq_varid(ncid, trange_error_name_2, trange_error_varid_2), "trange error var inq error", error)
 
 
    Call check (nf90_inq_varid(ncid, autoc_name, autoc_varid), "autoc var inq error", error)
    Call check (nf90_inq_varid(ncid, tpc_name, tpc_varid), "tpc var inq error", error)
 
    Call check (nf90_inq_varid(ncid, y_mean_name, ymean_varid), "ymean var inq error", error)
    Call check (nf90_inq_varid(ncid, y_std_name, ystd_varid), "ystd var inq error", error)
    Call check (nf90_inq_varid(ncid, y_stdall_name, ystdall_varid), "ystd_all var inq error", error)
    Call check (nf90_inq_varid(ncid, y_max_name, ymax_varid), "ymax var inq error", error)
    Call check (nf90_inq_varid(ncid, y_min_name, ymin_varid), "ymin var inq error", error)
    If (error /= 0) Return
 
    Call check (nf90_Inquire_Dimension(ncid, x_dimid, len=file_nx), "x dim len error", error)
    Call check (nf90_Inquire_Dimension(ncid, y_dimid, len=file_ny), "y dim len error", error)
    Call check (nf90_Inquire_Dimension(ncid, time_dimid, len=file_ntimes), "time dim len error", error)
    If (error /= 0) Return
 
    If (nx /= file_nx .Or. ny /= file_ny) Then
      Print *, "Error dimensions in output file do not match current run."
      error = 1
      Return
    End If
 
    Allocate (file_times(file_ntimes))
    Call check (nf90_get_var(ncid, time_varid, file_times), "error getting file times list", error)
    If (error /= 0) Return
 
    If (file_times(1) > times(n_times)) Then !put data before everything in the file
      Print *, "Error cannot add data before data already in output file. (functionality still to be added)"
      error = 1
      Return
    Else
      If (file_times(file_ntimes) < times(1)) Then !put data after everything in the file
        trec = file_ntimes + 1
      Else ! at least some overlap
        Do i = 1, file_ntimes, 1
          If (file_times(1) == times(1)) Then
            trec = i
          End If
        End Do
        If (trec == 0) Then
          Print *, "Error, confusion over data output record location."
          error = 1
          Return
        Else
          Print *, "WARNING, overwriting data in output file, record ", trec, " to ", trec + n_times - 1
        End If
      End If
    End If
 
  End If
 
  count1 (1) = n_times
  start1 (1) = trec
  Call check (nf90_put_var(ncid, time_varid, times, start=start1, count=count1), "put times error", error)
  If (error /= 0) Return
 
 
  !correlation variables
  Call check (nf90_put_var(ncid, autoc_varid, mean_autocorr, start=start1, count=count1), "put mean autocorrelation err&
 &or", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, tpc_varid, mean_tp_corr, start=start1, count=count1), "put mean t_p correlation error",&
 &  error)
  If (error /= 0) Return
 
 
 
 
  !3-d variables
  count3 = (/ inx, iny, n_times /)
  start3 = (/ 1, 1, trec /)
  Call check (nf90_put_var(ncid, pcp_varid, real(pcp, kind(DP)), start=start3, count=count3), "put pcp error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, pop_varid, real(pop, kind(DP)), start=start3, count=count3), "put pop error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, pcp_error_varid, real(pcperror, kind(DP)), start=start3, count=count3), "put pcp_error&
 & error", error)
  If (error /= 0) Return
 
 
  Call check (nf90_put_var(ncid, tmean_varid, real(tmean, kind(DP)), start=start3, count=count3), "put tmean error", &
 & error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, tmean_error_varid, real(tmean_error, kind(DP)), start=start3, count=count3), "put tmea&
 &n_error error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, trange_varid, real(trange, kind(DP)), start=start3, count=count3), "put trange error", &
 & error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, trange_error_varid, real(trange_error, kind(DP)), start=start3, count=count3), "put tr&
 &ange_error error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, pcp_varid_2, real(pcp_2, kind(DP)), start=start3, count=count3), "put pcp error", &
 & error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, pop_varid_2, real(pop_2, kind(DP)), start=start3, count=count3), "put pop error", &
 & error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, pcp_error_varid_2, real(pcperror_2, kind(DP)), start=start3, count=count3), "put pcp_e&
 &rror error", error)
  If (error /= 0) Return
 
 
  Call check (nf90_put_var(ncid, tmean_varid_2, real(tmean_2, kind(DP)), start=start3, count=count3), "put tmean error",&
 &  error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, tmean_error_varid_2, real(tmean_error_2, kind(DP)), start=start3, count=count3), "put &
 &tmean_error error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, trange_varid_2, real(trange_2, kind(DP)), start=start3, count=count3), "put trange err&
 &or", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, trange_error_varid_2, real(trange_error_2, kind(DP)), start=start3, count=count3), "pu&
 &t trange_error error", error)
  If (error /= 0) Return
 
 
 
!transformed mean,std variables, min & max of normalized y
  Call check (nf90_put_var(ncid, ymean_varid, y_mean, start=start3, count=count3), "put ymean error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, ystd_varid, y_std, start=start3, count=count3), "put ystd error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, ystdall_varid, y_std_all, start=start3, count=count3), "put ystd_all error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, ymax_varid, y_min, start=start3, count=count3), "put ymax error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, ymin_varid, y_max, start=start3, count=count3), "put ymin error", error)
  If (error /= 0) Return
 
  Call check (nf90_close(ncid), "closing file error", error)
 
 
Contains
  Subroutine check (status, info, error)
    Integer, Intent (In) :: status
    Character (Len=*), Intent (In) :: info
    Integer, Intent (Out) :: error
 
    If (status /= nf90_noerr) Then
      Print *, trim (info) // ": " // trim (nf90_strerror(status))
      error = 1
    End If
  End Subroutine check
End Subroutine save_precip
