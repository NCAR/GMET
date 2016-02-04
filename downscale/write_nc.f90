 
SUBROUTINE save_coefficients (n_vars, var_names, coefs, startdate, enddate, times, site_list, station_var, stnid, &
& stnlat, stnlon, stnalt, forecast, file, error)
  USE netcdf
  USE type
  IMPLICIT NONE
 
  CHARACTER (LEN=100), INTENT (IN) :: var_names (:), stnid (:)
  INTEGER, INTENT (IN) :: n_vars, forecast
  CHARACTER (LEN=100), INTENT (IN) :: startdate, enddate, station_var
  CHARACTER (LEN=500), INTENT (IN) :: file, site_list
  REAL (DP), INTENT (IN) :: stnlat (:), stnlon (:), stnalt (:)
  REAL (DP), INTENT (IN) :: coefs (:, :, :)
  REAL (DP), INTENT (IN) :: times (:)
 
  INTEGER, INTENT (OUT) :: error
 
  ! Dimension names
  CHARACTER (LEN=*), PARAMETER :: STN_NAME = "station"
  CHARACTER (LEN=*), PARAMETER :: VAR_NAME = "variable"
  CHARACTER (LEN=*), PARAMETER :: CHAR_NAME = "string"
  CHARACTER (LEN=*), PARAMETER :: TIME_NAME = "time"
  CHARACTER (LEN=*), PARAMETER :: REC_NAME = "run"
  ! Variable Names
  CHARACTER (LEN=*), PARAMETER :: VARS_NAME = "variable_name"
  CHARACTER (LEN=*), PARAMETER :: STN_FILE_NAME = "station_file_name"
  CHARACTER (LEN=*), PARAMETER :: STN_ID_NAME = "station_id"
  CHARACTER (LEN=*), PARAMETER :: STN_LAT_NAME = "station_latitude"
  CHARACTER (LEN=*), PARAMETER :: STN_LON_NAME = "station_longitude"
  CHARACTER (LEN=*), PARAMETER :: STN_ALT_NAME = "station_altitude"
  CHARACTER (LEN=*), PARAMETER :: FORECAST_NAME = "forecast_hr"
  CHARACTER (LEN=*), PARAMETER :: STARTDATE_NAME = "run_start_date"
  CHARACTER (LEN=*), PARAMETER :: ENDDATE_NAME = "run_end_date"
  CHARACTER (LEN=*), PARAMETER :: TIME_VAR_NAME = "time"
  CHARACTER (LEN=*), PARAMETER :: COEFS_NAME = "coefficient"
  CHARACTER (LEN=*), PARAMETER :: CONSTANT_NAME = "constant_term"
  ! Units
  CHARACTER (LEN=*), PARAMETER :: UNITS = "units"
  CHARACTER (LEN=*), PARAMETER :: COEFS_UNITS = "linear_equation_coefficient"
  CHARACTER (LEN=*), PARAMETER :: LAT_UNITS = "degrees_north"
  CHARACTER (LEN=*), PARAMETER :: LON_UNITS = "degrees_east"
  CHARACTER (LEN=*), PARAMETER :: ALT_UNITS = "feet"
  CHARACTER (LEN=*), PARAMETER :: FORECAST_UNITS = "hours"
  CHARACTER (LEN=*), PARAMETER :: DATE_UNITS = "YYYYMMDD"
  CHARACTER (LEN=*), PARAMETER :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  CHARACTER (LEN=*), PARAMETER :: FILL = "_FillValue"
 
  INTEGER :: n_stns, n_chars, n_times
  INTEGER :: ncid, stn_dimid, vars_dimid, char_dimid, time_dimid, rec_dimid
  INTEGER :: coefs_varid, time_varid, vars_varid, forecast_varid, startdate_varid, enddate_varid
  INTEGER :: stn_id_varid, stn_lat_varid, stn_lon_varid, stn_alt_varid
  INTEGER :: charids (2)
  INTEGER :: dimids (4)
  INTEGER :: count4 (4), start4 (4), count3 (3), start3 (3), count2 (2), start2 (2), count1 (1), start1 (1)
  INTEGER :: nstns, nvars, nrecs, ntimes, rec, i
 
  CHARACTER (LEN=100) :: startdates, enddates
  INTEGER, ALLOCATABLE :: forecasts (:)
 
  INTEGER :: forecast_arr (1)
 
  n_chars = 100
  n_stns = size (stnlat)
  n_times = size (times)
 
  error = nf90_open (file, nf90_write, ncid)
  IF (error /= nf90_noerr) THEN
   error = 0
     ! Create the file.
   CALL check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
   IF (error /= 0) RETURN
 
     ! Define the dimensions.
   CALL check (nf90_def_dim(ncid, STN_NAME, n_stns, stn_dimid), "station dim def error", error)
   CALL check (nf90_def_dim(ncid, VAR_NAME, n_vars+1, vars_dimid), "variable dim def error", error)
   CALL check (nf90_def_dim(ncid, CHAR_NAME, n_chars, char_dimid), "char dim def error", error)
   CALL check (nf90_def_dim(ncid, TIME_NAME, n_times, time_dimid), "time dim def error", error)
   CALL check (nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid), "rec dim def error", error)
   IF (error /= 0) RETURN
 
     ! Define the variables.
   charids = (/ char_dimid, vars_dimid /)
   CALL check (nf90_def_var(ncid, VARS_NAME, NF90_CHAR, charids, vars_varid), "variable_name var def error", error)
   charids = (/ char_dimid, stn_dimid /)
   CALL check (nf90_def_var(ncid, STN_ID_NAME, NF90_CHAR, charids, stn_id_varid), "station_id var def error", error)
   CALL check (nf90_def_var(ncid, STN_LAT_NAME, NF90_DOUBLE, stn_dimid, stn_lat_varid), "station_latitude var def error&
  &", error)
   CALL check (nf90_def_var(ncid, STN_LON_NAME, NF90_DOUBLE, stn_dimid, stn_lon_varid), "station_longitude var def erro&
  &r", error)
   CALL check (nf90_def_var(ncid, STN_ALT_NAME, NF90_DOUBLE, stn_dimid, stn_alt_varid), "station_altitude var def error&
  &", error)
   IF (error /= 0) RETURN
 
   charids = (/ char_dimid, rec_dimid /)
   CALL check (nf90_def_var(ncid, STARTDATE_NAME, NF90_CHAR, charids, startdate_varid), "start_date var def error", &
  & error)
   CALL check (nf90_def_var(ncid, ENDDATE_NAME, NF90_CHAR, charids, enddate_varid), "end_date var def error", error)
   CALL check (nf90_def_var(ncid, FORECAST_NAME, NF90_INT, rec_dimid, forecast_varid), "forecast var def error", error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_def_var(ncid, TIME_VAR_NAME, NF90_DOUBLE, time_dimid, time_varid), "time var def error", error)
   dimids = (/ stn_dimid, time_dimid, vars_dimid, rec_dimid /)
   CALL check (nf90_def_var(ncid, COEFS_NAME, NF90_DOUBLE, dimids, coefs_varid), "coefficient var def error", error)
   IF (error /= 0) RETURN
 
     ! Add attributes.
   CALL check (nf90_put_att(ncid, stn_lat_varid, UNITS, LAT_UNITS), "station_lat units attribute error", error)
   CALL check (nf90_put_att(ncid, stn_lon_varid, UNITS, LON_UNITS), "station_lon units attribute error", error)
   CALL check (nf90_put_att(ncid, stn_alt_varid, UNITS, ALT_UNITS), "station_alt units attribute error", error)
   CALL check (nf90_put_att(ncid, startdate_varid, UNITS, DATE_UNITS), "start_date units attribute error", error)
   CALL check (nf90_put_att(ncid, enddate_varid, UNITS, DATE_UNITS), "end_date units attribute error", error)
   CALL check (nf90_put_att(ncid, forecast_varid, UNITS, FORECAST_UNITS), "forecase units attribute error", error)
   CALL check (nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)
   CALL check (nf90_put_att(ncid, coefs_varid, UNITS, COEFS_UNITS), "coefficient units attribute error", error)
   CALL check (nf90_put_att(ncid, NF90_GLOBAL, STN_FILE_NAME, site_list), "station_file_name global attribute error", &
  & error)
 
     ! End define mode.
   CALL check (nf90_enddef(ncid), "end define mode error", error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_put_var(ncid, stn_lat_varid, stnlat), "put staion_lat error", error)
   CALL check (nf90_put_var(ncid, stn_lon_varid, stnlon), "put staion_lon error", error)
   CALL check (nf90_put_var(ncid, stn_alt_varid, stnalt), "put staion_alt error", error)
   CALL check (nf90_put_var(ncid, time_varid, times), "put times error", error)
 
   count2 = (/ 1, 1 /)
   start2 = (/ 1, 1 /)
   DO i = 1, n_stns, 1
    count2 (1) = len (trim(stnid(i)))
    start2 (2) = i
    CALL check (nf90_put_var(ncid, stn_id_varid, stnid(i), start=start2, count=count2), "put staion_id error", error)
    IF (error /= 0) RETURN
   END DO
 
   count2 (1) = len (trim(CONSTANT_NAME))
   start2 (2) = 1
   CALL check (nf90_put_var(ncid, vars_varid, CONSTANT_NAME, start=start2, count=count2), "put variable_name error", &
  & error)
   IF (error /= 0) RETURN
   DO i = 1, n_vars, 1
    count2 (1) = len (trim(var_names(i)))
    start2 (2) = i + 1
    CALL check (nf90_put_var(ncid, vars_varid, var_names(i), start=start2, count=count2), "put variable_name error", &
   & error)
    IF (error /= 0) RETURN
   END DO
 
   rec = 1
   nrecs = 0
 
  ELSE
 
     ! File already exists, get dim and var ids
   CALL check (nf90_inq_dimid(ncid, STN_NAME, stn_dimid), "station dim inq error", error)
   CALL check (nf90_inq_dimid(ncid, VAR_NAME, vars_dimid), "variable dim inq error", error)
   CALL check (nf90_inq_dimid(ncid, CHAR_NAME, char_dimid), "char dim inq error", error)
   CALL check (nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
   CALL check (nf90_inq_dimid(ncid, REC_NAME, rec_dimid), "run dim inq error", error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_inq_varid(ncid, VARS_NAME, vars_varid), "variable_name var inq error", error)
   CALL check (nf90_inq_varid(ncid, STN_ID_NAME, stn_id_varid), "station_id var inq error", error)
   CALL check (nf90_inq_varid(ncid, STN_LAT_NAME, stn_lat_varid), "station_latitude var inq error", error)
   CALL check (nf90_inq_varid(ncid, STN_LON_NAME, stn_lon_varid), "station_longitude var inq error", error)
   CALL check (nf90_inq_varid(ncid, STN_ALT_NAME, stn_alt_varid), "station_altitude var inq error", error)
   CALL check (nf90_inq_varid(ncid, STARTDATE_NAME, startdate_varid), "start_date var inq error", error)
   CALL check (nf90_inq_varid(ncid, ENDDATE_NAME, enddate_varid), "end_date var inq error", error)
   CALL check (nf90_inq_varid(ncid, FORECAST_NAME, forecast_varid), "forecast var inq error", error)
   CALL check (nf90_inq_varid(ncid, TIME_VAR_NAME, time_varid), "time var inq error", error)
   CALL check (nf90_inq_varid(ncid, COEFS_NAME, coefs_varid), "coefficient var inq error", error)
   IF (error /= 0) RETURN
 
     ! Verify Dimensions match
   CALL check (nf90_Inquire_Dimension(ncid, stn_dimid, len=nstns), "station dim len error", error)
   CALL check (nf90_Inquire_Dimension(ncid, vars_dimid, len=nvars), "variable dim len error", error)
   CALL check (nf90_Inquire_Dimension(ncid, time_dimid, len=ntimes), "time dim len error", error)
   CALL check (nf90_Inquire_Dimension(ncid, rec_dimid, len=nrecs), "run dim len error", error)
   IF (error /= 0) RETURN
 
   IF (n_stns /= nstns .OR. nvars /= n_vars+1 .OR. ntimes /= n_times) THEN
    PRINT *, "Error dimensions in output file do not match current run."
    error = 1
    RETURN
   END IF
 
   ALLOCATE (forecasts(nrecs))
 
   CALL check (nf90_get_var(ncid, forecast_varid, forecasts), "error getting file forecasts", error)
   IF (error /= 0) RETURN
 
   rec = nrecs + 1
   DO i = 1, nrecs, 1
    count2 = (/ 100, 1 /)
    start2 = (/ 1, i /)
    CALL check (nf90_get_var(ncid, startdate_varid, startdates, start=start2, count=count2), "error getting file startd&
   &ates", error)
    CALL check (nf90_get_var(ncid, enddate_varid, enddates, start=start2, count=count2), "error getting file enddates", &
   & error)
    IF (forecasts(i) == forecast .AND. startdates(1:8) == startdate(1:8) .AND. enddates(1:8) == enddate(1:8)) THEN
     PRINT *, "WARNING, overwriting data in output file, record ", i
     rec = i
    END IF
   END DO
 
 
  END IF
 
  count2 = (/ len (trim(startdate)), 1 /)
  start2 = (/ 1, rec /)
  CALL check (nf90_put_var(ncid, startdate_varid, startdate, start=start2, count=count2), "put start_date error", &
 & error)
  IF (error /= 0) RETURN
 
  count2 = (/ len (trim(enddate)), 1 /)
  start2 = (/ 1, rec /)
  CALL check (nf90_put_var(ncid, enddate_varid, enddate, start=start2, count=count2), "put start_date error", error)
  IF (error /= 0) RETURN
 
  count1 (1) = 1
  start1 (1) = rec
  forecast_arr = forecast
  CALL check (nf90_put_var(ncid, forecast_varid, forecast_arr, start=start1, count=count1), "put forecast error", &
 & error)
  IF (error /= 0) RETURN
 
  count4 = (/ n_stns, n_times, n_vars + 1, 1 /)
  start4 = (/ 1, 1, 1, rec /)
  forecast_arr = forecast
  CALL check (nf90_put_var(ncid, coefs_varid, coefs, start=start4, count=count4), "put coefficients error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_close(ncid), "closing file error", error)
 
CONTAINS
  SUBROUTINE check (status, info, error)
   INTEGER, INTENT (IN) :: status
   CHARACTER (LEN=*), INTENT (IN) :: info
   INTEGER, INTENT (OUT) :: error
 
   IF (status /= nf90_noerr) THEN
    PRINT *, trim (info) // ": " // trim (nf90_strerror(status))
    error = 1
   END IF
  END SUBROUTINE check
END SUBROUTINE save_coefficients
 
SUBROUTINE save_precip (pcp, pop, pcperror, tmean, tmean_error, trange, trange_error, nx, ny, grdlat, grdlon, grdalt, &
& times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, file, error, pcp_2, pop_2, pcperror_2, &
& tmean_2, tmean_error_2, trange_2, trange_error_2)
  USE netcdf
  USE type
  IMPLICIT NONE
 
  REAL (SP), INTENT (IN) :: pcp (:, :), pop (:, :), pcperror (:, :)
  REAL (SP), INTENT (IN) :: tmean (:, :), tmean_error (:, :), trange (:, :), trange_error (:, :)
 
  REAL (SP), INTENT (IN) :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
  REAL (SP), INTENT (IN) :: tmean_2 (:, :), tmean_error_2 (:, :), trange_2 (:, :), trange_error_2 (:, :)
 
  INTEGER (I4B), INTENT (IN) :: nx, ny
  REAL (DP), INTENT (IN) :: grdlat (:), grdlon (:), grdalt (:)
  REAL (DP), INTENT (IN) :: times (:)
  REAL (DP), INTENT (IN) :: mean_autocorr (:), mean_tp_corr (:)
 
  REAL (DP), INTENT (IN) :: y_mean (:, :), y_std (:, :), y_min (:, :), y_max (:, :), y_std_all (:, :)
  CHARACTER (LEN=500), INTENT (IN) :: file
  INTEGER, INTENT (OUT) :: error
 
 
  ! Dimension names
  CHARACTER (LEN=*), PARAMETER :: Y_NAME = "y"
  CHARACTER (LEN=*), PARAMETER :: X_NAME = "x"
  CHARACTER (LEN=*), PARAMETER :: TIME_NAME = "time"
 
  ! Variable Names
  CHARACTER (LEN=*), PARAMETER :: LAT_NAME = "latitude"
  CHARACTER (LEN=*), PARAMETER :: LON_NAME = "longitude"
  CHARACTER (LEN=*), PARAMETER :: autoc_name = "auto_corr"
  CHARACTER (LEN=*), PARAMETER :: tpc_name = "tp_corr"
  CHARACTER (LEN=*), PARAMETER :: ALT_NAME = "altitude"
  CHARACTER (LEN=*), PARAMETER :: PCP_NAME = "pcp"
  CHARACTER (LEN=*), PARAMETER :: POP_NAME = "pop"
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_NAME = "pcp_error"
  CHARACTER (LEN=*), PARAMETER :: tmean_name = "tmean"
  CHARACTER (LEN=*), PARAMETER :: tmean_error_name = "tmean_error"
  CHARACTER (LEN=*), PARAMETER :: trange_name = "trange"
  CHARACTER (LEN=*), PARAMETER :: trange_error_name = "trange_error"
  CHARACTER (LEN=*), PARAMETER :: y_mean_name = "ymean"
  CHARACTER (LEN=*), PARAMETER :: y_std_name = "ystd"
  CHARACTER (LEN=*), PARAMETER :: y_stdall_name = "ystd_all"
  CHARACTER (LEN=*), PARAMETER :: y_min_name = "ymax"
  CHARACTER (LEN=*), PARAMETER :: y_max_name = "ymin"
  CHARACTER (LEN=*), PARAMETER :: PCP_NAME_2 = "pcp_2"
  CHARACTER (LEN=*), PARAMETER :: POP_NAME_2 = "pop_2"
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_NAME_2 = "pcp_error_2"
  CHARACTER (LEN=*), PARAMETER :: tmean_name_2 = "tmean_2"
  CHARACTER (LEN=*), PARAMETER :: tmean_error_name_2 = "tmean_error_2"
  CHARACTER (LEN=*), PARAMETER :: trange_name_2 = "trange_2"
  CHARACTER (LEN=*), PARAMETER :: trange_error_name_2 = "trange_error_2"
 
 
  CHARACTER (LEN=*), PARAMETER :: LONG_NAME = "long_name"
  CHARACTER (LEN=*), PARAMETER :: PCP_LONG_NAME = "estimated precip in normal space"
  CHARACTER (LEN=*), PARAMETER :: POP_LONG_NAME = "probability of precipitation occurrence"
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_LONG_NAME = "error in estimated precip"
  CHARACTER (LEN=*), PARAMETER :: tmean_long_name = "estimated daily mean temperature"
  CHARACTER (LEN=*), PARAMETER :: tmean_error_long_name = "error in estimated daily mean temp"
  CHARACTER (LEN=*), PARAMETER :: trange_long_name = "estimated diurnal range"
  CHARACTER (LEN=*), PARAMETER :: trange_error_long_name = "error in estimated diurnal range"
  CHARACTER (LEN=*), PARAMETER :: autoc_long_name = "Lag-1 autocorrelation of temperature"
  CHARACTER (LEN=*), PARAMETER :: tpc_long_name = "Correlation of diurnal range and precipitation"
  CHARACTER (LEN=*), PARAMETER :: y_mean_long_name = "mean of transformed non-zero precip"
  CHARACTER (LEN=*), PARAMETER :: y_std_long_name = "std. dev. of transformed non-zero precip"
  CHARACTER (LEN=*), PARAMETER :: y_max_long_name = "max of normalized transformed non-zero precip"
  CHARACTER (LEN=*), PARAMETER :: y_min_long_name = "min of normalized transformed non-zero precip"
  CHARACTER (LEN=*), PARAMETER :: PCP_LONG_NAME_2 = "estimated precip in normal space (no slope)"
  CHARACTER (LEN=*), PARAMETER :: POP_LONG_NAME_2 = "probability of precipitation occurrence (no slope)"
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_LONG_NAME_2 = "error in estimated precip (no slope)"
  CHARACTER (LEN=*), PARAMETER :: tmean_long_name_2 = "estimated daily mean temperature (no slope)"
  CHARACTER (LEN=*), PARAMETER :: tmean_error_long_name_2 = "error in estimated daily mean temp (no slope)"
  CHARACTER (LEN=*), PARAMETER :: trange_long_name_2 = "estimated diurnal range (no slope)"
  CHARACTER (LEN=*), PARAMETER :: trange_error_long_name_2 = "error in estimated diurnal range (no slope)"
 
 
  ! Units
  CHARACTER (LEN=*), PARAMETER :: UNITS = "units"
  CHARACTER (LEN=*), PARAMETER :: PCP_UNITS = ""
  CHARACTER (LEN=*), PARAMETER :: POP_UNITS = ""
  CHARACTER (LEN=*), PARAMETER :: autoc_units = ""
  CHARACTER (LEN=*), PARAMETER :: tpc_units = ""
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_UNITS = ""
  CHARACTER (LEN=*), PARAMETER :: tmean_units = "deg_C"
  CHARACTER (LEN=*), PARAMETER :: trange_units = "deg_C"
  CHARACTER (LEN=*), PARAMETER :: tmean_error_units = "deg_C"
  CHARACTER (LEN=*), PARAMETER :: trange_error_units = "deg_C"
  CHARACTER (LEN=*), PARAMETER :: y_mean_units = ""
  CHARACTER (LEN=*), PARAMETER :: y_std_units = ""
  CHARACTER (LEN=*), PARAMETER :: y_max_units = ""
  CHARACTER (LEN=*), PARAMETER :: y_min_units = ""
 
  CHARACTER (LEN=*), PARAMETER :: LAT_UNITS = "degrees_north"
  CHARACTER (LEN=*), PARAMETER :: LON_UNITS = "degrees_east"
  CHARACTER (LEN=*), PARAMETER :: ALT_UNITS = "meters"
  CHARACTER (LEN=*), PARAMETER :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  CHARACTER (LEN=*), PARAMETER :: FILL = "_FillValue"
 
  REAL (DP), ALLOCATABLE :: file_times (:)
 
  INTEGER :: n_chars, n_times, inx, iny
  INTEGER :: ncid, x_dimid, y_dimid, time_dimid
  INTEGER :: lat_varid, lon_varid, autoc_varid, alt_varid, time_varid, pcp_varid, pop_varid, pcp_error_varid, tpc_varid
  INTEGER :: tmean_varid, tmean_error_varid, trange_varid, trange_error_varid
  INTEGER :: ymean_varid, ystd_varid, ymax_varid, ymin_varid, ystdall_varid
  INTEGER :: count1 (1), start1 (1), count2 (2), start2 (2), count3 (3), start3 (3), dimids2 (2), dimids3 (3)
  INTEGER :: trec, nrecs, file_nx, file_ny, file_ntimes, i
 
  INTEGER :: pcp_varid_2, pop_varid_2, pcp_error_varid_2
  INTEGER :: tmean_varid_2, tmean_error_varid_2, trange_varid_2, trange_error_varid_2
 
  trec = 0
  n_chars = 100
  n_times = size (times)
  inx = nx
  iny = ny
 
  IF (size(grdlat) /= inx*iny) THEN
   PRINT *, "Error "
  END IF
 
  error = nf90_open (file, nf90_write, ncid)
  IF (error /= nf90_noerr) THEN
   error = 0
     ! Create the file.
   CALL check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
   IF (error /= 0) RETURN
 
     ! Define the dimensions.
   CALL check (nf90_def_dim(ncid, Y_NAME, iny, y_dimid), "y dim def error", error)
   CALL check (nf90_def_dim(ncid, X_NAME, inx, x_dimid), "x dim def error", error)
   CALL check (nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), "time dim def error", error)
   IF (error /= 0) RETURN
 
     ! Define the variables.
   dimids2 = (/ x_dimid, y_dimid /)
   CALL check (nf90_def_var(ncid, LAT_NAME, NF90_DOUBLE, dimids2, lat_varid), "lat var def error", error)
   CALL check (nf90_def_var(ncid, LON_NAME, NF90_DOUBLE, dimids2, lon_varid), "lon var def error", error)
   CALL check (nf90_def_var(ncid, ALT_NAME, NF90_DOUBLE, dimids2, alt_varid), "lon var def error", error)
   CALL check (nf90_def_var(ncid, TIME_NAME, NF90_DOUBLE, time_dimid, time_varid), "time var def error", error)
 
     !correlation variables
   CALL check (nf90_def_var(ncid, autoc_name, NF90_DOUBLE, time_dimid, autoc_varid), "auto correlaion var def error", &
  & error)
   CALL check (nf90_def_var(ncid, tpc_name, NF90_DOUBLE, time_dimid, tpc_varid), "tp correlaion var def error", error)
 
 
 
   IF (error /= 0) RETURN
 
   dimids3 = (/ x_dimid, y_dimid, time_dimid /)
   CALL check (nf90_def_var(ncid, PCP_NAME, NF90_DOUBLE, dimids3, pcp_varid), "pcp var def error", error)
   CALL check (nf90_def_var(ncid, POP_NAME, NF90_DOUBLE, dimids3, pop_varid), "pop var def error", error)
   CALL check (nf90_def_var(ncid, PCP_ERROR_NAME, NF90_DOUBLE, dimids3, pcp_error_varid), "pcp_error var def error", &
  & error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_def_var(ncid, tmean_name, NF90_DOUBLE, dimids3, tmean_varid), "tmean var def error", error)
   CALL check (nf90_def_var(ncid, tmean_error_name, NF90_DOUBLE, dimids3, tmean_error_varid), "tmean error var def erro&
  &r", error)
   CALL check (nf90_def_var(ncid, trange_name, NF90_DOUBLE, dimids3, trange_varid), "trange var def error", error)
   CALL check (nf90_def_var(ncid, trange_error_name, NF90_DOUBLE, dimids3, trange_error_varid), "trange error var def e&
  &rror", error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_def_var(ncid, PCP_NAME_2, NF90_DOUBLE, dimids3, pcp_varid_2), "pcp var def error", error)
   CALL check (nf90_def_var(ncid, POP_NAME_2, NF90_DOUBLE, dimids3, pop_varid_2), "pop var def error", error)
   CALL check (nf90_def_var(ncid, PCP_ERROR_NAME_2, NF90_DOUBLE, dimids3, pcp_error_varid_2), "pcp_error var def error",&
  &  error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_def_var(ncid, tmean_name_2, NF90_DOUBLE, dimids3, tmean_varid_2), "tmean var def error", error)
   CALL check (nf90_def_var(ncid, tmean_error_name_2, NF90_DOUBLE, dimids3, tmean_error_varid_2), "tmean error var def &
  &error", error)
   CALL check (nf90_def_var(ncid, trange_name_2, NF90_DOUBLE, dimids3, trange_varid_2), "trange var def error", error)
   CALL check (nf90_def_var(ncid, trange_error_name_2, NF90_DOUBLE, dimids3, trange_error_varid_2), "trange error var d&
  &ef error", error)
   IF (error /= 0) RETURN
 
 
     !transformed mean, std variables, min, max of normalized y
   CALL check (nf90_def_var(ncid, y_mean_name, NF90_DOUBLE, dimids3, ymean_varid), "y_mean var def error", error)
   CALL check (nf90_def_var(ncid, y_std_name, NF90_DOUBLE, dimids3, ystd_varid), "y_std var def error", error)
   CALL check (nf90_def_var(ncid, y_max_name, NF90_DOUBLE, dimids3, ymax_varid), "y_max var def error", error)
   CALL check (nf90_def_var(ncid, y_min_name, NF90_DOUBLE, dimids3, ymin_varid), "y_min var def error", error)
   CALL check (nf90_def_var(ncid, y_stdall_name, NF90_DOUBLE, dimids3, ystdall_varid), "y_std_all var def error", &
  & error)
 
    ! Add attributes.
 
     !long names
   CALL check (nf90_put_att(ncid, pcp_varid, LONG_NAME, PCP_LONG_NAME), "pcp long_name attribute error", error)
   CALL check (nf90_put_att(ncid, pop_varid, LONG_NAME, POP_LONG_NAME), "pcp long_name attribute error", error)
   CALL check (nf90_put_att(ncid, pcp_error_varid, LONG_NAME, PCP_ERROR_LONG_NAME), "pcp_error long_name attribute erro&
  &r", error)
 
   CALL check (nf90_put_att(ncid, tmean_varid, LONG_NAME, tmean_long_name), "tmean long_name attribute error", error)
   CALL check (nf90_put_att(ncid, tmean_error_varid, LONG_NAME, tmean_error_long_name), "tmean long_name attribute erro&
  &r", error)
   CALL check (nf90_put_att(ncid, trange_varid, LONG_NAME, trange_long_name), "trange long_name attribute error", &
  & error)
   CALL check (nf90_put_att(ncid, trange_error_varid, LONG_NAME, trange_error_long_name), "trange long_name attribute e&
  &rror", error)
 
   CALL check (nf90_put_att(ncid, pcp_varid_2, LONG_NAME, PCP_LONG_NAME_2), "pcp long_name attribute error", error)
   CALL check (nf90_put_att(ncid, pop_varid_2, LONG_NAME, POP_LONG_NAME_2), "pcp long_name attribute error", error)
   CALL check (nf90_put_att(ncid, pcp_error_varid_2, LONG_NAME, PCP_ERROR_LONG_NAME_2), "pcp_error long_name attribute &
  &error", error)
 
   CALL check (nf90_put_att(ncid, tmean_varid_2, LONG_NAME, tmean_long_name_2), "tmean long_name attribute error", &
  & error)
   CALL check (nf90_put_att(ncid, tmean_error_varid_2, LONG_NAME, tmean_error_long_name_2), "tmean long_name attribute &
  &error", error)
   CALL check (nf90_put_att(ncid, trange_varid_2, LONG_NAME, trange_long_name_2), "trange long_name attribute error", &
  & error)
   CALL check (nf90_put_att(ncid, trange_error_varid_2, LONG_NAME, trange_error_long_name_2), "trange long_name attribu&
  &te error", error)
 
 
     !correlation variables
   CALL check (nf90_put_att(ncid, autoc_varid, LONG_NAME, autoc_long_name), "auto_corr long_name attribute error", &
  & error)
   CALL check (nf90_put_att(ncid, tpc_varid, LONG_NAME, tpc_long_name), "tp_corr long_name attribute error", error)
 
     !transformed mean, std variables, min, max of normalized y
   CALL check (nf90_put_att(ncid, ymean_varid, LONG_NAME, y_mean_long_name), "ymean long_name attribute error", error)
   CALL check (nf90_put_att(ncid, ystd_varid, LONG_NAME, y_std_long_name), "ystd long_name attribute error", error)
   CALL check (nf90_put_att(ncid, ystdall_varid, LONG_NAME, y_std_long_name), "ystd_all long_name attribute error", &
  & error)
   CALL check (nf90_put_att(ncid, ymax_varid, LONG_NAME, y_max_long_name), "ymax long_name attribute error", error)
   CALL check (nf90_put_att(ncid, ymin_varid, LONG_NAME, y_min_long_name), "ymin long_name attribute error", error)
     !units
   CALL check (nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS), "lat units attribute error", error)
   CALL check (nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS), "lon units attribute error", error)
   CALL check (nf90_put_att(ncid, alt_varid, UNITS, ALT_UNITS), "alt units attribute error", error)
   CALL check (nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)
 
 
   CALL check (nf90_put_att(ncid, pcp_varid, UNITS, PCP_UNITS), "pcp units attribute error", error)
   CALL check (nf90_put_att(ncid, pop_varid, UNITS, POP_UNITS), "pcp units attribute error", error)
   CALL check (nf90_put_att(ncid, pcp_error_varid, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)
 
   CALL check (nf90_put_att(ncid, tmean_varid, UNITS, tmean_units), "tmean units attribute error", error)
   CALL check (nf90_put_att(ncid, tmean_error_varid, UNITS, tmean_error_units), "tmean_error units attribute error", &
  & error)
   CALL check (nf90_put_att(ncid, trange_varid, UNITS, trange_units), "trange units attribute error", error)
   CALL check (nf90_put_att(ncid, trange_error_varid, UNITS, trange_error_units), "trange_error units attribute error", &
  & error)
 
   CALL check (nf90_put_att(ncid, pcp_varid_2, UNITS, PCP_UNITS), "pcp units attribute error", error)
   CALL check (nf90_put_att(ncid, pop_varid_2, UNITS, POP_UNITS), "pcp units attribute error", error)
   CALL check (nf90_put_att(ncid, pcp_error_varid_2, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)
 
   CALL check (nf90_put_att(ncid, tmean_varid_2, UNITS, tmean_units), "tmean units attribute error", error)
   CALL check (nf90_put_att(ncid, tmean_error_varid_2, UNITS, tmean_error_units), "tmean_error units attribute error", &
  & error)
   CALL check (nf90_put_att(ncid, trange_varid_2, UNITS, trange_units), "trange units attribute error", error)
   CALL check (nf90_put_att(ncid, trange_error_varid_2, UNITS, trange_error_units), "trange_error units attribute error&
  &", error)
 
 
     !correlation variables
   CALL check (nf90_put_att(ncid, autoc_varid, UNITS, autoc_units), "auto correlation units attribute error", error)
   CALL check (nf90_put_att(ncid, tpc_varid, UNITS, tpc_units), "tp correlation units attribute error", error)
 
     !transformed mean,std variables, min, max of normalized y
   CALL check (nf90_put_att(ncid, ymean_varid, UNITS, y_mean_units), "ymean units attribute error", error)
   CALL check (nf90_put_att(ncid, ystd_varid, UNITS, y_std_units), "ystd units attribute error", error)
   CALL check (nf90_put_att(ncid, ystdall_varid, UNITS, y_std_units), "ystd_all units attribute error", error)
   CALL check (nf90_put_att(ncid, ymax_varid, UNITS, y_max_units), "ymax units attribute error", error)
   CALL check (nf90_put_att(ncid, ymin_varid, UNITS, y_min_units), "ymin units attribute error", error)
 
   IF (error /= 0) RETURN
 
 
 
 
     ! End define mode.
   CALL check (nf90_enddef(ncid), "end define mode error", error)
   IF (error /= 0) RETURN
 
   count2 = (/ inx, iny /)
   start2 = (/ 1, 1 /)
 
   CALL check (nf90_put_var(ncid, lat_varid, grdlat, start=start2, count=count2), "put lat error", error)
   CALL check (nf90_put_var(ncid, lon_varid, grdlon, start=start2, count=count2), "put lon error", error)
   CALL check (nf90_put_var(ncid, alt_varid, grdalt, start=start2, count=count2), "put alt error", error)
 
   trec = 1
   nrecs = 0
 
 
  ELSE
 
     ! File already exists, get dim and var ids
   CALL check (nf90_inq_dimid(ncid, X_NAME, x_dimid), "x dim inq error", error)
   CALL check (nf90_inq_dimid(ncid, Y_NAME, y_dimid), "y dim inq error", error)
   CALL check (nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_inq_varid(ncid, LAT_NAME, lat_varid), "lat var inq error", error)
   CALL check (nf90_inq_varid(ncid, LON_NAME, lon_varid), "lon var inq error", error)
   CALL check (nf90_inq_varid(ncid, TIME_NAME, time_varid), "time var inq error", error)
   CALL check (nf90_inq_varid(ncid, PCP_NAME, pcp_varid), "pcp var inq error", error)
   CALL check (nf90_inq_varid(ncid, POP_NAME, pop_varid), "pop var inq error", error)
   CALL check (nf90_inq_varid(ncid, PCP_ERROR_NAME, pcp_error_varid), "pcp_error var inq error", error)
   CALL check (nf90_inq_varid(ncid, tmean_name, tmean_varid), "tmean var inq error", error)
   CALL check (nf90_inq_varid(ncid, tmean_error_name, tmean_error_varid), "tmean error var inq error", error)
   CALL check (nf90_inq_varid(ncid, trange_name, trange_varid), "trange var inq error", error)
   CALL check (nf90_inq_varid(ncid, trange_error_name, trange_error_varid), "trange error var inq error", error)
 
   CALL check (nf90_inq_varid(ncid, PCP_NAME_2, pcp_varid_2), "pcp var inq error", error)
   CALL check (nf90_inq_varid(ncid, POP_NAME_2, pop_varid_2), "pop var inq error", error)
   CALL check (nf90_inq_varid(ncid, PCP_ERROR_NAME_2, pcp_error_varid_2), "pcp_error var inq error", error)
   CALL check (nf90_inq_varid(ncid, tmean_name_2, tmean_varid_2), "tmean var inq error", error)
   CALL check (nf90_inq_varid(ncid, tmean_error_name_2, tmean_error_varid_2), "tmean error var inq error", error)
   CALL check (nf90_inq_varid(ncid, trange_name_2, trange_varid_2), "trange var inq error", error)
   CALL check (nf90_inq_varid(ncid, trange_error_name_2, trange_error_varid_2), "trange error var inq error", error)
 
 
   CALL check (nf90_inq_varid(ncid, autoc_name, autoc_varid), "autoc var inq error", error)
   CALL check (nf90_inq_varid(ncid, tpc_name, tpc_varid), "tpc var inq error", error)
 
   CALL check (nf90_inq_varid(ncid, y_mean_name, ymean_varid), "ymean var inq error", error)
   CALL check (nf90_inq_varid(ncid, y_std_name, ystd_varid), "ystd var inq error", error)
   CALL check (nf90_inq_varid(ncid, y_stdall_name, ystdall_varid), "ystd_all var inq error", error)
   CALL check (nf90_inq_varid(ncid, y_max_name, ymax_varid), "ymax var inq error", error)
   CALL check (nf90_inq_varid(ncid, y_min_name, ymin_varid), "ymin var inq error", error)
   IF (error /= 0) RETURN
 
   CALL check (nf90_Inquire_Dimension(ncid, x_dimid, len=file_nx), "x dim len error", error)
   CALL check (nf90_Inquire_Dimension(ncid, y_dimid, len=file_ny), "y dim len error", error)
   CALL check (nf90_Inquire_Dimension(ncid, time_dimid, len=file_ntimes), "time dim len error", error)
   IF (error /= 0) RETURN
 
   IF (nx /= file_nx .OR. ny /= file_ny) THEN
    PRINT *, "Error dimensions in output file do not match current run."
    error = 1
    RETURN
   END IF
 
   ALLOCATE (file_times(file_ntimes))
   CALL check (nf90_get_var(ncid, time_varid, file_times), "error getting file times list", error)
   IF (error /= 0) RETURN
 
   IF (file_times(1) > times(n_times)) THEN !put data before everything in the file
    PRINT *, "Error cannot add data before data already in output file. (functionality still to be added)"
    error = 1
    RETURN
   ELSE
    IF (file_times(file_ntimes) < times(1)) THEN !put data after everything in the file
     trec = file_ntimes + 1
    ELSE ! at least some overlap
     DO i = 1, file_ntimes, 1
      IF (file_times(1) == times(1)) THEN
       trec = i
      END IF
     END DO
     IF (trec == 0) THEN
      PRINT *, "Error, confusion over data output record location."
      error = 1
      RETURN
     ELSE
      PRINT *, "WARNING, overwriting data in output file, record ", trec, " to ", trec + n_times - 1
     END IF
    END IF
   END IF
 
  END IF
 
  count1 (1) = n_times
  start1 (1) = trec
  CALL check (nf90_put_var(ncid, time_varid, times, start=start1, count=count1), "put times error", error)
  IF (error /= 0) RETURN
 
 
  !correlation variables
  CALL check (nf90_put_var(ncid, autoc_varid, mean_autocorr, start=start1, count=count1), "put mean autocorrelation err&
 &or", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, tpc_varid, mean_tp_corr, start=start1, count=count1), "put mean t_p correlation error",&
 &  error)
  IF (error /= 0) RETURN
 
 
 
 
  !3-d variables
  count3 = (/ inx, iny, n_times /)
  start3 = (/ 1, 1, trec /)
  CALL check (nf90_put_var(ncid, pcp_varid, real(pcp, kind(DP)), start=start3, count=count3), "put pcp error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, pop_varid, real(pop, kind(DP)), start=start3, count=count3), "put pop error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, pcp_error_varid, real(pcperror, kind(DP)), start=start3, count=count3), "put pcp_error&
 & error", error)
  IF (error /= 0) RETURN
 
 
  CALL check (nf90_put_var(ncid, tmean_varid, real(tmean, kind(DP)), start=start3, count=count3), "put tmean error", &
 & error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, tmean_error_varid, real(tmean_error, kind(DP)), start=start3, count=count3), "put tmea&
 &n_error error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, trange_varid, real(trange, kind(DP)), start=start3, count=count3), "put trange error", &
 & error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, trange_error_varid, real(trange_error, kind(DP)), start=start3, count=count3), "put tr&
 &ange_error error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, pcp_varid_2, real(pcp_2, kind(DP)), start=start3, count=count3), "put pcp error", &
 & error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, pop_varid_2, real(pop_2, kind(DP)), start=start3, count=count3), "put pop error", &
 & error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, pcp_error_varid_2, real(pcperror_2, kind(DP)), start=start3, count=count3), "put pcp_e&
 &rror error", error)
  IF (error /= 0) RETURN
 
 
  CALL check (nf90_put_var(ncid, tmean_varid_2, real(tmean_2, kind(DP)), start=start3, count=count3), "put tmean error",&
 &  error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, tmean_error_varid_2, real(tmean_error_2, kind(DP)), start=start3, count=count3), "put &
 &tmean_error error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, trange_varid_2, real(trange_2, kind(DP)), start=start3, count=count3), "put trange err&
 &or", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, trange_error_varid_2, real(trange_error_2, kind(DP)), start=start3, count=count3), "pu&
 &t trange_error error", error)
  IF (error /= 0) RETURN
 
 
 
!transformed mean,std variables, min & max of normalized y
  CALL check (nf90_put_var(ncid, ymean_varid, y_mean, start=start3, count=count3), "put ymean error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, ystd_varid, y_std, start=start3, count=count3), "put ystd error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, ystdall_varid, y_std_all, start=start3, count=count3), "put ystd_all error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, ymax_varid, y_min, start=start3, count=count3), "put ymax error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, ymin_varid, y_max, start=start3, count=count3), "put ymin error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_close(ncid), "closing file error", error)
 
 
CONTAINS
  SUBROUTINE check (status, info, error)
   INTEGER, INTENT (IN) :: status
   CHARACTER (LEN=*), INTENT (IN) :: info
   INTEGER, INTENT (OUT) :: error
 
   IF (status /= nf90_noerr) THEN
    PRINT *, trim (info) // ": " // trim (nf90_strerror(status))
    error = 1
   END IF
  END SUBROUTINE check
END SUBROUTINE save_precip
