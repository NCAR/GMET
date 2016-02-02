
subroutine save_coefficients(n_vars, var_names, coefs, startdate, enddate, times, &
     site_list, station_var, stnid, stnlat, stnlon, stnalt, forecast, file, error)
  use netcdf
  use type
  implicit none

  character (len = 100), intent(in) :: var_names(:), stnid(:)
  integer, intent(in) :: n_vars, forecast
  character (len = 100), intent(in) :: startdate, enddate, station_var
  character (len = 500), intent(in) :: file, site_list
  real(DP), intent(in) :: stnlat(:), stnlon(:), stnalt(:)
  real(DP), intent(in) :: coefs(:,:,:)
  real(DP), intent(in) :: times(:)

  integer, intent(out) :: error

  ! Dimension names
  character (len = *), parameter :: STN_NAME = "station"
  character (len = *), parameter :: VAR_NAME = "variable"
  character (len = *), parameter :: CHAR_NAME = "string"
  character (len = *), parameter :: TIME_NAME = "time"
  character (len = *), parameter :: REC_NAME = "run"
  ! Variable Names
  character (len = *), parameter :: VARS_NAME = "variable_name"
  character (len = *), parameter :: STN_FILE_NAME = "station_file_name"
  character (len = *), parameter :: STN_ID_NAME = "station_id"
  character (len = *), parameter :: STN_LAT_NAME = "station_latitude"
  character (len = *), parameter :: STN_LON_NAME = "station_longitude"
  character (len = *), parameter :: STN_ALT_NAME = "station_altitude"
  character (len = *), parameter :: FORECAST_NAME = "forecast_hr"
  character (len = *), parameter :: STARTDATE_NAME = "run_start_date"
  character (len = *), parameter :: ENDDATE_NAME = "run_end_date"
  character (len = *), parameter :: TIME_VAR_NAME = "time"
  character (len = *), parameter :: COEFS_NAME = "coefficient"
  character (len = *), parameter :: CONSTANT_NAME = "constant_term"
  ! Units
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: COEFS_UNITS = "linear_equation_coefficient"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: ALT_UNITS = "feet"
  character (len = *), parameter :: FORECAST_UNITS = "hours"
  character (len = *), parameter :: DATE_UNITS = "YYYYMMDD"
  character (len = *), parameter :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len = *), parameter :: FILL = "_FillValue"

  integer :: n_stns, n_chars, n_times
  integer :: ncid, stn_dimid, vars_dimid, char_dimid, time_dimid, rec_dimid
  integer :: coefs_varid, time_varid, vars_varid, forecast_varid, startdate_varid, enddate_varid
  integer :: stn_id_varid, stn_lat_varid, stn_lon_varid, stn_alt_varid
  integer :: charids(2)
  integer :: dimids(4)
  integer :: count4(4), start4(4), count3(3), start3(3), count2(2), start2(2), count1(1), start1(1)
  integer :: nstns, nvars, nrecs, ntimes, rec, i

  character (len = 100) :: startdates, enddates
  integer, allocatable :: forecasts(:)

  integer :: forecast_arr(1)

  n_chars = 100
  n_stns = size(stnlat)
  n_times = size(times)

  error = nf90_open(file, nf90_write, ncid)
  if(error /= nf90_noerr) then
     error = 0
     ! Create the file. 
     call check( nf90_create(file, nf90_clobber, ncid), "File creation error", error)
     if(error /=0 ) return

     ! Define the dimensions. 
     call check( nf90_def_dim(ncid, STN_NAME, n_stns, stn_dimid), "station dim def error", error)
     call check( nf90_def_dim(ncid, VAR_NAME, n_vars+1, vars_dimid), "variable dim def error", error)
     call check( nf90_def_dim(ncid, CHAR_NAME, n_chars, char_dimid), "char dim def error", error)
     call check( nf90_def_dim(ncid, TIME_NAME, n_times, time_dimid), "time dim def error", error)
     call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid), "rec dim def error", error)
     if(error /=0 ) return

     ! Define the variables. 
     charids = (/ char_dimid, vars_dimid /)
     call check( nf90_def_var(ncid, VARS_NAME, NF90_CHAR, charids, vars_varid), &
          "variable_name var def error", error)
     charids = (/ char_dimid, stn_dimid /)
     call check( nf90_def_var(ncid, STN_ID_NAME, NF90_CHAR, charids, stn_id_varid), &
          "station_id var def error", error)
     call check( nf90_def_var(ncid, STN_LAT_NAME, NF90_DOUBLE, stn_dimid, stn_lat_varid), &
          "station_latitude var def error", error)
     call check( nf90_def_var(ncid, STN_LON_NAME, NF90_DOUBLE, stn_dimid, stn_lon_varid), &
          "station_longitude var def error", error)
     call check( nf90_def_var(ncid, STN_ALT_NAME, NF90_DOUBLE, stn_dimid, stn_alt_varid), &
          "station_altitude var def error", error)
     if(error /=0 ) return

     charids = (/ char_dimid, rec_dimid /)
     call check( nf90_def_var(ncid, STARTDATE_NAME, NF90_CHAR, charids, startdate_varid), &
          "start_date var def error", error)
     call check( nf90_def_var(ncid, ENDDATE_NAME, NF90_CHAR, charids, enddate_varid), &
          "end_date var def error", error)
     call check( nf90_def_var(ncid, FORECAST_NAME, NF90_INT, rec_dimid, forecast_varid), &
          "forecast var def error", error)
     if(error /=0 ) return

     call check( nf90_def_var(ncid, TIME_VAR_NAME, NF90_DOUBLE, time_dimid, time_varid), &
          "time var def error", error)
     dimids = (/ stn_dimid, time_dimid, vars_dimid, rec_dimid /)
     call check( nf90_def_var(ncid, COEFS_NAME, NF90_DOUBLE, dimids, coefs_varid), &
          "coefficient var def error", error)
     if(error /=0 ) return

     ! Add attributes. 
     call check( nf90_put_att(ncid, stn_lat_varid, UNITS, LAT_UNITS), "station_lat units attribute error", error)
     call check( nf90_put_att(ncid, stn_lon_varid, UNITS, LON_UNITS), "station_lon units attribute error", error)
     call check( nf90_put_att(ncid, stn_alt_varid, UNITS, ALT_UNITS), "station_alt units attribute error", error)
     call check( nf90_put_att(ncid, startdate_varid, UNITS, DATE_UNITS), "start_date units attribute error", error)
     call check( nf90_put_att(ncid, enddate_varid, UNITS, DATE_UNITS), "end_date units attribute error", error)
     call check( nf90_put_att(ncid, forecast_varid, UNITS, FORECAST_UNITS), "forecase units attribute error", error )
     call check( nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)
     call check( nf90_put_att(ncid, coefs_varid, UNITS, COEFS_UNITS), "coefficient units attribute error", error)
     call check( nf90_put_att(ncid, NF90_GLOBAL, STN_FILE_NAME, site_list), &
          "station_file_name global attribute error", error)

     ! End define mode.
     call check( nf90_enddef(ncid), "end define mode error", error)
     if(error /=0 ) return

     call check( nf90_put_var(ncid, stn_lat_varid, stnlat), "put staion_lat error", error) 
     call check( nf90_put_var(ncid, stn_lon_varid, stnlon), "put staion_lon error", error) 
     call check( nf90_put_var(ncid, stn_alt_varid, stnalt), "put staion_alt error", error) 
     call check( nf90_put_var(ncid, time_varid, times), "put times error", error) 

     count2 = (/ 1, 1 /)
     start2 = (/ 1, 1 /)
     do i = 1, n_stns, 1
        count2(1) = len(trim(stnid(i)))
        start2(2) = i
        call check( nf90_put_var(ncid, stn_id_varid, stnid(i), start = start2, count = count2), &
             "put staion_id error", error)
        if(error /=0 ) return
     end do

     count2(1) = len(trim(CONSTANT_NAME))
     start2(2) = 1
     call check( nf90_put_var(ncid, vars_varid, CONSTANT_NAME, start = start2, count = count2), &
          "put variable_name error", error)
     if(error /=0 ) return
     do i = 1, n_vars, 1
        count2(1) = len(trim(var_names(i)))
        start2(2) = i+1
        call check( nf90_put_var(ncid, vars_varid, var_names(i), start = start2, count = count2), &
             "put variable_name error", error)
        if(error /=0 ) return
     end do

     rec = 1
     nrecs = 0

  else

     ! File already exists, get dim and var ids
     call check( nf90_inq_dimid(ncid, STN_NAME, stn_dimid), "station dim inq error", error)
     call check( nf90_inq_dimid(ncid, VAR_NAME, vars_dimid), "variable dim inq error", error)
     call check( nf90_inq_dimid(ncid, CHAR_NAME, char_dimid), "char dim inq error", error)
     call check( nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
     call check( nf90_inq_dimid(ncid, REC_NAME, rec_dimid), "run dim inq error", error)
     if(error /=0 ) return

     call check( nf90_inq_varid(ncid, VARS_NAME, vars_varid), "variable_name var inq error", error)
     call check( nf90_inq_varid(ncid, STN_ID_NAME, stn_id_varid), "station_id var inq error", error)
     call check( nf90_inq_varid(ncid, STN_LAT_NAME, stn_lat_varid), "station_latitude var inq error", error)
     call check( nf90_inq_varid(ncid, STN_LON_NAME, stn_lon_varid), "station_longitude var inq error", error)
     call check( nf90_inq_varid(ncid, STN_ALT_NAME, stn_alt_varid), "station_altitude var inq error", error)
     call check( nf90_inq_varid(ncid, STARTDATE_NAME, startdate_varid), "start_date var inq error", error)
     call check( nf90_inq_varid(ncid, ENDDATE_NAME, enddate_varid), "end_date var inq error", error)
     call check( nf90_inq_varid(ncid, FORECAST_NAME, forecast_varid), "forecast var inq error", error)
     call check( nf90_inq_varid(ncid, TIME_VAR_NAME, time_varid), "time var inq error", error)
     call check( nf90_inq_varid(ncid, COEFS_NAME, coefs_varid), "coefficient var inq error", error)
     if(error /= 0 ) return

     ! Verify Dimensions match
     call check( nf90_Inquire_Dimension(ncid, stn_dimid, len = nstns), "station dim len error", error)
     call check( nf90_inquire_dimension(ncid, vars_dimid, len = nvars), "variable dim len error", error)
     call check( nf90_inquire_dimension(ncid, time_dimid, len = ntimes), "time dim len error", error)
     call check( nf90_inquire_dimension(ncid, rec_dimid, len = nrecs), "run dim len error", error)
     if(error /= 0 ) return

     if (n_stns /= nstns .or. nvars /= n_vars+1 .or. ntimes /= n_times) then
        print *, "Error dimensions in output file do not match current run."
        error = 1
        return
     endif

     allocate(forecasts(nrecs))

     call check( nf90_get_var(ncid, forecast_varid, forecasts), "error getting file forecasts", error)
     if(error /= 0 ) return

     rec = nrecs + 1
     do i = 1, nrecs, 1
        count2 = (/ 100, 1 /)
        start2 = (/ 1, i /)
        call check( nf90_get_var(ncid, startdate_varid, startdates, start = start2, count = count2), &
             "error getting file startdates", error)
        call check( nf90_get_var(ncid, enddate_varid, enddates, start = start2, count = count2), &
             "error getting file enddates", error)
        if (forecasts(i) == forecast .and. startdates(1:8) == startdate(1:8) .and. enddates(1:8) == enddate(1:8)) then
           print *, "WARNING, overwriting data in output file, record ", i
           rec = i
        endif
     end do


  endif

  count2 = (/ len(trim(startdate)), 1 /)
  start2 = (/ 1, rec /)
  call check( nf90_put_var(ncid, startdate_varid, startdate, start = start2, count = count2), &
       "put start_date error", error)
  if(error /=0 ) return

  count2 = (/ len(trim(enddate)), 1 /)
  start2 = (/ 1, rec /)
  call check( nf90_put_var(ncid, enddate_varid, enddate, start = start2, count = count2), &
       "put start_date error", error)
  if(error /=0 ) return

  count1(1) = 1
  start1(1) = rec
  forecast_arr = forecast
  call check( nf90_put_var(ncid, forecast_varid, forecast_arr, start = start1, count = count1), &
       "put forecast error", error)
  if(error /=0 ) return

  count4 = (/ n_stns, n_times, n_vars+1, 1 /)
  start4 = (/ 1, 1, 1, rec /)
  forecast_arr = forecast
  call check( nf90_put_var(ncid, coefs_varid, coefs, start = start4, count = count4), &
       "put coefficients error", error)
  if(error /=0 ) return

  call check( nf90_close(ncid), "closing file error", error)

contains
  subroutine check(status, info, error)
    integer, intent (in) :: status
    character (len = *), intent(in) :: info
    integer, intent(out) :: error

    if(status /= nf90_noerr) then
       print *, trim(info)//": "// trim(nf90_strerror(status))
       error = 1
    end if
  end subroutine check  
end subroutine save_coefficients

subroutine save_precip(pcp, pop, pcperror, tmean, tmean_error, trange, trange_error, &
		       nx, ny, grdlat, grdlon, grdalt, Times, mean_autocorr, mean_tp_corr, &
                       y_mean, y_std,y_std_all,y_min,y_max,file, error, &
		  pcp_2, pop_2, pcperror_2, tmean_2, tmean_error_2, trange_2, trange_error_2)
  use netcdf
  use type
  implicit none

  real(SP), intent(in) :: pcp(:,:), pop(:,:), pcperror(:,:)
  real(SP), intent(in) :: tmean(:,:),tmean_error(:,:),trange(:,:),trange_error(:,:)

  real(SP), intent(in) :: pcp_2(:,:), pop_2(:,:), pcperror_2(:,:)
  real(SP), intent(in) :: tmean_2(:,:),tmean_error_2(:,:),trange_2(:,:),trange_error_2(:,:)

  integer(I4B), intent(in) :: nx, ny
  real(DP), intent(in) :: grdlat(:), grdlon(:), grdalt(:)
  real(DP), intent(in) :: Times(:)
  real(DP), intent(in) :: mean_autocorr(:),mean_tp_corr(:)

  real(DP), intent(in) :: y_mean(:,:), y_std(:,:),y_min(:,:),y_max(:,:),y_std_all(:,:)
  character (len = 500), intent(in) :: file
  integer, intent(out) :: error


  ! Dimension names
  character (len = *), parameter :: Y_NAME = "y"
  character (len = *), parameter :: X_NAME = "x"
  character (len = *), parameter :: TIME_NAME = "time"

  ! Variable Names
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: autoc_name = "auto_corr"
  character (len = *), parameter :: tpc_name = "tp_corr"
  character (len = *), parameter :: ALT_NAME = "altitude"
  character (len = *), parameter :: PCP_NAME = "pcp"
  character (len = *), parameter :: POP_NAME = "pop"
  character (len = *), parameter :: PCP_ERROR_NAME = "pcp_error"
  character (len = *), parameter :: tmean_name = "tmean"
  character (len = *), parameter :: tmean_error_name = "tmean_error"
  character (len = *), parameter :: trange_name = "trange"
  character (len = *), parameter :: trange_error_name = "trange_error"
  character (len = *), parameter :: y_mean_name = "ymean"
  character (len = *), parameter :: y_std_name = "ystd"
  character (len = *), parameter :: y_stdall_name = "ystd_all"
  character (len = *), parameter :: y_min_name = "ymax"
  character (len = *), parameter :: y_max_name = "ymin"
  character (len = *), parameter :: PCP_NAME_2 = "pcp_2"
  character (len = *), parameter :: POP_NAME_2 = "pop_2"
  character (len = *), parameter :: PCP_ERROR_NAME_2 = "pcp_error_2"
  character (len = *), parameter :: tmean_name_2 = "tmean_2"
  character (len = *), parameter :: tmean_error_name_2 = "tmean_error_2"
  character (len = *), parameter :: trange_name_2 = "trange_2"
  character (len = *), parameter :: trange_error_name_2 = "trange_error_2"


  character (len = *), parameter :: LONG_NAME = "long_name"
  character (len = *), parameter :: PCP_LONG_NAME = "estimated precip in normal space"
  character (len = *), parameter :: POP_LONG_NAME = "probability of precipitation occurrence"
  character (len = *), parameter :: PCP_ERROR_LONG_NAME = "error in estimated precip"
  character (len = *), parameter :: tmean_long_name = "estimated daily mean temperature"
  character (len = *), parameter :: tmean_error_long_name = "error in estimated daily mean temp"
  character (len = *), parameter :: trange_long_name = "estimated diurnal range"
  character (len = *), parameter :: trange_error_long_name = "error in estimated diurnal range"
  character (len = *), parameter :: autoc_long_name = "Lag-1 autocorrelation of temperature"
  character (len = *), parameter :: tpc_long_name   = "Correlation of diurnal range and precipitation"
  character (len = *), parameter :: y_mean_long_name = "mean of transformed non-zero precip"
  character (len = *), parameter :: y_std_long_name = "std. dev. of transformed non-zero precip"
  character (len = *), parameter :: y_max_long_name = "max of normalized transformed non-zero precip"
  character (len = *), parameter :: y_min_long_name = "min of normalized transformed non-zero precip"
  character (len = *), parameter :: PCP_LONG_NAME_2 = "estimated precip in normal space (no slope)"
  character (len = *), parameter :: POP_LONG_NAME_2 = "probability of precipitation occurrence (no slope)"
  character (len = *), parameter :: PCP_ERROR_LONG_NAME_2 = "error in estimated precip (no slope)"
  character (len = *), parameter :: tmean_long_name_2 = "estimated daily mean temperature (no slope)"
  character (len = *), parameter :: tmean_error_long_name_2 = "error in estimated daily mean temp (no slope)"
  character (len = *), parameter :: trange_long_name_2 = "estimated diurnal range (no slope)"
  character (len = *), parameter :: trange_error_long_name_2 = "error in estimated diurnal range (no slope)"
 

  ! Units
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PCP_UNITS = ""
  character (len = *), parameter :: POP_UNITS = ""
  character (len = *), parameter :: autoc_units = ""
  character (len = *), parameter :: tpc_units = ""
  character (len = *), parameter :: PCP_ERROR_UNITS = ""
  character (len = *), parameter :: tmean_units = "deg_C"
  character (len = *), parameter :: trange_units = "deg_C"
  character (len = *), parameter :: tmean_error_units = "deg_C"
  character (len = *), parameter :: trange_error_units = "deg_C"
  character (len = *), parameter :: y_mean_units = ""
  character (len = *), parameter :: y_std_units = ""
  character (len = *), parameter :: y_max_units = ""
  character (len = *), parameter :: y_min_units = ""

  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: ALT_UNITS = "meters"
  character (len = *), parameter :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len = *), parameter :: FILL = "_FillValue"

  real(DP), allocatable :: file_times(:)

  integer :: n_chars, n_times, inx, iny
  integer :: ncid, x_dimid, y_dimid, time_dimid
  integer :: lat_varid, lon_varid, autoc_varid, alt_varid, time_varid, pcp_varid, pop_varid, pcp_error_varid,tpc_varid
  integer :: tmean_varid, tmean_error_varid, trange_varid, trange_error_varid
  integer :: ymean_varid, ystd_varid, ymax_varid, ymin_varid,ystdall_varid
  integer :: count1(1), start1(1), count2(2), start2(2), count3(3), start3(3), dimids2(2), dimids3(3)
  integer :: trec, nrecs, file_nx, file_ny, file_ntimes, i

  integer :: pcp_varid_2, pop_varid_2, pcp_error_varid_2
  integer :: tmean_varid_2, tmean_error_varid_2, trange_varid_2, trange_error_varid_2

  trec = 0
  n_chars = 100
  n_times = size(Times)
  inx = nx
  iny = ny

  if(size(grdlat) /= inx * iny) then
     print *, "Error "
  endif

  error = nf90_open(file, nf90_write, ncid)
  if(error /= nf90_noerr) then
     error = 0
     ! Create the file. 
     call check( nf90_create(file, nf90_clobber, ncid), "File creation error", error)
     if(error /=0 ) return

     ! Define the dimensions. 
     call check( nf90_def_dim(ncid, Y_NAME, iny, y_dimid), "y dim def error", error)
     call check( nf90_def_dim(ncid, X_NAME, inx, x_dimid), "x dim def error", error)
     call check( nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), "time dim def error", error)
     if(error /=0 ) return

     ! Define the variables. 
     dimids2 = (/ x_dimid, y_dimid /)
     call check( nf90_def_var(ncid, LAT_NAME, NF90_DOUBLE, dimids2, lat_varid), &
          "lat var def error", error)
     call check( nf90_def_var(ncid, LON_NAME, NF90_DOUBLE, dimids2, lon_varid), &
          "lon var def error", error)
     call check( nf90_def_var(ncid, ALT_NAME, NF90_DOUBLE, dimids2, alt_varid), &
          "lon var def error", error)
     call check( nf90_def_var(ncid, TIME_NAME, NF90_DOUBLE, time_dimid, time_varid), &
          "time var def error", error)

     !correlation variables
     call check( nf90_def_var(ncid, autoc_name, nf90_double,time_dimid,autoc_varid), &
	  "auto correlaion var def error",error)
     call check( nf90_def_var(ncid, tpc_name, nf90_double,time_dimid,tpc_varid), &
	  "tp correlaion var def error",error)



     if(error /=0 ) return

     dimids3 = (/ x_dimid, y_dimid, time_dimid /)
     call check( nf90_def_var(ncid, PCP_NAME, NF90_DOUBLE, dimids3, pcp_varid), &
          "pcp var def error", error)
     call check( nf90_def_var(ncid, POP_NAME, NF90_DOUBLE, dimids3, pop_varid), &
          "pop var def error", error)
     call check( nf90_def_var(ncid, PCP_ERROR_NAME, NF90_DOUBLE, dimids3, pcp_error_varid), &
          "pcp_error var def error", error)
     if(error /=0 ) return

     call check( nf90_def_var(ncid, tmean_name, NF90_DOUBLE, dimids3, tmean_varid), &
          "tmean var def error", error)
     call check( nf90_def_var(ncid, tmean_error_name, NF90_DOUBLE, dimids3, tmean_error_varid), &
          "tmean error var def error", error)
     call check( nf90_def_var(ncid, trange_name, NF90_DOUBLE, dimids3, trange_varid), &
          "trange var def error", error)
     call check( nf90_def_var(ncid, trange_error_name, NF90_DOUBLE, dimids3, trange_error_varid), &
          "trange error var def error", error)
     if(error /=0 ) return

     call check( nf90_def_var(ncid, PCP_NAME_2, NF90_DOUBLE, dimids3, pcp_varid_2), &
          "pcp var def error", error)
     call check( nf90_def_var(ncid, POP_NAME_2, NF90_DOUBLE, dimids3, pop_varid_2), &
          "pop var def error", error)
     call check( nf90_def_var(ncid, PCP_ERROR_NAME_2, NF90_DOUBLE, dimids3, pcp_error_varid_2), &
          "pcp_error var def error", error)
     if(error /=0 ) return

     call check( nf90_def_var(ncid, tmean_name_2, NF90_DOUBLE, dimids3, tmean_varid_2), &
          "tmean var def error", error)
     call check( nf90_def_var(ncid, tmean_error_name_2, NF90_DOUBLE, dimids3, tmean_error_varid_2), &
          "tmean error var def error", error)
     call check( nf90_def_var(ncid, trange_name_2, NF90_DOUBLE, dimids3, trange_varid_2), &
          "trange var def error", error)
     call check( nf90_def_var(ncid, trange_error_name_2, NF90_DOUBLE, dimids3, trange_error_varid_2), &
          "trange error var def error", error)
     if(error /=0 ) return


     !transformed mean, std variables, min, max of normalized y
     call check( nf90_def_var(ncid, y_mean_name, nf90_double,dimids3,ymean_varid), &
	  "y_mean var def error",error)
     call check( nf90_def_var(ncid, y_std_name, nf90_double,dimids3,ystd_varid), &
	  "y_std var def error",error)
     call check( nf90_def_var(ncid, y_max_name, nf90_double,dimids3,ymax_varid), &
	  "y_max var def error",error)
     call check( nf90_def_var(ncid, y_min_name, nf90_double,dimids3,ymin_varid), &
	  "y_min var def error",error)
     call check( nf90_def_var(ncid, y_stdall_name, nf90_double,dimids3,ystdall_varid), &
	  "y_std_all var def error",error)

    ! Add attributes. 

     !long names
     call check( nf90_put_att(ncid, pcp_varid, LONG_NAME, PCP_LONG_NAME), "pcp long_name attribute error", error)
     call check( nf90_put_att(ncid, pop_varid, LONG_NAME, POP_LONG_NAME), "pcp long_name attribute error", error)
     call check( nf90_put_att(ncid, pcp_error_varid, LONG_NAME, PCP_ERROR_LONG_NAME), "pcp_error long_name attribute error", error)

     call check( nf90_put_att(ncid, tmean_varid, LONG_NAME, tmean_LONG_NAME), "tmean long_name attribute error", error)
     call check( nf90_put_att(ncid, tmean_error_varid, LONG_NAME, tmean_error_LONG_NAME), &
                 "tmean long_name attribute error", error)
     call check( nf90_put_att(ncid, trange_varid, LONG_NAME, trange_LONG_NAME), "trange long_name attribute error", error)
     call check( nf90_put_att(ncid, trange_error_varid, LONG_NAME, trange_error_LONG_NAME), &
                 "trange long_name attribute error", error)

     call check( nf90_put_att(ncid, pcp_varid_2, LONG_NAME, PCP_LONG_NAME_2), "pcp long_name attribute error", error)
     call check( nf90_put_att(ncid, pop_varid_2, LONG_NAME, POP_LONG_NAME_2), "pcp long_name attribute error", error)
     call check( nf90_put_att(ncid, pcp_error_varid_2, LONG_NAME, PCP_ERROR_LONG_NAME_2), "pcp_error long_name attribute error", error)

     call check( nf90_put_att(ncid, tmean_varid_2, LONG_NAME, tmean_LONG_NAME_2), "tmean long_name attribute error", error)
     call check( nf90_put_att(ncid, tmean_error_varid_2, LONG_NAME, tmean_error_LONG_NAME_2), &
                 "tmean long_name attribute error", error)
     call check( nf90_put_att(ncid, trange_varid_2, LONG_NAME, trange_LONG_NAME_2), "trange long_name attribute error", error)
     call check( nf90_put_att(ncid, trange_error_varid_2, LONG_NAME, trange_error_LONG_NAME_2), &
                 "trange long_name attribute error", error)


     !correlation variables
     call check( nf90_put_att(ncid, autoc_varid, LONG_NAME, autoc_LONG_NAME), &
                 "auto_corr long_name attribute error", error)
     call check( nf90_put_att(ncid, tpc_varid, LONG_NAME, tpc_long_NAME), &
                 "tp_corr long_name attribute error", error)

     !transformed mean, std variables, min, max of normalized y
     call check( nf90_put_att(ncid, ymean_varid, LONG_NAME, y_mean_long_name), &
                 "ymean long_name attribute error", error)
     call check( nf90_put_att(ncid, ystd_varid, LONG_NAME, y_std_long_name), &
                 "ystd long_name attribute error", error)
     call check( nf90_put_att(ncid, ystdall_varid, LONG_NAME, y_std_long_name), &
                 "ystd_all long_name attribute error", error)
     call check( nf90_put_att(ncid, ymax_varid, LONG_NAME, y_max_long_name), &
                 "ymax long_name attribute error", error)
     call check( nf90_put_att(ncid, ymin_varid, LONG_NAME, y_min_long_name), &
                 "ymin long_name attribute error", error)
     !units
     call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS), "lat units attribute error", error)
     call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS), "lon units attribute error", error)
     call check( nf90_put_att(ncid, alt_varid, UNITS, ALT_UNITS), "alt units attribute error", error)
     call check( nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)

   
     call check( nf90_put_att(ncid, pcp_varid, UNITS, PCP_UNITS), "pcp units attribute error", error)
     call check( nf90_put_att(ncid, pop_varid, UNITS, POP_UNITS), "pcp units attribute error", error)
     call check( nf90_put_att(ncid, pcp_error_varid, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)

     call check( nf90_put_att(ncid, tmean_varid, UNITS, tmean_UNITS), "tmean units attribute error", error)
     call check( nf90_put_att(ncid, tmean_error_varid, UNITS, tmean_error_UNITS), &
                 "tmean_error units attribute error", error)
     call check( nf90_put_att(ncid, trange_varid, UNITS, trange_UNITS), "trange units attribute error", error)
     call check( nf90_put_att(ncid, trange_error_varid, UNITS, trange_error_UNITS), &
                 "trange_error units attribute error", error)

     call check( nf90_put_att(ncid, pcp_varid_2, UNITS, PCP_UNITS), "pcp units attribute error", error)
     call check( nf90_put_att(ncid, pop_varid_2, UNITS, POP_UNITS), "pcp units attribute error", error)
     call check( nf90_put_att(ncid, pcp_error_varid_2, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)

     call check( nf90_put_att(ncid, tmean_varid_2, UNITS, tmean_UNITS), "tmean units attribute error", error)
     call check( nf90_put_att(ncid, tmean_error_varid_2, UNITS, tmean_error_UNITS), &
                 "tmean_error units attribute error", error)
     call check( nf90_put_att(ncid, trange_varid_2, UNITS, trange_UNITS), "trange units attribute error", error)
     call check( nf90_put_att(ncid, trange_error_varid_2, UNITS, trange_error_UNITS), &
                 "trange_error units attribute error", error)


     !correlation variables
     call check( nf90_put_att(ncid, autoc_varid,units,autoc_units), "auto correlation units attribute error",error)
     call check( nf90_put_att(ncid, tpc_varid,units, tpc_units), "tp correlation units attribute error",error)

     !transformed mean,std variables, min, max of normalized y
     call check( nf90_put_att(ncid, ymean_varid,units,y_mean_units), "ymean units attribute error",error)
     call check( nf90_put_att(ncid, ystd_varid,units, y_std_units), "ystd units attribute error",error)
     call check( nf90_put_att(ncid, ystdall_varid,units, y_std_units), "ystd_all units attribute error",error)
     call check( nf90_put_att(ncid, ymax_varid,units,y_max_units), "ymax units attribute error",error)
     call check( nf90_put_att(ncid, ymin_varid,units, y_min_units), "ymin units attribute error",error)

     if(error /=0 ) return




     ! End define mode.
     call check( nf90_enddef(ncid), "end define mode error", error)
     if(error /=0 ) return

     count2 = (/ inx, iny /)
     start2 = (/ 1, 1 /)

     call check( nf90_put_var(ncid, lat_varid, grdlat, start = start2, count = count2), "put lat error", error) 
     call check( nf90_put_var(ncid, lon_varid, grdlon, start = start2, count = count2), "put lon error", error) 
     call check( nf90_put_var(ncid, alt_varid, grdalt, start = start2, count = count2), "put alt error", error)

     trec = 1
     nrecs = 0


  else

     ! File already exists, get dim and var ids
     call check( nf90_inq_dimid(ncid, X_NAME, x_dimid), "x dim inq error", error)
     call check( nf90_inq_dimid(ncid, Y_NAME, y_dimid), "y dim inq error", error)
     call check( nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
     if(error /=0 ) return

     call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid), "lat var inq error", error)
     call check( nf90_inq_varid(ncid, LON_NAME, lon_varid), "lon var inq error", error)
     call check( nf90_inq_varid(ncid, TIME_NAME, time_varid), "time var inq error", error)
     call check( nf90_inq_varid(ncid, PCP_NAME, pcp_varid), "pcp var inq error", error)
     call check( nf90_inq_varid(ncid, POP_NAME, pop_varid), "pop var inq error", error)
     call check( nf90_inq_varid(ncid, PCP_ERROR_NAME, pcp_error_varid), "pcp_error var inq error", error)
     call check( nf90_inq_varid(ncid, tmean_NAME, tmean_varid), "tmean var inq error", error)
     call check( nf90_inq_varid(ncid, tmean_error_NAME, tmean_error_varid), "tmean error var inq error", error)
     call check( nf90_inq_varid(ncid, trange_NAME, trange_varid), "trange var inq error", error)
     call check( nf90_inq_varid(ncid, trange_error_NAME, trange_error_varid), "trange error var inq error", error)

     call check( nf90_inq_varid(ncid, PCP_NAME_2, pcp_varid_2), "pcp var inq error", error)
     call check( nf90_inq_varid(ncid, POP_NAME_2, pop_varid_2), "pop var inq error", error)
     call check( nf90_inq_varid(ncid, PCP_ERROR_NAME_2, pcp_error_varid_2), "pcp_error var inq error", error)
     call check( nf90_inq_varid(ncid, tmean_NAME_2, tmean_varid_2), "tmean var inq error", error)
     call check( nf90_inq_varid(ncid, tmean_error_NAME_2, tmean_error_varid_2), "tmean error var inq error", error)
     call check( nf90_inq_varid(ncid, trange_NAME_2, trange_varid_2), "trange var inq error", error)
     call check( nf90_inq_varid(ncid, trange_error_NAME_2, trange_error_varid_2), "trange error var inq error", error)


     call check( nf90_inq_varid(ncid, autoc_NAME, autoc_varid), "autoc var inq error", error)
     call check( nf90_inq_varid(ncid, tpc_NAME, tpc_varid), "tpc var inq error", error)

     call check( nf90_inq_varid(ncid, y_mean_NAME, ymean_varid), "ymean var inq error", error)
     call check( nf90_inq_varid(ncid, y_std_name, ystd_varid), "ystd var inq error", error)
     call check( nf90_inq_varid(ncid, y_stdall_name, ystdall_varid), "ystd_all var inq error", error)
     call check( nf90_inq_varid(ncid, y_max_NAME, ymax_varid), "ymax var inq error", error)
     call check( nf90_inq_varid(ncid, y_min_name, ymin_varid), "ymin var inq error", error)
     if(error /= 0 ) return

     call check( nf90_Inquire_Dimension(ncid, x_dimid, len = file_nx), "x dim len error", error)
     call check( nf90_Inquire_Dimension(ncid, y_dimid, len = file_ny), "y dim len error", error)
     call check( nf90_Inquire_Dimension(ncid, time_dimid, len = file_ntimes), "time dim len error", error)
     if(error /= 0 ) return

     if (nx /= file_nx .or. ny /= file_ny) then
        print *, "Error dimensions in output file do not match current run."
        error = 1
        return
     endif

     allocate(file_times(file_ntimes))
     call check( nf90_get_var(ncid, time_varid, file_times), "error getting file times list", error)
     if(error /= 0 ) return

     if(file_times(1) > Times(n_times)) then  !put data before everything in the file
        print *, "Error cannot add data before data already in output file. (functionality still to be added)"
        error = 1
        return
     else
        if(file_times(file_ntimes) < Times(1)) then !put data after everything in the file
           trec = file_ntimes+1
        else  ! at least some overlap
           do i = 1, file_ntimes, 1
              if(file_times(1) == Times(1)) then
                 trec = i
              endif
           end do
           if(trec == 0) then
              print *, "Error, confusion over data output record location."
              error = 1
              return
           else
              print *, "WARNING, overwriting data in output file, record ", trec, " to ", trec + n_times -1
           endif
        endif
     endif

  endif

  count1(1) = n_times
  start1(1) = trec
  call check( nf90_put_var(ncid, time_varid, times, start = start1, count = count1), &
       "put times error", error)
  if(error /=0 ) return


  !correlation variables
  call check( nf90_put_var(ncid, autoc_varid, mean_autocorr, start = start1, count = count1), &
       "put mean autocorrelation error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, tpc_varid, mean_tp_corr, start = start1, count = count1), &
       "put mean t_p correlation error", error)
  if(error /=0 ) return

  


  !3-d variables
  count3 = (/ inx, iny, n_times /)
  start3 = (/ 1, 1, trec /)
  call check( nf90_put_var(ncid, pcp_varid, real(pcp,kind(dp)), start = start3, count = count3), &
       "put pcp error", error)
  if(error /=0 ) return
  
  call check( nf90_put_var(ncid, pop_varid, real(pop,kind(dp)), start = start3, count = count3), &
       "put pop error", error)
  if(error /=0 ) return
  
  call check( nf90_put_var(ncid, pcp_error_varid, real(pcperror,kind(dp)), start = start3, count = count3), &
       "put pcp_error error", error)
  if(error /=0 ) return
  

  call check( nf90_put_var(ncid, tmean_varid, real(tmean,kind(dp)), start = start3, count = count3), &
       "put tmean error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, tmean_error_varid, real(tmean_error,kind(dp)), start = start3, count = count3), &
       "put tmean_error error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, trange_varid, real(trange,kind(dp)), start = start3, count = count3), &
       "put trange error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, trange_error_varid, real(trange_error,kind(dp)), start = start3, count = count3), &
       "put trange_error error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, pcp_varid_2, real(pcp_2,kind(dp)), start = start3, count = count3), &
       "put pcp error", error)
  if(error /=0 ) return
  
  call check( nf90_put_var(ncid, pop_varid_2, real(pop_2,kind(dp)), start = start3, count = count3), &
       "put pop error", error)
  if(error /=0 ) return
  
  call check( nf90_put_var(ncid, pcp_error_varid_2, real(pcperror_2,kind(dp)), start = start3, count = count3), &
       "put pcp_error error", error)
  if(error /=0 ) return
  

  call check( nf90_put_var(ncid, tmean_varid_2, real(tmean_2,kind(dp)), start = start3, count = count3), &
       "put tmean error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, tmean_error_varid_2, real(tmean_error_2,kind(dp)), start = start3, count = count3), &
       "put tmean_error error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, trange_varid_2, real(trange_2,kind(dp)), start = start3, count = count3), &
       "put trange error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, trange_error_varid_2, real(trange_error_2,kind(dp)), start = start3, count = count3), &
       "put trange_error error", error)
  if(error /=0 ) return



!transformed mean,std variables, min & max of normalized y
  call check( nf90_put_var(ncid, ymean_varid, y_mean, start = start3, count = count3), &
       "put ymean error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, ystd_varid, y_std, start = start3, count = count3), &
       "put ystd error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, ystdall_varid, y_std_all, start = start3, count = count3), &
       "put ystd_all error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, ymax_varid, y_min, start = start3, count = count3), &
       "put ymax error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, ymin_varid, y_max, start = start3, count = count3), &
       "put ymin error", error)
  if(error /=0 ) return

  call check( nf90_close(ncid), "closing file error", error)
  

contains
  subroutine check(status, info, error)
    integer, intent (in) :: status
    character (len = *), intent(in) :: info
    integer, intent(out) :: error

    if(status /= nf90_noerr) then
       print *, trim(info)//": "// trim(nf90_strerror(status))
       error = 1
    end if
  end subroutine check  
end subroutine save_precip
