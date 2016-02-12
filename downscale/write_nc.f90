subroutine save_coefficients (n_vars, var_names, coefs, startdate, enddate, times, site_list, &
& station_var, stnid, stnlat, stnlon, stnalt, forecast, file, error)
  use netcdf
  use type
  implicit none
 
  character (len=100), intent (in) :: var_names (:), stnid (:)
  integer, intent (in) :: n_vars, forecast
  character (len=100), intent (in) :: startdate, enddate, station_var
  character (len=500), intent (in) :: file, site_list
  real (dp), intent (in) :: stnlat (:), stnlon (:), stnalt (:)
  real (dp), intent (in) :: coefs (:, :, :)
  real (dp), intent (in) :: times (:)
 
  integer, intent (out) :: error
 
  ! Dimension names
  character (len=*), parameter :: stn_name = "station"
  character (len=*), parameter :: var_name = "variable"
  character (len=*), parameter :: char_name = "string"
  character (len=*), parameter :: time_name = "time"
  character (len=*), parameter :: rec_name = "run"
  ! Variable Names
  character (len=*), parameter :: vars_name = "variable_name"
  character (len=*), parameter :: stn_file_name = "station_file_name"
  character (len=*), parameter :: stn_id_name = "station_id"
  character (len=*), parameter :: stn_lat_name = "station_latitude"
  character (len=*), parameter :: stn_lon_name = "station_longitude"
  character (len=*), parameter :: stn_alt_name = "station_altitude"
  character (len=*), parameter :: forecast_name = "forecast_hr"
  character (len=*), parameter :: startdate_name = "run_start_date"
  character (len=*), parameter :: enddate_name = "run_end_date"
  character (len=*), parameter :: time_var_name = "time"
  character (len=*), parameter :: coefs_name = "coefficient"
  character (len=*), parameter :: constant_name = "constant_term"
  ! Units
  character (len=*), parameter :: units = "units"
  character (len=*), parameter :: coefs_units = "linear_equation_coefficient"
  character (len=*), parameter :: lat_units = "degrees_north"
  character (len=*), parameter :: lon_units = "degrees_east"
  character (len=*), parameter :: alt_units = "feet"
  character (len=*), parameter :: forecast_units = "hours"
  character (len=*), parameter :: date_units = "YYYYMMDD"
  character (len=*), parameter :: time_units = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len=*), parameter :: fill = "_FillValue"
 
  integer :: n_stns, n_chars, n_times
  integer :: ncid, stn_dimid, vars_dimid, char_dimid, time_dimid, rec_dimid
  integer :: coefs_varid, time_varid, vars_varid, forecast_varid, startdate_varid, enddate_varid
  integer :: stn_id_varid, stn_lat_varid, stn_lon_varid, stn_alt_varid
  integer :: charids (2)
  integer :: dimids (4)
  integer :: count4 (4), start4 (4), count3 (3), start3 (3), count2 (2), start2 (2), count1 (1), &
 & start1 (1)
  integer :: nstns, nvars, nrecs, ntimes, rec, i
 
  character (len=100) :: startdates, enddates
  integer, allocatable :: forecasts (:)
 
  integer :: forecast_arr (1)
 
  n_chars = 100
  n_stns = size (stnlat)
  n_times = size (times)
 
  error = nf90_open (file, nf90_write, ncid)
  if (error /= nf90_noerr) then
    error = 0
     ! Create the file.
    call check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
    if (error /= 0) return
 
     ! Define the dimensions.
    call check (nf90_def_dim(ncid, stn_name, n_stns, stn_dimid), "station dim def error", error)
    call check (nf90_def_dim(ncid, var_name, n_vars+1, vars_dimid), "variable dim def error", &
   & error)
    call check (nf90_def_dim(ncid, char_name, n_chars, char_dimid), "char dim def error", error)
    call check (nf90_def_dim(ncid, time_name, n_times, time_dimid), "time dim def error", error)
    call check (nf90_def_dim(ncid, rec_name, nf90_unlimited, rec_dimid), "rec dim def error", &
   & error)
    if (error /= 0) return
 
     ! Define the variables.
    charids = (/ char_dimid, vars_dimid /)
    call check (nf90_def_var(ncid, vars_name, nf90_char, charids, vars_varid), "variable_name var d&
   &ef error", error)
    charids = (/ char_dimid, stn_dimid /)
    call check (nf90_def_var(ncid, stn_id_name, nf90_char, charids, stn_id_varid), "station_id var &
   &def error", error)
    call check (nf90_def_var(ncid, stn_lat_name, nf90_double, stn_dimid, stn_lat_varid), "station_l&
   &atitude var def error", error)
    call check (nf90_def_var(ncid, stn_lon_name, nf90_double, stn_dimid, stn_lon_varid), "station_l&
   &ongitude var def error", error)
    call check (nf90_def_var(ncid, stn_alt_name, nf90_double, stn_dimid, stn_alt_varid), "station_a&
   &ltitude var def error", error)
    if (error /= 0) return
 
    charids = (/ char_dimid, rec_dimid /)
    call check (nf90_def_var(ncid, startdate_name, nf90_char, charids, startdate_varid), "start_dat&
   &e var def error", error)
    call check (nf90_def_var(ncid, enddate_name, nf90_char, charids, enddate_varid), "end_date var &
   &def error", error)
    call check (nf90_def_var(ncid, forecast_name, nf90_int, rec_dimid, forecast_varid), "forecast v&
   &ar def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, time_var_name, nf90_double, time_dimid, time_varid), "time var d&
   &ef error", error)
    dimids = (/ stn_dimid, time_dimid, vars_dimid, rec_dimid /)
    call check (nf90_def_var(ncid, coefs_name, nf90_double, dimids, coefs_varid), "coefficient var &
   &def error", error)
    if (error /= 0) return
 
     ! Add attributes.
    call check (nf90_put_att(ncid, stn_lat_varid, units, lat_units), "station_lat units attribute e&
   &rror", error)
    call check (nf90_put_att(ncid, stn_lon_varid, units, lon_units), "station_lon units attribute e&
   &rror", error)
    call check (nf90_put_att(ncid, stn_alt_varid, units, alt_units), "station_alt units attribute e&
   &rror", error)
    call check (nf90_put_att(ncid, startdate_varid, units, date_units), "start_date units attribute&
   & error", error)
    call check (nf90_put_att(ncid, enddate_varid, units, date_units), "end_date units attribute err&
   &or", error)
    call check (nf90_put_att(ncid, forecast_varid, units, forecast_units), "forecase units attribut&
   &e error", error)
    call check (nf90_put_att(ncid, time_varid, units, time_units), "time units attribute error", &
   & error)
    call check (nf90_put_att(ncid, coefs_varid, units, coefs_units), "coefficient units attribute e&
   &rror", error)
    call check (nf90_put_att(ncid, nf90_global, stn_file_name, site_list), "station_file_name globa&
   &l attribute error", error)
 
     ! End define mode.
    call check (nf90_enddef(ncid), "end define mode error", error)
    if (error /= 0) return
 
    call check (nf90_put_var(ncid, stn_lat_varid, stnlat), "put staion_lat error", error)
    call check (nf90_put_var(ncid, stn_lon_varid, stnlon), "put staion_lon error", error)
    call check (nf90_put_var(ncid, stn_alt_varid, stnalt), "put staion_alt error", error)
    call check (nf90_put_var(ncid, time_varid, times), "put times error", error)
 
    count2 = (/ 1, 1 /)
    start2 = (/ 1, 1 /)
    do i = 1, n_stns, 1
      count2 (1) = len (trim(stnid(i)))
      start2 (2) = i
      call check (nf90_put_var(ncid, stn_id_varid, stnid(i), start=start2, count=count2), "put stai&
     &on_id error", error)
      if (error /= 0) return
    end do
 
    count2 (1) = len (trim(constant_name))
    start2 (2) = 1
    call check (nf90_put_var(ncid, vars_varid, constant_name, start=start2, count=count2), "put var&
   &iable_name error", error)
    if (error /= 0) return
    do i = 1, n_vars, 1
      count2 (1) = len (trim(var_names(i)))
      start2 (2) = i + 1
      call check (nf90_put_var(ncid, vars_varid, var_names(i), start=start2, count=count2), "put va&
     &riable_name error", error)
      if (error /= 0) return
    end do
 
    rec = 1
    nrecs = 0
 
  else
 
     ! File already exists, get dim and var ids
    call check (nf90_inq_dimid(ncid, stn_name, stn_dimid), "station dim inq error", error)
    call check (nf90_inq_dimid(ncid, var_name, vars_dimid), "variable dim inq error", error)
    call check (nf90_inq_dimid(ncid, char_name, char_dimid), "char dim inq error", error)
    call check (nf90_inq_dimid(ncid, time_name, time_dimid), "time dim inq error", error)
    call check (nf90_inq_dimid(ncid, rec_name, rec_dimid), "run dim inq error", error)
    if (error /= 0) return
 
    call check (nf90_inq_varid(ncid, vars_name, vars_varid), "variable_name var inq error", error)
    call check (nf90_inq_varid(ncid, stn_id_name, stn_id_varid), "station_id var inq error", error)
    call check (nf90_inq_varid(ncid, stn_lat_name, stn_lat_varid), "station_latitude var inq error",&
   &  error)
    call check (nf90_inq_varid(ncid, stn_lon_name, stn_lon_varid), "station_longitude var inq error&
   &", error)
    call check (nf90_inq_varid(ncid, stn_alt_name, stn_alt_varid), "station_altitude var inq error",&
   &  error)
    call check (nf90_inq_varid(ncid, startdate_name, startdate_varid), "start_date var inq error", &
   & error)
    call check (nf90_inq_varid(ncid, enddate_name, enddate_varid), "end_date var inq error", error)
    call check (nf90_inq_varid(ncid, forecast_name, forecast_varid), "forecast var inq error", &
   & error)
    call check (nf90_inq_varid(ncid, time_var_name, time_varid), "time var inq error", error)
    call check (nf90_inq_varid(ncid, coefs_name, coefs_varid), "coefficient var inq error", error)
    if (error /= 0) return
 
     ! Verify Dimensions match
    call check (nf90_inquire_dimension(ncid, stn_dimid, len=nstns), "station dim len error", error)
    call check (nf90_inquire_dimension(ncid, vars_dimid, len=nvars), "variable dim len error", &
   & error)
    call check (nf90_inquire_dimension(ncid, time_dimid, len=ntimes), "time dim len error", error)
    call check (nf90_inquire_dimension(ncid, rec_dimid, len=nrecs), "run dim len error", error)
    if (error /= 0) return
 
    if (n_stns /= nstns .or. nvars /= n_vars+1 .or. ntimes /= n_times) then
      print *, "Error dimensions in output file do not match current run."
      error = 1
      return
    end if
 
    allocate (forecasts(nrecs))
 
    call check (nf90_get_var(ncid, forecast_varid, forecasts), "error getting file forecasts", &
   & error)
    if (error /= 0) return
 
    rec = nrecs + 1
    do i = 1, nrecs, 1
      count2 = (/ 100, 1 /)
      start2 = (/ 1, i /)
      call check (nf90_get_var(ncid, startdate_varid, startdates, start=start2, count=count2), "err&
     &or getting file startdates", error)
      call check (nf90_get_var(ncid, enddate_varid, enddates, start=start2, count=count2), "error g&
     &etting file enddates", error)
      if (forecasts(i) == forecast .and. startdates(1:8) == startdate(1:8) .and. enddates(1:8) == &
     & enddate(1:8)) then
        print *, "WARNING, overwriting data in output file, record ", i
        rec = i
      end if
    end do
 
 
  end if
 
  count2 = (/ len (trim(startdate)), 1 /)
  start2 = (/ 1, rec /)
  call check (nf90_put_var(ncid, startdate_varid, startdate, start=start2, count=count2), "put star&
 &t_date error", error)
  if (error /= 0) return
 
  count2 = (/ len (trim(enddate)), 1 /)
  start2 = (/ 1, rec /)
  call check (nf90_put_var(ncid, enddate_varid, enddate, start=start2, count=count2), "put start_da&
 &te error", error)
  if (error /= 0) return
 
  count1 (1) = 1
  start1 (1) = rec
  forecast_arr = forecast
  call check (nf90_put_var(ncid, forecast_varid, forecast_arr, start=start1, count=count1), "put fo&
 &recast error", error)
  if (error /= 0) return
 
  count4 = (/ n_stns, n_times, n_vars + 1, 1 /)
  start4 = (/ 1, 1, 1, rec /)
  forecast_arr = forecast
  call check (nf90_put_var(ncid, coefs_varid, coefs, start=start4, count=count4), "put coefficients&
 & error", error)
  if (error /= 0) return
 
  call check (nf90_close(ncid), "closing file error", error)
 
contains
  subroutine check (status, info, error)
    integer, intent (in) :: status
    character (len=*), intent (in) :: info
    integer, intent (out) :: error
 
    if (status /= nf90_noerr) then
      print *, trim (info) // ": " // trim (nf90_strerror(status))
      error = 1
    end if
  end subroutine check
end subroutine save_coefficients
 
!modified AJN Sept 2013
!subroutine save_precip(pcp, pop, pcperror, tmean, tmean_error, trange, trange_error, &
! AWW-Feb2016 renamed
subroutine save_forcing_regression (pcp, pop, pcperror, tmean, tmean_error, trange, trange_error, &
& nx, ny, grdlat, grdlon, grdalt, times, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, &
& y_min, y_max, file, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_error_2, trange_2, &
& trange_error_2)
  use netcdf
  use type
  implicit none
 
  real (sp), intent (in) :: pcp (:, :), pop (:, :), pcperror (:, :)
  real (sp), intent (in) :: tmean (:, :), tmean_error (:, :), trange (:, :), trange_error (:, :)
 
  real (sp), intent (in) :: pcp_2 (:, :), pop_2 (:, :), pcperror_2 (:, :)
  real (sp), intent (in) :: tmean_2 (:, :), tmean_error_2 (:, :), trange_2 (:, :), trange_error_2 &
 & (:, :)
 
  integer (i4b), intent (in) :: nx, ny
  real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
  real (dp), intent (in) :: times (:)
  real (dp), intent (in) :: mean_autocorr (:), mean_tp_corr (:)
 
  real (dp), intent (in) :: y_mean (:, :), y_std (:, :), y_min (:, :), y_max (:, :), y_std_all (:, &
 & :)
  character (len=500), intent (in) :: file
  integer, intent (out) :: error
 
  ! Dimension names
  character (len=*), parameter :: y_name = "y"
  character (len=*), parameter :: x_name = "x"
  character (len=*), parameter :: time_name = "time"
 
  ! Variable Names
  character (len=*), parameter :: lat_name = "latitude"
  character (len=*), parameter :: lon_name = "longitude"
  character (len=*), parameter :: autoc_name = "auto_corr"
  character (len=*), parameter :: tpc_name = "tp_corr"
  character (len=*), parameter :: alt_name = "altitude"
  character (len=*), parameter :: pcp_name = "pcp"
  character (len=*), parameter :: pop_name = "pop"
  character (len=*), parameter :: pcp_error_name = "pcp_error"
  character (len=*), parameter :: tmean_name = "tmean"
  character (len=*), parameter :: tmean_error_name = "tmean_error"
  character (len=*), parameter :: trange_name = "trange"
  character (len=*), parameter :: trange_error_name = "trange_error"
  character (len=*), parameter :: y_mean_name = "ymean"
  character (len=*), parameter :: y_std_name = "ystd"
  character (len=*), parameter :: y_stdall_name = "ystd_all"
  character (len=*), parameter :: y_min_name = "ymax"
  character (len=*), parameter :: y_max_name = "ymin"
  character (len=*), parameter :: pcp_name_2 = "pcp_2"
  character (len=*), parameter :: pop_name_2 = "pop_2"
  character (len=*), parameter :: pcp_error_name_2 = "pcp_error_2"
  character (len=*), parameter :: tmean_name_2 = "tmean_2"
  character (len=*), parameter :: tmean_error_name_2 = "tmean_error_2"
  character (len=*), parameter :: trange_name_2 = "trange_2"
  character (len=*), parameter :: trange_error_name_2 = "trange_error_2"
 
 
  character (len=*), parameter :: long_name = "long_name"
  character (len=*), parameter :: pcp_long_name = "estimated precip in normal space"
  character (len=*), parameter :: pop_long_name = "probability of precipitation occurrence"
  character (len=*), parameter :: pcp_error_long_name = "error in estimated precip"
  character (len=*), parameter :: tmean_long_name = "estimated daily mean temperature"
  character (len=*), parameter :: tmean_error_long_name = "error in estimated daily mean temp"
  character (len=*), parameter :: trange_long_name = "estimated diurnal range"
  character (len=*), parameter :: trange_error_long_name = "error in estimated diurnal range"
  character (len=*), parameter :: autoc_long_name = "Lag-1 autocorrelation of temperature"
  character (len=*), parameter :: tpc_long_name = "Correlation of diurnal range and precipitation"
  character (len=*), parameter :: y_mean_long_name = "mean of transformed non-zero precip"
  character (len=*), parameter :: y_std_long_name = "std. dev. of transformed non-zero precip"
  character (len=*), parameter :: y_max_long_name = "max of normalized transformed non-zero precip"
  character (len=*), parameter :: y_min_long_name = "min of normalized transformed non-zero precip"
  character (len=*), parameter :: pcp_long_name_2 = "estimated precip in normal space (no slope)"
  character (len=*), parameter :: pop_long_name_2 = "probability of precipitation occurrence (no sl&
 &ope)"
  character (len=*), parameter :: pcp_error_long_name_2 = "error in estimated precip (no slope)"
  character (len=*), parameter :: tmean_long_name_2 = "estimated daily mean temperature (no slope)"
  character (len=*), parameter :: tmean_error_long_name_2 = "error in estimated daily mean temp (no&
 & slope)"
  character (len=*), parameter :: trange_long_name_2 = "estimated diurnal range (no slope)"
  character (len=*), parameter :: trange_error_long_name_2 = "error in estimated diurnal range (no &
 &slope)"
 
 
  ! Units
  character (len=*), parameter :: units = "units"
  character (len=*), parameter :: pcp_units = ""
  character (len=*), parameter :: pop_units = ""
  character (len=*), parameter :: autoc_units = ""
  character (len=*), parameter :: tpc_units = ""
  character (len=*), parameter :: pcp_error_units = ""
  character (len=*), parameter :: tmean_units = "deg_C"
  character (len=*), parameter :: trange_units = "deg_C"
  character (len=*), parameter :: tmean_error_units = "deg_C"
  character (len=*), parameter :: trange_error_units = "deg_C"
  character (len=*), parameter :: y_mean_units = ""
  character (len=*), parameter :: y_std_units = ""
  character (len=*), parameter :: y_max_units = ""
  character (len=*), parameter :: y_min_units = ""
 
  character (len=*), parameter :: lat_units = "degrees_north"
  character (len=*), parameter :: lon_units = "degrees_east"
  character (len=*), parameter :: alt_units = "meters"
  character (len=*), parameter :: time_units = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len=*), parameter :: fill = "_FillValue"
 
  real (dp), allocatable :: file_times (:)
 
  integer :: n_chars, n_times, inx, iny
  integer :: ncid, x_dimid, y_dimid, time_dimid
  integer :: lat_varid, lon_varid, autoc_varid, alt_varid, time_varid, pcp_varid, pop_varid, &
 & pcp_error_varid, tpc_varid
  integer :: tmean_varid, tmean_error_varid, trange_varid, trange_error_varid
  integer :: ymean_varid, ystd_varid, ymax_varid, ymin_varid, ystdall_varid
  integer :: count1 (1), start1 (1), count2 (2), start2 (2), count3 (3), start3 (3), dimids2 (2), &
 & dimids3 (3)
  integer :: trec, nrecs, file_nx, file_ny, file_ntimes, i
 
  integer :: pcp_varid_2, pop_varid_2, pcp_error_varid_2
  integer :: tmean_varid_2, tmean_error_varid_2, trange_varid_2, trange_error_varid_2
 
  trec = 0
  n_chars = 100
  n_times = size (times)
  inx = nx
  iny = ny
 
  if (size(grdlat) /= inx*iny) then
    print *, "Error "
  end if
 
  error = nf90_open (file, nf90_write, ncid)
  if (error /= nf90_noerr) then
    error = 0
     ! Create the file.
    call check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
    if (error /= 0) return
 
     ! Define the dimensions.
    call check (nf90_def_dim(ncid, y_name, iny, y_dimid), "y dim def error", error)
    call check (nf90_def_dim(ncid, x_name, inx, x_dimid), "x dim def error", error)
    call check (nf90_def_dim(ncid, time_name, nf90_unlimited, time_dimid), "time dim def error", &
   & error)
    if (error /= 0) return
 
     ! Define the variables.
    dimids2 = (/ x_dimid, y_dimid /)
    call check (nf90_def_var(ncid, lat_name, nf90_double, dimids2, lat_varid), "lat var def error", &
   & error)
    call check (nf90_def_var(ncid, lon_name, nf90_double, dimids2, lon_varid), "lon var def error", &
   & error)
    call check (nf90_def_var(ncid, alt_name, nf90_double, dimids2, alt_varid), "lon var def error", &
   & error)
    call check (nf90_def_var(ncid, time_name, nf90_double, time_dimid, time_varid), "time var def e&
   &rror", error)
 
     !correlation variables
    call check (nf90_def_var(ncid, autoc_name, nf90_double, time_dimid, autoc_varid), "auto correla&
   &ion var def error", error)
    call check (nf90_def_var(ncid, tpc_name, nf90_double, time_dimid, tpc_varid), "tp correlaion va&
   &r def error", error)
 
 
 
    if (error /= 0) return
 
    dimids3 = (/ x_dimid, y_dimid, time_dimid /)
    call check (nf90_def_var(ncid, pcp_name, nf90_double, dimids3, pcp_varid), "pcp var def error", &
   & error)
    call check (nf90_def_var(ncid, pop_name, nf90_double, dimids3, pop_varid), "pop var def error", &
   & error)
    call check (nf90_def_var(ncid, pcp_error_name, nf90_double, dimids3, pcp_error_varid), "pcp_err&
   &or var def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, tmean_name, nf90_double, dimids3, tmean_varid), "tmean var def e&
   &rror", error)
    call check (nf90_def_var(ncid, tmean_error_name, nf90_double, dimids3, tmean_error_varid), "tme&
   &an error var def error", error)
    call check (nf90_def_var(ncid, trange_name, nf90_double, dimids3, trange_varid), "trange var de&
   &f error", error)
    call check (nf90_def_var(ncid, trange_error_name, nf90_double, dimids3, trange_error_varid), "t&
   &range error var def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, pcp_name_2, nf90_double, dimids3, pcp_varid_2), "pcp var def err&
   &or", error)
    call check (nf90_def_var(ncid, pop_name_2, nf90_double, dimids3, pop_varid_2), "pop var def err&
   &or", error)
    call check (nf90_def_var(ncid, pcp_error_name_2, nf90_double, dimids3, pcp_error_varid_2), "pcp&
   &_error var def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, tmean_name_2, nf90_double, dimids3, tmean_varid_2), "tmean var d&
   &ef error", error)
    call check (nf90_def_var(ncid, tmean_error_name_2, nf90_double, dimids3, tmean_error_varid_2), &
   & "tmean error var def error", error)
    call check (nf90_def_var(ncid, trange_name_2, nf90_double, dimids3, trange_varid_2), "trange va&
   &r def error", error)
    call check (nf90_def_var(ncid, trange_error_name_2, nf90_double, dimids3, &
   & trange_error_varid_2), "trange error var def error", error)
    if (error /= 0) return
 
 
     !transformed mean, std variables, min, max of normalized y
    call check (nf90_def_var(ncid, y_mean_name, nf90_double, dimids3, ymean_varid), "y_mean var def&
   & error", error)
    call check (nf90_def_var(ncid, y_std_name, nf90_double, dimids3, ystd_varid), "y_std var def er&
   &ror", error)
    call check (nf90_def_var(ncid, y_max_name, nf90_double, dimids3, ymax_varid), "y_max var def er&
   &ror", error)
    call check (nf90_def_var(ncid, y_min_name, nf90_double, dimids3, ymin_varid), "y_min var def er&
   &ror", error)
    call check (nf90_def_var(ncid, y_stdall_name, nf90_double, dimids3, ystdall_varid), "y_std_all &
   &var def error", error)
 
    ! Add attributes.
 
     !long names
    call check (nf90_put_att(ncid, pcp_varid, long_name, pcp_long_name), "pcp long_name attribute e&
   &rror", error)
    call check (nf90_put_att(ncid, pop_varid, long_name, pop_long_name), "pcp long_name attribute e&
   &rror", error)
    call check (nf90_put_att(ncid, pcp_error_varid, long_name, pcp_error_long_name), "pcp_error lon&
   &g_name attribute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid, long_name, tmean_long_name), "tmean long_name attri&
   &bute error", error)
    call check (nf90_put_att(ncid, tmean_error_varid, long_name, tmean_error_long_name), "tmean lon&
   &g_name attribute error", error)
    call check (nf90_put_att(ncid, trange_varid, long_name, trange_long_name), "trange long_name at&
   &tribute error", error)
    call check (nf90_put_att(ncid, trange_error_varid, long_name, trange_error_long_name), "trange &
   &long_name attribute error", error)
 
    call check (nf90_put_att(ncid, pcp_varid_2, long_name, pcp_long_name_2), "pcp long_name attribu&
   &te error", error)
    call check (nf90_put_att(ncid, pop_varid_2, long_name, pop_long_name_2), "pcp long_name attribu&
   &te error", error)
    call check (nf90_put_att(ncid, pcp_error_varid_2, long_name, pcp_error_long_name_2), "pcp_error&
   & long_name attribute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid_2, long_name, tmean_long_name_2), "tmean long_name a&
   &ttribute error", error)
    call check (nf90_put_att(ncid, tmean_error_varid_2, long_name, tmean_error_long_name_2), "tmean&
   & long_name attribute error", error)
    call check (nf90_put_att(ncid, trange_varid_2, long_name, trange_long_name_2), "trange long_nam&
   &e attribute error", error)
    call check (nf90_put_att(ncid, trange_error_varid_2, long_name, trange_error_long_name_2), "tra&
   &nge long_name attribute error", error)
 
 
     !correlation variables
    call check (nf90_put_att(ncid, autoc_varid, long_name, autoc_long_name), "auto_corr long_name a&
   &ttribute error", error)
    call check (nf90_put_att(ncid, tpc_varid, long_name, tpc_long_name), "tp_corr long_name attribu&
   &te error", error)
 
     !transformed mean, std variables, min, max of normalized y
    call check (nf90_put_att(ncid, ymean_varid, long_name, y_mean_long_name), "ymean long_name attr&
   &ibute error", error)
    call check (nf90_put_att(ncid, ystd_varid, long_name, y_std_long_name), "ystd long_name attribu&
   &te error", error)
    call check (nf90_put_att(ncid, ystdall_varid, long_name, y_std_long_name), "ystd_all long_name &
   &attribute error", error)
    call check (nf90_put_att(ncid, ymax_varid, long_name, y_max_long_name), "ymax long_name attribu&
   &te error", error)
    call check (nf90_put_att(ncid, ymin_varid, long_name, y_min_long_name), "ymin long_name attribu&
   &te error", error)
     !units
    call check (nf90_put_att(ncid, lat_varid, units, lat_units), "lat units attribute error", &
   & error)
    call check (nf90_put_att(ncid, lon_varid, units, lon_units), "lon units attribute error", &
   & error)
    call check (nf90_put_att(ncid, alt_varid, units, alt_units), "alt units attribute error", &
   & error)
    call check (nf90_put_att(ncid, time_varid, units, time_units), "time units attribute error", &
   & error)
 
 
    call check (nf90_put_att(ncid, pcp_varid, units, pcp_units), "pcp units attribute error", &
   & error)
    call check (nf90_put_att(ncid, pop_varid, units, pop_units), "pcp units attribute error", &
   & error)
    call check (nf90_put_att(ncid, pcp_error_varid, units, pcp_error_units), "pcp_error units attri&
   &bute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid, units, tmean_units), "tmean units attribute error", &
   & error)
    call check (nf90_put_att(ncid, tmean_error_varid, units, tmean_error_units), "tmean_error units&
   & attribute error", error)
    call check (nf90_put_att(ncid, trange_varid, units, trange_units), "trange units attribute erro&
   &r", error)
    call check (nf90_put_att(ncid, trange_error_varid, units, trange_error_units), "trange_error un&
   &its attribute error", error)
 
    call check (nf90_put_att(ncid, pcp_varid_2, units, pcp_units), "pcp units attribute error", &
   & error)
    call check (nf90_put_att(ncid, pop_varid_2, units, pop_units), "pcp units attribute error", &
   & error)
    call check (nf90_put_att(ncid, pcp_error_varid_2, units, pcp_error_units), "pcp_error units att&
   &ribute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid_2, units, tmean_units), "tmean units attribute error&
   &", error)
    call check (nf90_put_att(ncid, tmean_error_varid_2, units, tmean_error_units), "tmean_error uni&
   &ts attribute error", error)
    call check (nf90_put_att(ncid, trange_varid_2, units, trange_units), "trange units attribute er&
   &ror", error)
    call check (nf90_put_att(ncid, trange_error_varid_2, units, trange_error_units), "trange_error &
   &units attribute error", error)
 
 
     !correlation variables
    call check (nf90_put_att(ncid, autoc_varid, units, autoc_units), "auto correlation units attrib&
   &ute error", error)
    call check (nf90_put_att(ncid, tpc_varid, units, tpc_units), "tp correlation units attribute er&
   &ror", error)
 
     !transformed mean,std variables, min, max of normalized y
    call check (nf90_put_att(ncid, ymean_varid, units, y_mean_units), "ymean units attribute error",&
   &  error)
    call check (nf90_put_att(ncid, ystd_varid, units, y_std_units), "ystd units attribute error", &
   & error)
    call check (nf90_put_att(ncid, ystdall_varid, units, y_std_units), "ystd_all units attribute er&
   &ror", error)
    call check (nf90_put_att(ncid, ymax_varid, units, y_max_units), "ymax units attribute error", &
   & error)
    call check (nf90_put_att(ncid, ymin_varid, units, y_min_units), "ymin units attribute error", &
   & error)
 
    if (error /= 0) return
 
 
 
 
     ! End define mode.
    call check (nf90_enddef(ncid), "end define mode error", error)
    if (error /= 0) return
 
    count2 = (/ inx, iny /)
    start2 = (/ 1, 1 /)
 
    call check (nf90_put_var(ncid, lat_varid, grdlat, start=start2, count=count2), "put lat error", &
   & error)
    call check (nf90_put_var(ncid, lon_varid, grdlon, start=start2, count=count2), "put lon error", &
   & error)
    call check (nf90_put_var(ncid, alt_varid, grdalt, start=start2, count=count2), "put alt error", &
   & error)
 
    trec = 1
    nrecs = 0
 
 
  else
 
     ! File already exists, get dim and var ids
    call check (nf90_inq_dimid(ncid, x_name, x_dimid), "x dim inq error", error)
    call check (nf90_inq_dimid(ncid, y_name, y_dimid), "y dim inq error", error)
    call check (nf90_inq_dimid(ncid, time_name, time_dimid), "time dim inq error", error)
    if (error /= 0) return
 
    call check (nf90_inq_varid(ncid, lat_name, lat_varid), "lat var inq error", error)
    call check (nf90_inq_varid(ncid, lon_name, lon_varid), "lon var inq error", error)
    call check (nf90_inq_varid(ncid, time_name, time_varid), "time var inq error", error)
    call check (nf90_inq_varid(ncid, pcp_name, pcp_varid), "pcp var inq error", error)
    call check (nf90_inq_varid(ncid, pop_name, pop_varid), "pop var inq error", error)
    call check (nf90_inq_varid(ncid, pcp_error_name, pcp_error_varid), "pcp_error var inq error", &
   & error)
    call check (nf90_inq_varid(ncid, tmean_name, tmean_varid), "tmean var inq error", error)
    call check (nf90_inq_varid(ncid, tmean_error_name, tmean_error_varid), "tmean error var inq err&
   &or", error)
    call check (nf90_inq_varid(ncid, trange_name, trange_varid), "trange var inq error", error)
    call check (nf90_inq_varid(ncid, trange_error_name, trange_error_varid), "trange error var inq &
   &error", error)
 
    call check (nf90_inq_varid(ncid, pcp_name_2, pcp_varid_2), "pcp var inq error", error)
    call check (nf90_inq_varid(ncid, pop_name_2, pop_varid_2), "pop var inq error", error)
    call check (nf90_inq_varid(ncid, pcp_error_name_2, pcp_error_varid_2), "pcp_error var inq error&
   &", error)
    call check (nf90_inq_varid(ncid, tmean_name_2, tmean_varid_2), "tmean var inq error", error)
    call check (nf90_inq_varid(ncid, tmean_error_name_2, tmean_error_varid_2), "tmean error var inq&
   & error", error)
    call check (nf90_inq_varid(ncid, trange_name_2, trange_varid_2), "trange var inq error", error)
    call check (nf90_inq_varid(ncid, trange_error_name_2, trange_error_varid_2), "trange error var &
   &inq error", error)
 
 
    call check (nf90_inq_varid(ncid, autoc_name, autoc_varid), "autoc var inq error", error)
    call check (nf90_inq_varid(ncid, tpc_name, tpc_varid), "tpc var inq error", error)
 
    call check (nf90_inq_varid(ncid, y_mean_name, ymean_varid), "ymean var inq error", error)
    call check (nf90_inq_varid(ncid, y_std_name, ystd_varid), "ystd var inq error", error)
    call check (nf90_inq_varid(ncid, y_stdall_name, ystdall_varid), "ystd_all var inq error", &
   & error)
    call check (nf90_inq_varid(ncid, y_max_name, ymax_varid), "ymax var inq error", error)
    call check (nf90_inq_varid(ncid, y_min_name, ymin_varid), "ymin var inq error", error)
    if (error /= 0) return
 
    call check (nf90_inquire_dimension(ncid, x_dimid, len=file_nx), "x dim len error", error)
    call check (nf90_inquire_dimension(ncid, y_dimid, len=file_ny), "y dim len error", error)
    call check (nf90_inquire_dimension(ncid, time_dimid, len=file_ntimes), "time dim len error", &
   & error)
    if (error /= 0) return
 
    if (nx /= file_nx .or. ny /= file_ny) then
      print *, "Error dimensions in output file do not match current run."
      error = 1
      return
    end if
 
    allocate (file_times(file_ntimes))
    call check (nf90_get_var(ncid, time_varid, file_times), "error getting file times list", error)
    if (error /= 0) return
 
    if (file_times(1) > times(n_times)) then !put data before everything in the file
      print *, "Error cannot add data before data already in output file. (functionality still to b&
     &e added)"
      error = 1
      return
    else
      if (file_times(file_ntimes) < times(1)) then !put data after everything in the file
        trec = file_ntimes + 1
      else ! at least some overlap
        do i = 1, file_ntimes, 1
          if (file_times(1) == times(1)) then
            trec = i
          end if
        end do
        if (trec == 0) then
          print *, "Error, confusion over data output record location."
          error = 1
          return
        else
          print *, "WARNING, overwriting data in output file, record ", trec, " to ", trec + &
         & n_times - 1
        end if
      end if
    end if
 
  end if
 
  count1 (1) = n_times
  start1 (1) = trec
  call check (nf90_put_var(ncid, time_varid, times, start=start1, count=count1), "put times error", &
 & error)
  if (error /= 0) return
 
 
  !correlation variables
  call check (nf90_put_var(ncid, autoc_varid, mean_autocorr, start=start1, count=count1), "put mean&
 & autocorrelation error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, tpc_varid, mean_tp_corr, start=start1, count=count1), "put mean t_&
 &p correlation error", error)
  if (error /= 0) return
 
 
 
 
  !3-d variables
  count3 = (/ inx, iny, n_times /)
  start3 = (/ 1, 1, trec /)
  call check (nf90_put_var(ncid, pcp_varid, real(pcp, kind(dp)), start=start3, count=count3), "put &
 &pcp error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, pop_varid, real(pop, kind(dp)), start=start3, count=count3), "put &
 &pop error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, pcp_error_varid, real(pcperror, kind(dp)), start=start3, &
 & count=count3), "put pcp_error error", error)
  if (error /= 0) return
 
 
  call check (nf90_put_var(ncid, tmean_varid, real(tmean, kind(dp)), start=start3, count=count3), "&
 &put tmean error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, tmean_error_varid, real(tmean_error, kind(dp)), start=start3, &
 & count=count3), "put tmean_error error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, trange_varid, real(trange, kind(dp)), start=start3, count=count3), &
 & "put trange error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, trange_error_varid, real(trange_error, kind(dp)), start=start3, &
 & count=count3), "put trange_error error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, pcp_varid_2, real(pcp_2, kind(dp)), start=start3, count=count3), "&
 &put pcp error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, pop_varid_2, real(pop_2, kind(dp)), start=start3, count=count3), "&
 &put pop error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, pcp_error_varid_2, real(pcperror_2, kind(dp)), start=start3, &
 & count=count3), "put pcp_error error", error)
  if (error /= 0) return
 
 
  call check (nf90_put_var(ncid, tmean_varid_2, real(tmean_2, kind(dp)), start=start3, &
 & count=count3), "put tmean error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, tmean_error_varid_2, real(tmean_error_2, kind(dp)), start=start3, &
 & count=count3), "put tmean_error error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, trange_varid_2, real(trange_2, kind(dp)), start=start3, &
 & count=count3), "put trange error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, trange_error_varid_2, real(trange_error_2, kind(dp)), &
 & start=start3, count=count3), "put trange_error error", error)
  if (error /= 0) return
 
 
 
!transformed mean,std variables, min & max of normalized y
  call check (nf90_put_var(ncid, ymean_varid, y_mean, start=start3, count=count3), "put ymean error&
 &", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, ystd_varid, y_std, start=start3, count=count3), "put ystd error", &
 & error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, ystdall_varid, y_std_all, start=start3, count=count3), "put ystd_a&
 &ll error", error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, ymax_varid, y_min, start=start3, count=count3), "put ymax error", &
 & error)
  if (error /= 0) return
 
  call check (nf90_put_var(ncid, ymin_varid, y_max, start=start3, count=count3), "put ymin error", &
 & error)
  if (error /= 0) return
 
  call check (nf90_close(ncid), "closing file error", error)
 
 
contains
  subroutine check (status, info, error)
    integer, intent (in) :: status
    character (len=*), intent (in) :: info
    integer, intent (out) :: error
 
    if (status /= nf90_noerr) then
      print *, trim (info) // ": " // trim (nf90_strerror(status))
      error = 1
    end if
  end subroutine check
 
end subroutine save_forcing_regression
!end subroutine save_precip
