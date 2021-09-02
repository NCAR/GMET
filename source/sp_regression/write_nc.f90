! ===== Writing related subroutines ======

! Subroutine to save the 
subroutine save_coefficients (n_vars, var_names, coefs, startdate, enddate, times, site_list, &
& stnid, stnlat, stnlon, stnelev, forecast, file, error)
  use netcdf
  use type
  implicit none
 
  character (len=100), intent (in) :: var_names (:), stnid (:)
  integer, intent (in) :: n_vars, forecast
  character (len=100), intent (in) :: startdate, enddate
  character (len=500), intent (in) :: file, site_list
  real (dp), intent (in) :: stnlat (:), stnlon (:), stnelev (:)
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
  character (len=*), parameter :: stn_elev_name = "station_elevation"
  character (len=*), parameter :: forecast_name = "forecast_hr"   ! used?
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
  character (len=*), parameter :: elev_units = "feet"
  character (len=*), parameter :: forecast_units = "hours"
  character (len=*), parameter :: date_units = "YYYYMMDD"
  character (len=*), parameter :: time_units = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len=*), parameter :: fill = "_FillValue"
   
  integer :: n_stns, n_chars, n_times
  integer :: ncid, stn_dimid, vars_dimid, char_dimid, time_dimid, rec_dimid
  integer :: coefs_varid, time_varid, vars_varid, forecast_varid, startdate_varid, enddate_varid
  integer :: stn_id_varid, stn_lat_varid, stn_lon_varid, stn_elev_varid
  integer :: charids (2)
  integer :: dimids (4)
  integer :: count4 (4), start4 (4), count2 (2), start2 (2), count1 (1), start1 (1)
  integer :: nstns, nvars, nrecs, ntimes, rec, i
 
  character (len=100) :: startdates, enddates
  integer, allocatable :: forecasts (:)
  integer :: forecast_arr (1)
 
  n_chars = 100
  n_stns = size (stnlat)
  n_times = size (times)

  error = nf90_open (file, nf90_write, ncid)
  if (error /= nf90_noerr) then
    print*, 'creating new coefficients file'
    error = 0
    ! Create the file.
    call check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
    if (error /= 0) return
    print*, 'Created coefficients file: '//trim(file)
 
    ! Define the dimensions.
    call check (nf90_def_dim(ncid, stn_name, n_stns, stn_dimid), "station dim def error", error)
    call check (nf90_def_dim(ncid, var_name, n_vars+1, vars_dimid), "variable dim def error", error)
    call check (nf90_def_dim(ncid, char_name, n_chars, char_dimid), "char dim def error", error)
    call check (nf90_def_dim(ncid, time_name, n_times, time_dimid), "time dim def error", error)
    call check (nf90_def_dim(ncid, rec_name, nf90_unlimited, rec_dimid), "rec dim def error", error)
    if (error /= 0) return
 
     ! Define the variables.
    charids = (/ char_dimid, vars_dimid /)
    call check (nf90_def_var(ncid, vars_name, nf90_char, charids, vars_varid), "variable_name var def error", error)
    charids = (/ char_dimid, stn_dimid /)
    call check (nf90_def_var(ncid, stn_id_name, nf90_char, charids, stn_id_varid), "station_id var def error", error)
    call check (nf90_def_var(ncid, stn_lat_name, nf90_double, stn_dimid, stn_lat_varid), "station_latitude var def error", error)
    call check (nf90_def_var(ncid, stn_lon_name, nf90_double, stn_dimid, stn_lon_varid), "station_longitude var def error", error)
    call check (nf90_def_var(ncid, stn_elev_name, nf90_double, stn_dimid, stn_elev_varid), "station_elevation var def error", error)
    if (error /= 0) return
 
    charids = (/ char_dimid, rec_dimid /)
    call check (nf90_def_var(ncid, startdate_name, nf90_char, charids, startdate_varid), "start_date var def error", error)
    call check (nf90_def_var(ncid, enddate_name, nf90_char, charids, enddate_varid), "end_date var def error", error)
    call check (nf90_def_var(ncid, forecast_name, nf90_int, rec_dimid, forecast_varid), "forecast var def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, time_var_name, nf90_double, time_dimid, time_varid), "time var def error", error)
    dimids = (/ stn_dimid, time_dimid, vars_dimid, rec_dimid /)
    call check (nf90_def_var(ncid, coefs_name, nf90_double, dimids, coefs_varid), "coefficient var def error", error)
    if (error /= 0) return
 
     ! Add attributes.
    call check (nf90_put_att(ncid, stn_lat_varid, units, lat_units), "station_lat units attribute error", error)
    call check (nf90_put_att(ncid, stn_lon_varid, units, lon_units), "station_lon units attribute error", error)
    call check (nf90_put_att(ncid, stn_elev_varid, units, elev_units), "station_elevation units attribute error", error)
    call check (nf90_put_att(ncid, startdate_varid, units, date_units), "start_date units attribute error", error)
    call check (nf90_put_att(ncid, enddate_varid, units, date_units), "end_date units attribute error", error)
    call check (nf90_put_att(ncid, forecast_varid, units, forecast_units), "forecase units attribute error", error)
    call check (nf90_put_att(ncid, time_varid, units, time_units), "time units attribute error", error)
    call check (nf90_put_att(ncid, coefs_varid, units, coefs_units), "coefficient units attribute error", error)
    call check (nf90_put_att(ncid, nf90_global, stn_file_name, site_list), "station_file_name global attribute error", error)
 
     ! End define mode.
    call check (nf90_enddef(ncid), "end define mode error", error)
    if (error /= 0) return
 
    call check (nf90_put_var(ncid, stn_lat_varid, stnlat), "put station_lat error", error)
    call check (nf90_put_var(ncid, stn_lon_varid, stnlon), "put station_lon error", error)
    call check (nf90_put_var(ncid, stn_elev_varid, stnelev), "put station_elevation error", error)
    call check (nf90_put_var(ncid, time_varid, times), "put times error", error)
 
    count2 = (/ 1, 1 /)
    start2 = (/ 1, 1 /)
    do i = 1, n_stns, 1
      count2 (1) = len (trim(stnid(i)))
      start2 (2) = i
      call check (nf90_put_var(ncid, stn_id_varid, stnid(i), start=start2, count=count2), "put station_id error", error)
      if (error /= 0) return
    end do
 
    count2 (1) = len (trim(constant_name))
    start2 (2) = 1
    call check (nf90_put_var(ncid, vars_varid, constant_name, start=start2, count=count2), "put variable_name error", error)
    if (error /= 0) return
    do i = 1, n_vars, 1
      count2 (1) = len (trim(var_names(i)))
      start2 (2) = i + 1
      call check (nf90_put_var(ncid, vars_varid, var_names(i), start=start2, count=count2), "put variable_name error", error)
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
    call check (nf90_inq_varid(ncid, stn_lat_name, stn_lat_varid), "station_latitude var inq error", error)
    call check (nf90_inq_varid(ncid, stn_lon_name, stn_lon_varid), "station_longitude var inq error", error)
    call check (nf90_inq_varid(ncid, stn_elev_name, stn_elev_varid), "station_elevation var inq error", error)
    call check (nf90_inq_varid(ncid, startdate_name, startdate_varid), "start_date var inq error", error)
    call check (nf90_inq_varid(ncid, enddate_name, enddate_varid), "end_date var inq error", error)
    call check (nf90_inq_varid(ncid, forecast_name, forecast_varid), "forecast var inq error", error)
    call check (nf90_inq_varid(ncid, time_var_name, time_varid), "time var inq error", error)
    call check (nf90_inq_varid(ncid, coefs_name, coefs_varid), "coefficient var inq error", error)
    if (error /= 0) return
 
     ! Verify Dimensions match
    call check (nf90_inquire_dimension(ncid, stn_dimid, len=nstns), "station dim len error", error)
    call check (nf90_inquire_dimension(ncid, vars_dimid, len=nvars), "variable dim len error", error)
    call check (nf90_inquire_dimension(ncid, time_dimid, len=ntimes), "time dim len error", error)
    call check (nf90_inquire_dimension(ncid, rec_dimid, len=nrecs), "run dim len error", error)
    if (error /= 0) return
 
    if (n_stns /= nstns .or. nvars /= n_vars+1 .or. ntimes /= n_times) then
      print *, "Error dimensions in output file do not match current run."
      error = 1
      return
    end if
 
    allocate (forecasts(nrecs))
 
    call check (nf90_get_var(ncid, forecast_varid, forecasts), "error getting file forecasts", error)
    if (error /= 0) return
 
    rec = nrecs + 1
    do i = 1, nrecs, 1
      count2 = (/ 100, 1 /)
      start2 = (/ 1, i /)
      call check (nf90_get_var(ncid, startdate_varid, startdates, start=start2, count=count2), "error getting file startdates", error)
      call check (nf90_get_var(ncid, enddate_varid, enddates, start=start2, count=count2), "error getting file enddates", error)
      if (forecasts(i) == forecast .and. startdates(1:8) == startdate(1:8) .and. enddates(1:8) == enddate(1:8)) then
        print *, "WARNING, overwriting data in output file, record ", i
        rec = i
      end if
    end do
 
  end if   ! end IF clause on whether to open existing file or append
 
  count2 = (/ len (trim(startdate)), 1 /)
  start2 = (/ 1, rec /)
  call check (nf90_put_var(ncid, startdate_varid, startdate, start=start2, count=count2), "put start_date error", error)
  if (error /= 0) return
 
  count2 = (/ len (trim(enddate)), 1 /)
  start2 = (/ 1, rec /)
  call check (nf90_put_var(ncid, enddate_varid, enddate, start=start2, count=count2), "put end_date error", error)
  if (error /= 0) return
 
  count1 (1) = 1
  start1 (1) = rec
  forecast_arr = forecast
  call check (nf90_put_var(ncid, forecast_varid, forecast_arr, start=start1, count=count1), "put forecast error", error)
  if (error /= 0) return
 
  count4 = (/ n_stns, n_times, n_vars + 1, 1 /)
  start4 = (/ 1, 1, 1, rec /)
  forecast_arr = forecast
  call check (nf90_put_var(ncid, coefs_varid, coefs, start=start4, count=count4), "put coefficients error", error)
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
 
 
! ==== subroutine save_forcing_regression:  saves the forcing regression parameters for use by the SCRF 
!      program in generating ensemble forcings
subroutine save_forcing_regression (pcp, pop, pcperror, obs_max_pcp, tmean, tmean_error, trange, trange_error, &
& nx, ny, grdlat, grdlon, grdelev, times, mean_autocorr, mean_tp_corr, &
& file, error, pcp_2, pop_2, pcperror_2, tmean_2, tmean_error_2, trange_2, trange_error_2)

  use netcdf
  use type
  implicit none
 
  real (sp), intent (in) :: pcp(:,:), pop(:,:), pcperror(:,:)
  real (sp), intent (in) :: tmean(:,:), tmean_error(:,:), trange(:,:), trange_error(:,:)
  real (sp), intent (in) :: pcp_2(:,:), pop_2(:,:), pcperror_2(:,:)
  real (sp), intent (in) :: tmean_2(:,:), tmean_error_2(:,:), trange_2(:,:), trange_error_2(:,:)
 
  integer (i4b), intent (in) :: nx, ny
  real (dp), intent (in) :: grdlat(:), grdlon(:), grdelev(:)
  real (dp), intent (in) :: times(:)
  real (dp), intent (in) :: mean_autocorr(:), mean_tp_corr(:)
 
  real (dp), intent (in) :: obs_max_pcp(:, :)
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
  character (len=*), parameter :: elev_name = "elevation"
  character (len=*), parameter :: pcp_name = "pcp"
  character (len=*), parameter :: pop_name = "pop"
  character (len=*), parameter :: pcp_error_name = "pcp_error"
  character (len=*), parameter :: tmean_name = "tmean"
  character (len=*), parameter :: tmean_error_name = "tmean_error"
  character (len=*), parameter :: trange_name = "trange"
  character (len=*), parameter :: trange_error_name = "trange_error"
  character (len=*), parameter :: obs_max_pcp_name = "obs_max_pcp"
  character (len=*), parameter :: pcp_name_2 = "pcp_2"
  character (len=*), parameter :: pop_name_2 = "pop_2"
  character (len=*), parameter :: pcp_error_name_2 = "pcp_error_2"
  character (len=*), parameter :: tmean_name_2 = "tmean_2"
  character (len=*), parameter :: tmean_error_name_2 = "tmean_error_2"
  character (len=*), parameter :: trange_name_2 = "trange_2"
  character (len=*), parameter :: trange_error_name_2 = "trange_error_2"
 
  ! Long names 
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
  character (len=*), parameter :: obs_max_pcp_long_name = "Maximum obseved precipitation (transformed)"
  character (len=*), parameter :: pcp_long_name_2 = "estimated precip in normal space (no slope)"
  character (len=*), parameter :: pop_long_name_2 = "probability of precipitation occurrence (no slope)"
  character (len=*), parameter :: pcp_error_long_name_2 = "error in estimated precip (no slope)"
  character (len=*), parameter :: tmean_long_name_2 = "estimated daily mean temperature (no slope)"
  character (len=*), parameter :: tmean_error_long_name_2 = "error in estimated daily mean temp (no slope)"
  character (len=*), parameter :: trange_long_name_2 = "estimated diurnal range (no slope)"
  character (len=*), parameter :: trange_error_long_name_2 = "error in estimated diurnal range (no slope)"
 
  ! Units
  character (len=*), parameter :: units = "units"
  character (len=*), parameter :: pcp_units = "transformed precip"
  character (len=*), parameter :: pop_units = "unitless"
  character (len=*), parameter :: autoc_units = "unitless"
  character (len=*), parameter :: tpc_units = "unitless"
  character (len=*), parameter :: pcp_error_units = "transformed precip"
  character (len=*), parameter :: tmean_units = "deg_C"
  character (len=*), parameter :: trange_units = "deg_C"
  character (len=*), parameter :: tmean_error_units = "deg_C"
  character (len=*), parameter :: trange_error_units = "deg_C"
  character (len=*), parameter :: obs_max_pcp_units = "transformed precip"
 
  character (len=*), parameter :: lat_units = "degrees_north"
  character (len=*), parameter :: lon_units = "degrees_east"
  character (len=*), parameter :: elev_units = "meters"
  character (len=*), parameter :: time_units = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len=*), parameter :: fill = "_FillValue"
  
  integer, parameter           :: append_choice = 0  ! 0 = clobber existing output file
                                                     ! 1 = try to append to existing output file

  ! local variables
  real (dp), allocatable :: file_times (:)
 
  integer :: n_chars, n_times, inx, iny
  integer :: ncid, x_dimid, y_dimid, time_dimid
  integer :: lat_varid, lon_varid, autoc_varid, elev_varid, time_varid, pcp_varid, pop_varid, &
             & pcp_error_varid, tpc_varid
  integer :: tmean_varid, tmean_error_varid, trange_varid, trange_error_varid
  integer :: pcp_varid_2, pop_varid_2, pcp_error_varid_2, obs_max_pcp_varid
  integer :: tmean_varid_2, tmean_error_varid_2, trange_varid_2, trange_error_varid_2

  integer :: count1(1), start1(1), count2(2), start2(2), count3(3), start3(3), dimids2(2), dimids3(3)
  integer :: trec, nrecs, file_nx, file_ny, file_ntimes, i
 
  trec = 0
  n_chars = 100
  n_times = size (times)
  inx = nx
  iny = ny
 
  if (size(grdlat) /= inx*iny) then
    print *, "ERROR, gridlat is not the same size as inx*iny "
    print*, "These are, respectively: ", size(grdlat), inx, iny
  end if
 
  ! use parameter setting (above) to override file appending (if selected)
  ! could be handled better in nested if
  if(append_choice == 1) then
    print*, ' -- trying to append existing output file: ', file
    error = nf90_open (file, nf90_write, ncid)
  else
    ! don't try to open a file, and set error to anything but nf90_noerr (0)
    error = 1
  end if
  if (error /= nf90_noerr .or. append_choice == 0) then
    print*, ' -- creating new output file: ', file
    error = 0

    ! Create the file.
    call check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
    if (error /= 0) return
 
    ! Define the dimensions
    call check (nf90_def_dim(ncid, y_name, iny, y_dimid), "y dim def error", error)
    call check (nf90_def_dim(ncid, x_name, inx, x_dimid), "x dim def error", error)
    call check (nf90_def_dim(ncid, time_name, nf90_unlimited, time_dimid), "time dim def error", error)
    if (error /= 0) return
 
    ! Define the variables.
    dimids2 = (/ x_dimid, y_dimid /)

    ! AW
    call check (nf90_def_var(ncid, lat_name, nf90_float, dimids2, lat_varid), "lat var def error", error)
    call check (nf90_def_var(ncid, lon_name, nf90_float, dimids2, lon_varid), "lon var def error", error)
    call check (nf90_def_var(ncid, elev_name, nf90_float, dimids2, elev_varid), "lon var def error", error)
    call check (nf90_def_var(ncid, time_name, nf90_float, time_dimid, time_varid), "time var def error", error)
 
    ! Correlation variables
    call check (nf90_def_var(ncid, autoc_name, nf90_float, time_dimid, autoc_varid), "auto correlation var def error", error)
    call check (nf90_def_var(ncid, tpc_name, nf90_float, time_dimid, tpc_varid), "tp correlation var def error", error)
    if (error /= 0) return
 
    dimids3 = (/ x_dimid, y_dimid, time_dimid /)
    call check (nf90_def_var(ncid, pcp_name, nf90_float, dimids3, pcp_varid), "pcp var def error", error)
    call check (nf90_def_var(ncid, pop_name, nf90_float, dimids3, pop_varid), "pop var def error", error)
    call check (nf90_def_var(ncid, pcp_error_name, nf90_float, dimids3, pcp_error_varid), "pcp_error var def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, tmean_name, nf90_float, dimids3, tmean_varid), "tmean var def error", error)
    call check (nf90_def_var(ncid, tmean_error_name, nf90_float, dimids3, tmean_error_varid), "tmean error var def error", error)
    call check (nf90_def_var(ncid, trange_name, nf90_float, dimids3, trange_varid), "trange var def error", error)
    call check (nf90_def_var(ncid, trange_error_name, nf90_float, dimids3, trange_error_varid), "trange error var def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, pcp_name_2, nf90_float, dimids3, pcp_varid_2), "pcp var def error", error)
    call check (nf90_def_var(ncid, pop_name_2, nf90_float, dimids3, pop_varid_2), "pop var def error", error)
    call check (nf90_def_var(ncid, pcp_error_name_2, nf90_float, dimids3, pcp_error_varid_2), "pcp_error var def error", error)
    if (error /= 0) return
 
    call check (nf90_def_var(ncid, tmean_name_2, nf90_float, dimids3, tmean_varid_2), "tmean var def error", error)
    call check (nf90_def_var(ncid, tmean_error_name_2, nf90_float, dimids3, tmean_error_varid_2), "tmean error var def error", error)
    call check (nf90_def_var(ncid, trange_name_2, nf90_float, dimids3, trange_varid_2), "trange var def error", error)
    call check (nf90_def_var(ncid, trange_error_name_2, nf90_float, dimids3, trange_error_varid_2), "trange error var def error", error)
   call check (nf90_def_var(ncid, obs_max_pcp_name, nf90_double, dimids3, obs_max_pcp_varid), "obs_max_pcp var def error", error)
    if (error /= 0) return
 
    ! Add attributes.
 
     ! long names
    call check (nf90_put_att(ncid, pcp_varid, long_name, pcp_long_name), "pcp long_name attribute error", error)
    call check (nf90_put_att(ncid, pop_varid, long_name, pop_long_name), "pcp long_name attribute error", error)
    call check (nf90_put_att(ncid, pcp_error_varid, long_name, pcp_error_long_name), "pcp_error long_name attribute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid, long_name, tmean_long_name), "tmean long_name attribute error", error)
    call check (nf90_put_att(ncid, tmean_error_varid, long_name, tmean_error_long_name), "tmean long_name attribute error", error)
    call check (nf90_put_att(ncid, trange_varid, long_name, trange_long_name), "trange long_name attribute error", error)
    call check (nf90_put_att(ncid, trange_error_varid, long_name, trange_error_long_name), "trange long_name attribute error", error)
 
    call check (nf90_put_att(ncid, pcp_varid_2, long_name, pcp_long_name_2), "pcp long_name attribute error", error)
    call check (nf90_put_att(ncid, pop_varid_2, long_name, pop_long_name_2), "pcp long_name attribute error", error)
    call check (nf90_put_att(ncid, pcp_error_varid_2, long_name, pcp_error_long_name_2), "pcp_error long_name attribute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid_2, long_name, tmean_long_name_2), "tmean long_name attribute error", error)
    call check (nf90_put_att(ncid, tmean_error_varid_2, long_name, tmean_error_long_name_2), "tmean long_name attribute error", error)
    call check (nf90_put_att(ncid, trange_varid_2, long_name, trange_long_name_2), "trange long_name attribute error", error)
    call check (nf90_put_att(ncid, trange_error_varid_2, long_name, trange_error_long_name_2), "trange long_name attribute error", error)
 
    call check (nf90_put_att(ncid, obs_max_pcp_varid, long_name, obs_max_pcp_long_name), "obs_max_pcp long_name attribute error", error)

    ! correlation variables
    call check (nf90_put_att(ncid, autoc_varid, long_name, autoc_long_name), "auto_corr long_name attribute error", error)
    call check (nf90_put_att(ncid, tpc_varid, long_name, tpc_long_name), "tp_corr long_name attribute error", error)
   
    ! units
    call check (nf90_put_att(ncid, lat_varid, units, lat_units), "lat units attribute error", error)
    call check (nf90_put_att(ncid, lon_varid, units, lon_units), "lon units attribute error", error)
    call check (nf90_put_att(ncid, elev_varid, units, elev_units), "elev units attribute error", error)
    call check (nf90_put_att(ncid, time_varid, units, time_units), "time units attribute error", error)
 
    call check (nf90_put_att(ncid, pcp_varid, units, pcp_units), "pcp units attribute error", error)
    call check (nf90_put_att(ncid, pop_varid, units, pop_units), "pcp units attribute error", error)
    call check (nf90_put_att(ncid, pcp_error_varid, units, pcp_error_units), "pcp_error units attribute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid, units, tmean_units), "tmean units attribute error", error)
    call check (nf90_put_att(ncid, tmean_error_varid, units, tmean_error_units), "tmean_error units attribute error", error)
    call check (nf90_put_att(ncid, trange_varid, units, trange_units), "trange units attribute error", error)
    call check (nf90_put_att(ncid, trange_error_varid, units, trange_error_units), "trange_error units attribute error", error)
 
    call check (nf90_put_att(ncid, pcp_varid_2, units, pcp_units), "pcp units attribute error", error)
    call check (nf90_put_att(ncid, pop_varid_2, units, pop_units), "pcp units attribute error", error)
    call check (nf90_put_att(ncid, pcp_error_varid_2, units, pcp_error_units), "pcp_error units attribute error", error)
 
    call check (nf90_put_att(ncid, tmean_varid_2, units, tmean_units), "tmean units attribute error", error)
    call check (nf90_put_att(ncid, tmean_error_varid_2, units, tmean_error_units), "tmean_error units attribute error", error)
    call check (nf90_put_att(ncid, trange_varid_2, units, trange_units), "trange units attribute error", error)
    call check (nf90_put_att(ncid, trange_error_varid_2, units, trange_error_units), "trange_error units attribute error", error)
 
    ! correlation variables
    call check (nf90_put_att(ncid, autoc_varid, units, autoc_units), "auto correlation units attribute error", error)
    call check (nf90_put_att(ncid, tpc_varid, units, tpc_units), "tp correlation units attribute error", error)
 
    call check (nf90_put_att(ncid, obs_max_pcp_varid, units, obs_max_pcp_units), "obs_max_pcp units attribute error", error)
    if (error /= 0) return

    ! End define mode.
    call check (nf90_enddef(ncid), "end define mode error", error)
    if (error /= 0) return
 
    count2 = (/ inx, iny /)
    start2 = (/ 1, 1 /)
 
    call check (nf90_put_var(ncid, lat_varid, grdlat, start=start2, count=count2), "put lat grd error", error)
    call check (nf90_put_var(ncid, lon_varid, grdlon, start=start2, count=count2), "put lon grd error", error)
    call check (nf90_put_var(ncid, elev_varid, grdelev, start=start2, count=count2), "put elev grd error", error)
 
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
    call check (nf90_inq_varid(ncid, pcp_error_name, pcp_error_varid), "pcp_error var inq error", error)
    call check (nf90_inq_varid(ncid, tmean_name, tmean_varid), "tmean var inq error", error)
    call check (nf90_inq_varid(ncid, tmean_error_name, tmean_error_varid), "tmean error var inq error", error)
    call check (nf90_inq_varid(ncid, trange_name, trange_varid), "trange var inq error", error)
    call check (nf90_inq_varid(ncid, trange_error_name, trange_error_varid), "trange error var inq error", error)
 
    call check (nf90_inq_varid(ncid, pcp_name_2, pcp_varid_2), "pcp var inq error", error)
    call check (nf90_inq_varid(ncid, pop_name_2, pop_varid_2), "pop var inq error", error)
    call check (nf90_inq_varid(ncid, pcp_error_name_2, pcp_error_varid_2), "pcp_error var inq error", error)
    call check (nf90_inq_varid(ncid, tmean_name_2, tmean_varid_2), "tmean var inq error", error)
    call check (nf90_inq_varid(ncid, tmean_error_name_2, tmean_error_varid_2), "tmean error var inq error", error)
    call check (nf90_inq_varid(ncid, trange_name_2, trange_varid_2), "trange var inq error", error)
    call check (nf90_inq_varid(ncid, trange_error_name_2, trange_error_varid_2), "trange error var inq error", error)

    call check (nf90_inq_varid(ncid, autoc_name, autoc_varid), "autoc var inq error", error)
    call check (nf90_inq_varid(ncid, tpc_name, tpc_varid), "tpc var inq error", error)
    if (error /= 0) return
 
    call check (nf90_inq_varid(ncid, obs_max_pcp_name, obs_max_pcp_varid), "obs_max_pcp var inq error", error)
    if (error /= 0) return
 
    call check (nf90_inquire_dimension(ncid, x_dimid, len=file_nx), "x dim len error", error)
    call check (nf90_inquire_dimension(ncid, y_dimid, len=file_ny), "y dim len error", error)
    call check (nf90_inquire_dimension(ncid, time_dimid, len=file_ntimes), "time dim len error", error)
    if (error /= 0) return
 
    if (nx /= file_nx .or. ny /= file_ny) then
      print *, "ERROR: dimensions in output file do not match current run."
      error = 1
      return
    end if
 
    allocate (file_times(file_ntimes))
    call check (nf90_get_var(ncid, time_varid, file_times), "error getting file times list", error)
    if (error /= 0) return
 
    if (file_times(1) > times(n_times)) then ! put data before everything in the file
      print *, "ERROR: cannot add data before data already in output file. (functionality still to be added)"
      error = 1
      return
    else
      if (file_times(file_ntimes) < times(1)) then ! put data after everything in the file
        trec = file_ntimes + 1
      else ! at least some overlap
        do i = 1, file_ntimes, 1
          if (file_times(i) == times(1)) then
            trec = i
          end if
        end do
        if (trec == 0) then
          print *, "Error, confusion over data output record location."
          error = 1
          return
        else
          print *, "WARNING, overwriting data in output file, record ", trec, " to ", trec + n_times - 1
        end if
      end if
    end if
 
  end if
 
  count1 (1) = n_times
  start1 (1) = trec
  call check (nf90_put_var(ncid, time_varid, times, start=start1, count=count1), "put times error", error)
  if (error /= 0) return
 
  !correlation variables
  call check (nf90_put_var(ncid, autoc_varid, mean_autocorr, start=start1, count=count1), "put mean autocorrelation error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, tpc_varid, mean_tp_corr, start=start1, count=count1), "put mean t_p correlation error", error)
  if (error /= 0) return
 
  !3-d variables
  count3 = (/ inx, iny, n_times /)
  start3 = (/ 1, 1, trec /)
  call check (nf90_put_var(ncid, pcp_varid, real(pcp, kind(dp)), start=start3, count=count3), "put pcp error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, pop_varid, real(pop, kind(dp)), start=start3, count=count3), "put pop error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, pcp_error_varid, real(pcperror, kind(dp)), start=start3, count=count3), "put pcp_error error", error)
  if (error /= 0) return

  call check (nf90_put_var(ncid, obs_max_pcp_varid, obs_max_pcp, start=start3, count=count3), "put obs_max_pcp error", error)
  if (error /= 0) return

  call check (nf90_put_var(ncid, tmean_varid, real(tmean, kind(dp)), start=start3, count=count3), "put tmean error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, tmean_error_varid, real(tmean_error, kind(dp)), start=start3, count=count3), "put tmean_error error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, trange_varid, real(trange, kind(dp)), start=start3, count=count3), "put trange error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, trange_error_varid, real(trange_error, kind(dp)), start=start3, count=count3), "put trange_error error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, pcp_varid_2, real(pcp_2, kind(dp)), start=start3, count=count3), "put pcp error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, pop_varid_2, real(pop_2, kind(dp)), start=start3, count=count3), "put pop error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, pcp_error_varid_2, real(pcperror_2, kind(dp)), start=start3, count=count3), "put pcp_error error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, tmean_varid_2, real(tmean_2, kind(dp)), start=start3, count=count3), "put tmean error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, tmean_error_varid_2, real(tmean_error_2, kind(dp)), start=start3, count=count3), "put tmean_error error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, trange_varid_2, real(trange_2, kind(dp)), start=start3, count=count3), "put trange error", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, trange_error_varid_2, real(trange_error_2, kind(dp)), start=start3, count=count3), "put trange_error error", error)
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
