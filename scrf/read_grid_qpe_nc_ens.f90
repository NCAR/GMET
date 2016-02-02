subroutine read_grid_qpe_nc_ens(file_name, var_name, var, lats, lons, auto_corr, tp_corr, &
                            times, ntimes, error)
  use netcdf
  use nrtype
  implicit none
  
  character (len = *), intent(in) :: file_name
  character (len = *), intent(in) :: var_name
  real(DP), allocatable, intent(out) :: var(:, :, :)
  real(DP), allocatable, intent(out) :: lats(:,:), lons(:,:)
  real(DP), allocatable, intent(out) :: times(:)
  real(DP), allocatable, intent(out) :: auto_corr(:)
  real(DP), allocatable, intent(out) :: tp_corr(:)


  integer, intent(out) :: error
  integer(I4B),intent(out) :: ntimes

  integer(I4B) :: nlats, nlons
  integer(I4B) :: ncid, varid, latid, lonid, timeid,autoc_id,tpc_id
  integer(I4B) :: tmean_mean_id,tmean_stdev_id,trange_mean_id,trange_stdev_id
  integer(I4B) :: ndims, lon_dimid, lat_dimid, time_dimid
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: TIME_NAME = "time"
  character (len = *), parameter :: autoc_name = "auto_corr"
  character (len = *), parameter :: tpc_name = "tp_corr"


  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"

  integer(I4B), dimension(nf90_max_var_dims) :: dimIds


  ! Loop indices
  integer(I4B) :: lvl, lat, lon, rec, i, time_index

!  integer(I4B) :: start(3), count(3)

  ! To check the units attributes.
  character*80 pres_units_in, temp_units_in
  character*80 lat_units_in, lon_units_in

  error = 0
  time_index = 0

  ! Open the file. 
  call check( nf90_open(file_name, nf90_nowrite, ncid), "File open error", error)
  if(error /= 0) return

  ! Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, LAT_NAME, latid), "Latitutde name error", error)
  call check ( nf90_inquire_variable(ncid, latid, ndims = ndims, dimids = dimIds), "Latitutde inq error", error)
  if(error /= 0 .OR. ndims /= 2) return
  lat_dimid = dimIds(1)

  call check( nf90_inq_varid(ncid, LON_NAME, lonid), "Longitude name error", error)
  call check ( nf90_inquire_variable(ncid, lonid, ndims = ndims, dimids = dimIds), "Longitude inq error", error)
  if(error /= 0 .OR. ndims /= 2) return
  lon_dimid = dimIds(2)

  call check( nf90_inq_varid(ncid, TIME_NAME, timeid), "Time name error", error)
  call check ( nf90_inquire_variable(ncid, timeid, ndims = ndims, dimids = dimIds), "Time inq error", error)
  if(error /= 0 .OR. ndims /= 1) return

  call check( nf90_inq_varid(ncid, autoc_name, autoc_id), "Autocorrelation name error", error)
  call check ( nf90_inquire_variable(ncid, autoc_id, ndims = ndims, dimids = dimIds), "Autocorrelation inq error", error)
  if(error /= 0 .OR. ndims /= 1) return

  call check( nf90_inq_varid(ncid, tpc_name, tpc_id), "Temp-precip correlation name error", error)
  call check ( nf90_inquire_variable(ncid, tpc_id, ndims = ndims, dimids = dimIds), &
               "Temp-precip correlation inq error", error)
  if(error /= 0 .OR. ndims /= 1) return




  call check ( nf90_inquire_dimension(ncid, lat_dimid, len = nlats), "Latitutde dim error", error)
  call check ( nf90_inquire_dimension(ncid, lon_dimid, len = nlons), "Longitude dim error", error)
  if(error /=0 ) return

  allocate(lats(nlats,nlons), lons(nlats,nlons))

  ! Read the latitude and longitude data.
  call check( nf90_get_var(ncid, latid, lats), "Latitutde read error", error)
  call check( nf90_get_var(ncid, lonid, lons), "Longitude read error", error)

  call check( nf90_inq_varid(ncid, var_name, varid), "Variable name error", error)
  call check ( nf90_inquire_variable(ncid, varid, ndims = ndims, dimids = dimIds), "Variable inq error", error)

  if(error == 0 .AND. ndims == 3 .AND. dimIds(2) == latid .AND. dimIds(1) == lonid) then

     time_dimid = dimIds(3)
     
     call check ( nf90_inquire_dimension(ncid, time_dimid, len = ntimes), "Time dim error", error)

     allocate(times(ntimes))
     allocate(auto_corr(ntimes))
     allocate(tp_corr(ntimes))

     
     call check( nf90_get_var(ncid, timeid, times), "Time read error", error)
     if(error /= 0) return

     call check( nf90_get_var(ncid, autoc_id, auto_corr), "Autocorrelation read error", error)
     if(error /= 0) return

     call check( nf90_get_var(ncid, tpc_id, tp_corr), "Temp-precip correlation read error", error)
     if(error /= 0) return   


     allocate(var(nlats, nlons, ntimes))


     call check( nf90_get_var(ncid, varid, var), "Variable read error", error)

     if(error /= 0) then
        deallocate(var)
     endif
  endif

  ! Close the file. 
  i = nf90_close(ncid)

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

end subroutine read_grid_qpe_nc_ens

