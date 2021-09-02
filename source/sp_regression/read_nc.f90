subroutine read_nc_file (file_name, valid_time, var_name, var, lats, lons, error)
  use netcdf
  use type
  implicit none
 
  character (len=*), intent (in) :: file_name
  character (len=*), intent (in) :: var_name
  real (dp), intent (in) :: valid_time
  real (dp), allocatable, intent (out) :: var (:, :, :)
  real (dp), allocatable, intent (out) :: lats (:), lons (:)
  integer, intent (out) :: error
 
  integer (i4b) :: nlats, nlons, ntimes
  integer (i4b) :: ncid, varid, latid, lonid, timeid
  integer (i4b) :: ndims, lon_dimid, lat_dimid, time_dimid
  character (len=*), parameter :: lat_name = "latitude"
  character (len=*), parameter :: lon_name = "longitude"
  character (len=*), parameter :: time_name = "time"
 
  integer (i4b), dimension (nf90_max_var_dims) :: dimids
  real (dp), allocatable :: times (:), fcst_times (:)  ! modified AJN 2/4/2014
 
  ! Loop indices
  integer (i4b) :: i, time_index
 
  ! indicies for aggregating variables between forecast times
  integer (i4b) :: start_index
 
  integer (i4b) :: start (3), count (3)
 
  real (dp), allocatable :: accum_var (:, :, :)
 
  error = 0
  time_index = 0
 
  ! Open the file.
  call check (nf90_open(file_name, nf90_nowrite, ncid), "File open error", error)
  if (error /= 0) return
 
  ! Get the varids of the latitude and longitude coordinate variables.
  call check (nf90_inq_varid(ncid, lat_name, latid), "Latitude name error", error)
  call check (nf90_inquire_variable(ncid, latid, ndims=ndims, dimids=dimids), "Latitutde inq error",&
 &  error)
  if (error /= 0 .or. ndims /= 1) return
  lat_dimid = dimids (1)
 
  call check (nf90_inq_varid(ncid, lon_name, lonid), "Longitude name error", error)
  call check (nf90_inquire_variable(ncid, lonid, ndims=ndims, dimids=dimids), "Longitude inq error",&
 &  error)
  if (error /= 0 .or. ndims /= 1) return
  lon_dimid = dimids (1)
 
  call check (nf90_inq_varid(ncid, time_name, timeid), "Time name error", error)
  call check (nf90_inquire_variable(ncid, timeid, ndims=ndims, dimids=dimids), "Time inq error", &
 & error)
  if (error /= 0 .or. ndims /= 1) return
 
  call check (nf90_inquire_dimension(ncid, lat_dimid, len=nlats), "Latitutde dim error", error)
  call check (nf90_inquire_dimension(ncid, lon_dimid, len=nlons), "Longitude dim error", error)
  if (error /= 0) return
 
  allocate (lats(nlats), lons(nlons))
  ! Read the latitude and longitude data.
  call check (nf90_get_var(ncid, latid, lats), "Latitutde read error", error)
  call check (nf90_get_var(ncid, lonid, lons), "Longitude read error", error)
 
 
  call check (nf90_inq_varid(ncid, var_name, varid), "Variable name error", error)
  call check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), "Variable inq error", &
 & error)
 
  if (error == 0 .and. ndims == 3 .and. dimids(2) == latid .and. dimids(1) == lonid) then
 
    time_dimid = dimids (3)
 
    call check (nf90_inquire_dimension(ncid, time_dimid, len=ntimes), "Time dim error", error)
 
    allocate (times(ntimes))
    allocate (fcst_times(ntimes))
 
    call check (nf90_get_var(ncid, timeid, times), "Time read error", error)
    if (error /= 0) return
 
    fcst_times = times
 
    do i = 1, ntimes, 1
      if (times(i) == valid_time) time_index = i
    end do
    deallocate (times)
    if (time_index == 0) then
      print ("(A,F11.0,A)"), "Error: time (", valid_time, ") not found in this file. "
      error = 1
      deallocate (lats)
      deallocate (lons)
      return
    end if
 
    if (trim(var_name) == 'APCP_ens_mean_surface' .or. trim(var_name) == 'DSWRF_ens_mean_surface') &
   & then
       ! find fcst_time that is ~24 hours prior
       ! just find closest one that is either 24-hrs or less
 
      do i = 1, ntimes, 1
        if (fcst_times(i) >= valid_time-86400) then
          start_index = i
          print *, start_index, time_index
          exit
        end if
      end do
 
      start = (/ 1, 1, start_index /)
      count = (/ nlons, nlats, 1 /)
      allocate (var(nlons, nlats, 1))
      allocate (accum_var(nlons, nlats, 1))
      accum_var = 0.0
 
      do i = 1, (time_index-start_index), 1
        start = (/ 1, 1, start_index + i /)
        call check (nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)
        accum_var = accum_var + var
      end do
 
      var = accum_var
 
      if (error /= 0) then
        deallocate (accum_var)
      end if
 
    else
 
      start = (/ 1, 1, time_index /)
      count = (/ nlons, nlats, 1 /)
      allocate (var(nlons, nlats, 1))
 
      call check (nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)
 
    end if
 
    if (error /= 0) then
      deallocate (var)
    end if
 
    deallocate (fcst_times)
 
  end if
 
  ! Close the file.
  i = nf90_close (ncid)
 
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
end subroutine read_nc_file
 
