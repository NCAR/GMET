
subroutine read_nc_file(file_name, valid_time, var_name, var, lats, lons, error)
  use netcdf
  use type
  implicit none
  
  character (len = *), intent(in) :: file_name
  character (len = *), intent(in) :: var_name
  real(DP), intent(in) :: valid_time
  real(DP), allocatable, intent(out) :: var(:, :, :)
  real(DP), allocatable, intent(out) :: lats(:), lons(:)
  integer, intent(out) :: error

  integer(I4B) :: nlats, nlons, ntimes
  integer(I4B) :: ncid, varid, latid, lonid, timeid
  integer(I4B) :: ndims, lon_dimid, lat_dimid, time_dimid
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: TIME_NAME = "time"


  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"

  integer(I4B), dimension(nf90_max_var_dims) :: dimIds
  real(DP), allocatable :: times(:), fcst_times(:) 

  ! Loop indices
  integer(I4B) :: lvl, lat, lon, rec, i, time_index

  !indicies for aggregating variables between forecast times
  integer(I4B) :: start_index

  integer(I4B) :: start(3), count(3)

  ! To check the units attributes.
  character*80 pres_units_in, temp_units_in
  character*80 lat_units_in, lon_units_in

  real(DP), allocatable :: accum_var(:,:,:)

  error = 0
  time_index = 0

  ! Open the file. 
  call check( nf90_open(file_name, nf90_nowrite, ncid), "File open error", error)
  if(error /= 0) return

  ! Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, LAT_NAME, latid), "Latitutde name error", error)
  call check ( nf90_inquire_variable(ncid, latid, ndims = ndims, dimids = dimIds), "Latitutde inq error", error)
  if(error /= 0 .OR. ndims /= 1) return
  lat_dimid = dimIds(1)
  
  call check( nf90_inq_varid(ncid, LON_NAME, lonid), "Longitude name error", error)
  call check ( nf90_inquire_variable(ncid, lonid, ndims = ndims, dimids = dimIds), "Longitude inq error", error)
  if(error /= 0 .OR. ndims /= 1) return
  lon_dimid = dimIds(1)
  
  call check( nf90_inq_varid(ncid, TIME_NAME, timeid), "Time name error", error)
  call check ( nf90_inquire_variable(ncid, timeid, ndims = ndims, dimids = dimIds), "Time inq error", error)
  if(error /= 0 .OR. ndims /= 1) return

  call check ( nf90_inquire_dimension(ncid, lat_dimid, len = nlats), "Latitutde dim error", error)
  call check ( nf90_inquire_dimension(ncid, lon_dimid, len = nlons), "Longitude dim error", error)
  if(error /=0 ) return
  
  allocate(lats(nlats), lons(nlons))
  ! Read the latitude and longitude data.
  call check( nf90_get_var(ncid, latid, lats), "Latitutde read error", error)
  call check( nf90_get_var(ncid, lonid, lons), "Longitude read error", error)

  
  call check( nf90_inq_varid(ncid, var_name, varid), "Variable name error", error)
  call check ( nf90_inquire_variable(ncid, varid, ndims = ndims, dimids = dimIds), "Variable inq error", error)

  if(error == 0 .AND. ndims == 3 .AND. dimIds(2) == latid .AND. dimIds(1) == lonid) then

     time_dimid = dimIds(3)
     
     call check ( nf90_inquire_dimension(ncid, time_dimid, len = ntimes), "Time dim error", error)

     allocate(times(ntimes))
     allocate(fcst_times(ntimes))
     
     call check( nf90_get_var(ncid, timeid, times), "Time read error", error)
     if(error /= 0) return

     fcst_times = times

     do i = 1, ntimes, 1
        if(times(i) == valid_time) time_index = i
     enddo
     deallocate(times)
     if(time_index == 0) then
        print ("(AF11.0A))"), "Error: time (",valid_time, &
             ") not found in this file. "
        error = 1
        deallocate(lats)
        deallocate(lons)
        return
     endif

     if(trim(var_name) == 'APCP_ens_mean_surface' .or. trim(var_name) == 'DSWRF_ens_mean_surface') then
       !find fcst_time that is ~24 hours prior
       !just find closest one that is either 24-hrs or less
      
       do i = 1,ntimes,1
	 if(fcst_times(i) >= valid_time-86400) then
           start_index = i
	   exit
	 endif
       enddo

       start = (/ 1, 1, start_index /)
       count = (/ nlons, nlats,  1 /)
       allocate(var(nlons, nlats, 1))
       allocate(accum_var(nlons,nlats,1))
       accum_var = 0.0

       do i = 1,(time_index-start_index),1
	 start = (/ 1, 1, start_index+i /)
	 call check( nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)
	 accum_var = accum_var + var
       enddo

       var = accum_var

       if(error /= 0) then
         deallocate(accum_var)
       endif

     else

       start = (/ 1, 1, time_index /)
       count = (/ nlons, nlats,  1 /)
       allocate(var(nlons, nlats, 1))

       call check( nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)

     endif

     if(error /= 0) then
        deallocate(var)
     endif

     deallocate(fcst_times)

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
end subroutine read_nc_file

