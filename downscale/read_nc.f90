 
Subroutine read_nc_file (file_name, valid_time, var_name, var, lats, lons, error)
  Use netcdf
  Use type
  Implicit None
 
  Character (Len=*), Intent (In) :: file_name
  Character (Len=*), Intent (In) :: var_name
  Real (DP), Intent (In) :: valid_time
  Real (DP), Allocatable, Intent (Out) :: var (:, :, :)
  Real (DP), Allocatable, Intent (Out) :: lats (:), lons (:)
  Integer, Intent (Out) :: error
 
  Integer (I4B) :: nlats, nlons, ntimes
  Integer (I4B) :: ncid, varid, latid, lonid, timeid
  Integer (I4B) :: ndims, lon_dimid, lat_dimid, time_dimid
  Character (Len=*), Parameter :: LAT_NAME = "latitude"
  Character (Len=*), Parameter :: LON_NAME = "longitude"
  Character (Len=*), Parameter :: TIME_NAME = "time"
 
 
  ! We recommend that each variable carry a "units" attribute.
  Character (Len=*), Parameter :: UNITS = "units"
  Character (Len=*), Parameter :: LAT_UNITS = "degrees_north"
  Character (Len=*), Parameter :: LON_UNITS = "degrees_east"
 
  Integer (I4B), Dimension (nf90_max_var_dims) :: dimIds
  Real (DP), Allocatable :: times (:), fcst_times (:)
 
  ! Loop indices
  Integer (I4B) :: lvl, lat, lon, rec, i, time_index
 
  !indicies for aggregating variables between forecast times
  Integer (I4B) :: start_index
 
  Integer (I4B) :: start (3), count (3)
 
  ! To check the units attributes.
  Character * 80 pres_units_in, temp_units_in
  Character * 80 lat_units_in, lon_units_in
 
  Real (DP), Allocatable :: accum_var (:, :, :)
 
  error = 0
  time_index = 0
 
  ! Open the file.
  Call check (nf90_open(file_name, nf90_nowrite, ncid), "File open error", error)
  If (error /= 0) Return
 
  ! Get the varids of the latitude and longitude coordinate variables.
  Call check (nf90_inq_varid(ncid, LAT_NAME, latid), "Latitutde name error", error)
  Call check (nf90_inquire_variable(ncid, latid, ndims=ndims, dimIds=dimIds), "Latitutde inq error", error)
  If (error /= 0 .Or. ndims /= 1) Return
  lat_dimid = dimIds (1)
 
  Call check (nf90_inq_varid(ncid, LON_NAME, lonid), "Longitude name error", error)
  Call check (nf90_inquire_variable(ncid, lonid, ndims=ndims, dimIds=dimIds), "Longitude inq error", error)
  If (error /= 0 .Or. ndims /= 1) Return
  lon_dimid = dimIds (1)
 
  Call check (nf90_inq_varid(ncid, TIME_NAME, timeid), "Time name error", error)
  Call check (nf90_inquire_variable(ncid, timeid, ndims=ndims, dimIds=dimIds), "Time inq error", error)
  If (error /= 0 .Or. ndims /= 1) Return
 
  Call check (nf90_inquire_dimension(ncid, lat_dimid, len=nlats), "Latitutde dim error", error)
  Call check (nf90_inquire_dimension(ncid, lon_dimid, len=nlons), "Longitude dim error", error)
  If (error /= 0) Return
 
  Allocate (lats(nlats), lons(nlons))
  ! Read the latitude and longitude data.
  Call check (nf90_get_var(ncid, latid, lats), "Latitutde read error", error)
  Call check (nf90_get_var(ncid, lonid, lons), "Longitude read error", error)
 
 
  Call check (nf90_inq_varid(ncid, var_name, varid), "Variable name error", error)
  Call check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimIds=dimIds), "Variable inq error", error)
 
  If (error == 0 .And. ndims == 3 .And. dimIds(2) == latid .And. dimIds(1) == lonid) Then
 
    time_dimid = dimIds (3)
 
    Call check (nf90_inquire_dimension(ncid, time_dimid, len=ntimes), "Time dim error", error)
 
    Allocate (times(ntimes))
    Allocate (fcst_times(ntimes))
 
    Call check (nf90_get_var(ncid, timeid, times), "Time read error", error)
    If (error /= 0) Return
 
    fcst_times = times
 
    Do i = 1, ntimes, 1
      If (times(i) == valid_time) time_index = i
    End Do
    Deallocate (times)
    If (time_index == 0) Then
      Print ("(AF11.0A))"), "Error: time (", valid_time, ") not found in this file. "
      error = 1
      Deallocate (lats)
      Deallocate (lons)
      Return
    End If
 
    If (trim(var_name) == 'APCP_ens_mean_surface' .Or. trim(var_name) == 'DSWRF_ens_mean_surface') Then
       !find fcst_time that is ~24 hours prior
       !just find closest one that is either 24-hrs or less
 
      Do i = 1, ntimes, 1
        If (fcst_times(i) >= valid_time-86400) Then
          start_index = i
          Exit
        End If
      End Do
 
      start = (/ 1, 1, start_index /)
      count = (/ nlons, nlats, 1 /)
      Allocate (var(nlons, nlats, 1))
      Allocate (accum_var(nlons, nlats, 1))
      accum_var = 0.0
 
      Do i = 1, (time_index-start_index), 1
        start = (/ 1, 1, start_index + i /)
        Call check (nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)
        accum_var = accum_var + var
      End Do
 
      var = accum_var
 
      If (error /= 0) Then
        Deallocate (accum_var)
      End If
 
    Else
 
      start = (/ 1, 1, time_index /)
      count = (/ nlons, nlats, 1 /)
      Allocate (var(nlons, nlats, 1))
 
      Call check (nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)
 
    End If
 
    If (error /= 0) Then
      Deallocate (var)
    End If
 
    Deallocate (fcst_times)
 
  End If
 
  ! Close the file.
  i = nf90_close (ncid)
 
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
End Subroutine read_nc_file
 
