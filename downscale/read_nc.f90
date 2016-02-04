 
SUBROUTINE read_nc_file (file_name, valid_time, var_name, var, lats, lons, error)
  USE netcdf
  USE type
  IMPLICIT NONE
 
  CHARACTER (LEN=*), INTENT (IN) :: file_name
  CHARACTER (LEN=*), INTENT (IN) :: var_name
  REAL (DP), INTENT (IN) :: valid_time
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: var (:, :, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: lats (:), lons (:)
  INTEGER, INTENT (OUT) :: error
 
  INTEGER (I4B) :: nlats, nlons, ntimes
  INTEGER (I4B) :: ncid, varid, latid, lonid, timeid
  INTEGER (I4B) :: ndims, lon_dimid, lat_dimid, time_dimid
  CHARACTER (LEN=*), PARAMETER :: LAT_NAME = "latitude"
  CHARACTER (LEN=*), PARAMETER :: LON_NAME = "longitude"
  CHARACTER (LEN=*), PARAMETER :: TIME_NAME = "time"
 
 
  ! We recommend that each variable carry a "units" attribute.
  CHARACTER (LEN=*), PARAMETER :: UNITS = "units"
  CHARACTER (LEN=*), PARAMETER :: LAT_UNITS = "degrees_north"
  CHARACTER (LEN=*), PARAMETER :: LON_UNITS = "degrees_east"
 
  INTEGER (I4B), DIMENSION (nf90_max_var_dims) :: dimIds
  REAL (DP), ALLOCATABLE :: times (:), fcst_times (:)
 
  ! Loop indices
  INTEGER (I4B) :: lvl, lat, lon, rec, i, time_index
 
  !indicies for aggregating variables between forecast times
  INTEGER (I4B) :: start_index
 
  INTEGER (I4B) :: start (3), count (3)
 
  ! To check the units attributes.
  CHARACTER * 80 pres_units_in, temp_units_in
  CHARACTER * 80 lat_units_in, lon_units_in
 
  REAL (DP), ALLOCATABLE :: accum_var (:, :, :)
 
  error = 0
  time_index = 0
 
  ! Open the file.
  CALL check (nf90_open(file_name, nf90_nowrite, ncid), "File open error", error)
  IF (error /= 0) RETURN
 
  ! Get the varids of the latitude and longitude coordinate variables.
  CALL check (nf90_inq_varid(ncid, LAT_NAME, latid), "Latitutde name error", error)
  CALL check (nf90_inquire_variable(ncid, latid, ndims=ndims, dimIds=dimIds), "Latitutde inq error", error)
  IF (error /= 0 .OR. ndims /= 1) RETURN
  lat_dimid = dimIds (1)
 
  CALL check (nf90_inq_varid(ncid, LON_NAME, lonid), "Longitude name error", error)
  CALL check (nf90_inquire_variable(ncid, lonid, ndims=ndims, dimIds=dimIds), "Longitude inq error", error)
  IF (error /= 0 .OR. ndims /= 1) RETURN
  lon_dimid = dimIds (1)
 
  CALL check (nf90_inq_varid(ncid, TIME_NAME, timeid), "Time name error", error)
  CALL check (nf90_inquire_variable(ncid, timeid, ndims=ndims, dimIds=dimIds), "Time inq error", error)
  IF (error /= 0 .OR. ndims /= 1) RETURN
 
  CALL check (nf90_inquire_dimension(ncid, lat_dimid, len=nlats), "Latitutde dim error", error)
  CALL check (nf90_inquire_dimension(ncid, lon_dimid, len=nlons), "Longitude dim error", error)
  IF (error /= 0) RETURN
 
  ALLOCATE (lats(nlats), lons(nlons))
  ! Read the latitude and longitude data.
  CALL check (nf90_get_var(ncid, latid, lats), "Latitutde read error", error)
  CALL check (nf90_get_var(ncid, lonid, lons), "Longitude read error", error)
 
 
  CALL check (nf90_inq_varid(ncid, var_name, varid), "Variable name error", error)
  CALL check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimIds=dimIds), "Variable inq error", error)
 
  IF (error == 0 .AND. ndims == 3 .AND. dimIds(2) == latid .AND. dimIds(1) == lonid) THEN
 
   time_dimid = dimIds (3)
 
   CALL check (nf90_inquire_dimension(ncid, time_dimid, len=ntimes), "Time dim error", error)
 
   ALLOCATE (times(ntimes))
   ALLOCATE (fcst_times(ntimes))
 
   CALL check (nf90_get_var(ncid, timeid, times), "Time read error", error)
   IF (error /= 0) RETURN
 
   fcst_times = times
 
   DO i = 1, ntimes, 1
    IF (times(i) == valid_time) time_index = i
   END DO
   DEALLOCATE (times)
   IF (time_index == 0) THEN
    PRINT ("(AF11.0A))"), "Error: time (", valid_time, ") not found in this file. "
    error = 1
    DEALLOCATE (lats)
    DEALLOCATE (lons)
    RETURN
   END IF
 
   IF (trim(var_name) == 'APCP_ens_mean_surface' .OR. trim(var_name) == 'DSWRF_ens_mean_surface') THEN
       !find fcst_time that is ~24 hours prior
       !just find closest one that is either 24-hrs or less
 
    DO i = 1, ntimes, 1
     IF (fcst_times(i) >= valid_time-86400) THEN
      start_index = i
      EXIT
     END IF
    END DO
 
    start = (/ 1, 1, start_index /)
    count = (/ nlons, nlats, 1 /)
    ALLOCATE (var(nlons, nlats, 1))
    ALLOCATE (accum_var(nlons, nlats, 1))
    accum_var = 0.0
 
    DO i = 1, (time_index-start_index), 1
     start = (/ 1, 1, start_index + i /)
     CALL check (nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)
     accum_var = accum_var + var
    END DO
 
    var = accum_var
 
    IF (error /= 0) THEN
     DEALLOCATE (accum_var)
    END IF
 
   ELSE
 
    start = (/ 1, 1, time_index /)
    count = (/ nlons, nlats, 1 /)
    ALLOCATE (var(nlons, nlats, 1))
 
    CALL check (nf90_get_var(ncid, varid, var, start, count), "Variable read error", error)
 
   END IF
 
   IF (error /= 0) THEN
    DEALLOCATE (var)
   END IF
 
   DEALLOCATE (fcst_times)
 
  END IF
 
  ! Close the file.
  i = nf90_close (ncid)
 
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
END SUBROUTINE read_nc_file
 
