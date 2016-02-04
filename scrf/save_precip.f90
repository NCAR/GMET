SUBROUTINE save_precip (pcp, pop, pcperror, nx, ny, grdlat, grdlon, grdalt, Times, file, error)
  USE netcdf
  USE nrtype
  IMPLICIT NONE
 
  REAL (DP), INTENT (IN) :: pcp (:, :), pop (:, :), pcperror (:, :)
  INTEGER (I4B), INTENT (IN) :: nx, ny
  REAL (DP), INTENT (IN) :: grdlat (:), grdlon (:), grdalt (:)
  REAL (DP), INTENT (IN) :: Times (:)
  CHARACTER (LEN=500), INTENT (IN) :: file
  INTEGER, INTENT (OUT) :: error
 
 
  ! Dimension names
  CHARACTER (LEN=*), PARAMETER :: Y_NAME = "y"
  CHARACTER (LEN=*), PARAMETER :: X_NAME = "x"
  CHARACTER (LEN=*), PARAMETER :: TIME_NAME = "time"
 
  ! Variable Names
  CHARACTER (LEN=*), PARAMETER :: LAT_NAME = "latitude"
  CHARACTER (LEN=*), PARAMETER :: LON_NAME = "longitude"
  CHARACTER (LEN=*), PARAMETER :: ALT_NAME = "altitude"
  CHARACTER (LEN=*), PARAMETER :: PCP_NAME = "pcp"
  CHARACTER (LEN=*), PARAMETER :: POP_NAME = "pop"
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_NAME = "pcp_error"
 
  CHARACTER (LEN=*), PARAMETER :: LONG_NAME = "long_name"
  CHARACTER (LEN=*), PARAMETER :: PCP_LONG_NAME = "estimated precip in normal space"
  CHARACTER (LEN=*), PARAMETER :: POP_LONG_NAME = "probability of precipitation occurrence"
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_LONG_NAME = "error in estimated precip"
 
  ! Units
  CHARACTER (LEN=*), PARAMETER :: UNITS = "units"
  CHARACTER (LEN=*), PARAMETER :: PCP_UNITS = ""
  CHARACTER (LEN=*), PARAMETER :: POP_UNITS = ""
  CHARACTER (LEN=*), PARAMETER :: PCP_ERROR_UNITS = ""
  CHARACTER (LEN=*), PARAMETER :: LAT_UNITS = "degrees_north"
  CHARACTER (LEN=*), PARAMETER :: LON_UNITS = "degrees_east"
  CHARACTER (LEN=*), PARAMETER :: ALT_UNITS = "meters"
  CHARACTER (LEN=*), PARAMETER :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  CHARACTER (LEN=*), PARAMETER :: FILL = "_FillValue"
 
  REAL (DP), ALLOCATABLE :: file_times (:)
 
  INTEGER :: n_chars, n_times, inx, iny
  INTEGER :: ncid, x_dimid, y_dimid, time_dimid
  INTEGER :: lat_varid, lon_varid, alt_varid, time_varid, pcp_varid, pop_varid, pcp_error_varid
  INTEGER :: count1 (1), start1 (1), count2 (2), start2 (2), count3 (3), start3 (3), dimids2 (2), dimids3 (3)
  INTEGER :: trec, nrecs, file_nx, file_ny, file_ntimes, i
 
  trec = 0
  n_chars = 100
  n_times = size (Times)
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
   CALL check (nf90_def_dim(ncid, X_NAME, inx, x_dimid), "x dim def error", error)
   CALL check (nf90_def_dim(ncid, Y_NAME, iny, y_dimid), "y dim def error", error)
   CALL check (nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), "time dim def error", error)
   IF (error /= 0) RETURN
 
     ! Define the variables.
   dimids2 = (/ y_dimid, x_dimid /)
   CALL check (nf90_def_var(ncid, LAT_NAME, NF90_DOUBLE, dimids2, lat_varid), "lat var def error", error)
   CALL check (nf90_def_var(ncid, LON_NAME, NF90_DOUBLE, dimids2, lon_varid), "lon var def error", error)
   CALL check (nf90_def_var(ncid, ALT_NAME, NF90_DOUBLE, dimids2, alt_varid), "lon var def error", error)
   CALL check (nf90_def_var(ncid, TIME_NAME, NF90_DOUBLE, time_dimid, time_varid), "time var def error", error)
   IF (error /= 0) RETURN
 
   dimids3 = (/ y_dimid, x_dimid, time_dimid /)
   CALL check (nf90_def_var(ncid, PCP_NAME, NF90_DOUBLE, dimids3, pcp_varid), "pcp var def error", error)
   CALL check (nf90_def_var(ncid, POP_NAME, NF90_DOUBLE, dimids3, pop_varid), "pop var def error", error)
   CALL check (nf90_def_var(ncid, PCP_ERROR_NAME, NF90_DOUBLE, dimids3, pcp_error_varid), "pcp_error var def error", &
  & error)
   IF (error /= 0) RETURN
 
    ! Add attributes.
   CALL check (nf90_put_att(ncid, pcp_varid, LONG_NAME, PCP_LONG_NAME), "pcp long_name attribute error", error)
   CALL check (nf90_put_att(ncid, pop_varid, LONG_NAME, POP_LONG_NAME), "pcp long_name attribute error", error)
   CALL check (nf90_put_att(ncid, pcp_error_varid, LONG_NAME, PCP_ERROR_LONG_NAME), "pcp_error long_name attribute erro&
  &r", error)
 
   CALL check (nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS), "lat units attribute error", error)
   CALL check (nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS), "lon units attribute error", error)
   CALL check (nf90_put_att(ncid, alt_varid, UNITS, ALT_UNITS), "alt units attribute error", error)
   CALL check (nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)
 
   CALL check (nf90_put_att(ncid, pcp_varid, UNITS, PCP_UNITS), "pcp units attribute error", error)
   CALL check (nf90_put_att(ncid, pop_varid, UNITS, POP_UNITS), "pcp units attribute error", error)
   CALL check (nf90_put_att(ncid, pcp_error_varid, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)
   IF (error /= 0) RETURN
 
     ! End define mode.
   CALL check (nf90_enddef(ncid), "end define mode error", error)
   IF (error /= 0) RETURN
 
   count2 = (/ iny, inx /)
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
 
   IF (file_times(1) > Times(n_times)) THEN !put data before everything in the file
    PRINT *, "Error cannot add data before data already in output file. (functionality still to be added)"
    error = 1
    RETURN
   ELSE
    IF (file_times(file_ntimes) < Times(1)) THEN !put data after everything in the file
     trec = file_ntimes + 1
    ELSE ! at least some overlap
     DO i = 1, file_ntimes, 1
      IF (file_times(1) == Times(1)) THEN
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
  CALL check (nf90_put_var(ncid, time_varid, Times, start=start1, count=count1), "put times error", error)
  IF (error /= 0) RETURN
 
  count3 = (/ iny, inx, n_times /)
  start3 = (/ 1, 1, trec /)
  CALL check (nf90_put_var(ncid, pcp_varid, pcp, start=start3, count=count3), "put pcp error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, pop_varid, pop, start=start3, count=count3), "put pop error", error)
  IF (error /= 0) RETURN
 
  CALL check (nf90_put_var(ncid, pcp_error_varid, pcperror, start=start3, count=count3), "put pcp_error error", error)
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
