Subroutine save_precip (pcp, pop, pcperror, nx, ny, grdlat, grdlon, grdalt, Times, file, error)
  Use netcdf
  Use nrtype
  Implicit None
 
  Real (DP), Intent (In) :: pcp (:, :), pop (:, :), pcperror (:, :)
  Integer (I4B), Intent (In) :: nx, ny
  Real (DP), Intent (In) :: grdlat (:), grdlon (:), grdalt (:)
  Real (DP), Intent (In) :: Times (:)
  Character (Len=500), Intent (In) :: file
  Integer, Intent (Out) :: error
 
 
  ! Dimension names
  Character (Len=*), Parameter :: Y_NAME = "y"
  Character (Len=*), Parameter :: X_NAME = "x"
  Character (Len=*), Parameter :: TIME_NAME = "time"
 
  ! Variable Names
  Character (Len=*), Parameter :: LAT_NAME = "latitude"
  Character (Len=*), Parameter :: LON_NAME = "longitude"
  Character (Len=*), Parameter :: ALT_NAME = "altitude"
  Character (Len=*), Parameter :: PCP_NAME = "pcp"
  Character (Len=*), Parameter :: POP_NAME = "pop"
  Character (Len=*), Parameter :: PCP_ERROR_NAME = "pcp_error"
 
  Character (Len=*), Parameter :: LONG_NAME = "long_name"
  Character (Len=*), Parameter :: PCP_LONG_NAME = "estimated precip in normal space"
  Character (Len=*), Parameter :: POP_LONG_NAME = "probability of precipitation occurrence"
  Character (Len=*), Parameter :: PCP_ERROR_LONG_NAME = "error in estimated precip"
 
  ! Units
  Character (Len=*), Parameter :: UNITS = "units"
  Character (Len=*), Parameter :: PCP_UNITS = ""
  Character (Len=*), Parameter :: POP_UNITS = ""
  Character (Len=*), Parameter :: PCP_ERROR_UNITS = ""
  Character (Len=*), Parameter :: LAT_UNITS = "degrees_north"
  Character (Len=*), Parameter :: LON_UNITS = "degrees_east"
  Character (Len=*), Parameter :: ALT_UNITS = "meters"
  Character (Len=*), Parameter :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  Character (Len=*), Parameter :: FILL = "_FillValue"
 
  Real (DP), Allocatable :: file_times (:)
 
  Integer :: n_chars, n_times, inx, iny
  Integer :: ncid, x_dimid, y_dimid, time_dimid
  Integer :: lat_varid, lon_varid, alt_varid, time_varid, pcp_varid, pop_varid, pcp_error_varid
  Integer :: count1 (1), start1 (1), count2 (2), start2 (2), count3 (3), start3 (3), dimids2 (2), dimids3 (3)
  Integer :: trec, nrecs, file_nx, file_ny, file_ntimes, i
 
  trec = 0
  n_chars = 100
  n_times = size (Times)
  inx = nx
  iny = ny
 
  If (size(grdlat) /= inx*iny) Then
    Print *, "Error "
  End If
 
  error = nf90_open (file, nf90_write, ncid)
  If (error /= nf90_noerr) Then
    error = 0
     ! Create the file.
    Call check (nf90_create(file, nf90_clobber, ncid), "File creation error", error)
    If (error /= 0) Return
 
     ! Define the dimensions.
    Call check (nf90_def_dim(ncid, X_NAME, inx, x_dimid), "x dim def error", error)
    Call check (nf90_def_dim(ncid, Y_NAME, iny, y_dimid), "y dim def error", error)
    Call check (nf90_def_dim(ncid, TIME_NAME, NF90_UNLIMITED, time_dimid), "time dim def error", error)
    If (error /= 0) Return
 
     ! Define the variables.
    dimids2 = (/ y_dimid, x_dimid /)
    Call check (nf90_def_var(ncid, LAT_NAME, NF90_DOUBLE, dimids2, lat_varid), "lat var def error", error)
    Call check (nf90_def_var(ncid, LON_NAME, NF90_DOUBLE, dimids2, lon_varid), "lon var def error", error)
    Call check (nf90_def_var(ncid, ALT_NAME, NF90_DOUBLE, dimids2, alt_varid), "lon var def error", error)
    Call check (nf90_def_var(ncid, TIME_NAME, NF90_DOUBLE, time_dimid, time_varid), "time var def error", error)
    If (error /= 0) Return
 
    dimids3 = (/ y_dimid, x_dimid, time_dimid /)
    Call check (nf90_def_var(ncid, PCP_NAME, NF90_DOUBLE, dimids3, pcp_varid), "pcp var def error", error)
    Call check (nf90_def_var(ncid, POP_NAME, NF90_DOUBLE, dimids3, pop_varid), "pop var def error", error)
    Call check (nf90_def_var(ncid, PCP_ERROR_NAME, NF90_DOUBLE, dimids3, pcp_error_varid), "pcp_error var def error", &
   & error)
    If (error /= 0) Return
 
    ! Add attributes.
    Call check (nf90_put_att(ncid, pcp_varid, LONG_NAME, PCP_LONG_NAME), "pcp long_name attribute error", error)
    Call check (nf90_put_att(ncid, pop_varid, LONG_NAME, POP_LONG_NAME), "pcp long_name attribute error", error)
    Call check (nf90_put_att(ncid, pcp_error_varid, LONG_NAME, PCP_ERROR_LONG_NAME), "pcp_error long_name attribute err&
   &or", error)
 
    Call check (nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS), "lat units attribute error", error)
    Call check (nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS), "lon units attribute error", error)
    Call check (nf90_put_att(ncid, alt_varid, UNITS, ALT_UNITS), "alt units attribute error", error)
    Call check (nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)
 
    Call check (nf90_put_att(ncid, pcp_varid, UNITS, PCP_UNITS), "pcp units attribute error", error)
    Call check (nf90_put_att(ncid, pop_varid, UNITS, POP_UNITS), "pcp units attribute error", error)
    Call check (nf90_put_att(ncid, pcp_error_varid, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)
    If (error /= 0) Return
 
     ! End define mode.
    Call check (nf90_enddef(ncid), "end define mode error", error)
    If (error /= 0) Return
 
    count2 = (/ iny, inx /)
    start2 = (/ 1, 1 /)
 
    Call check (nf90_put_var(ncid, lat_varid, grdlat, start=start2, count=count2), "put lat error", error)
    Call check (nf90_put_var(ncid, lon_varid, grdlon, start=start2, count=count2), "put lon error", error)
    Call check (nf90_put_var(ncid, alt_varid, grdalt, start=start2, count=count2), "put alt error", error)
 
    trec = 1
    nrecs = 0
 
 
  Else
 
     ! File already exists, get dim and var ids
    Call check (nf90_inq_dimid(ncid, X_NAME, x_dimid), "x dim inq error", error)
    Call check (nf90_inq_dimid(ncid, Y_NAME, y_dimid), "y dim inq error", error)
    Call check (nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
    If (error /= 0) Return
 
    Call check (nf90_inq_varid(ncid, LAT_NAME, lat_varid), "lat var inq error", error)
    Call check (nf90_inq_varid(ncid, LON_NAME, lon_varid), "lon var inq error", error)
    Call check (nf90_inq_varid(ncid, TIME_NAME, time_varid), "time var inq error", error)
    Call check (nf90_inq_varid(ncid, PCP_NAME, pcp_varid), "pcp var inq error", error)
    Call check (nf90_inq_varid(ncid, POP_NAME, pop_varid), "pop var inq error", error)
    Call check (nf90_inq_varid(ncid, PCP_ERROR_NAME, pcp_error_varid), "pcp_error var inq error", error)
    If (error /= 0) Return
 
    Call check (nf90_Inquire_Dimension(ncid, x_dimid, len=file_nx), "x dim len error", error)
    Call check (nf90_Inquire_Dimension(ncid, y_dimid, len=file_ny), "y dim len error", error)
    Call check (nf90_Inquire_Dimension(ncid, time_dimid, len=file_ntimes), "time dim len error", error)
    If (error /= 0) Return
 
    If (nx /= file_nx .Or. ny /= file_ny) Then
      Print *, "Error dimensions in output file do not match current run."
      error = 1
      Return
    End If
 
    Allocate (file_times(file_ntimes))
    Call check (nf90_get_var(ncid, time_varid, file_times), "error getting file times list", error)
    If (error /= 0) Return
 
    If (file_times(1) > Times(n_times)) Then !put data before everything in the file
      Print *, "Error cannot add data before data already in output file. (functionality still to be added)"
      error = 1
      Return
    Else
      If (file_times(file_ntimes) < Times(1)) Then !put data after everything in the file
        trec = file_ntimes + 1
      Else ! at least some overlap
        Do i = 1, file_ntimes, 1
          If (file_times(1) == Times(1)) Then
            trec = i
          End If
        End Do
        If (trec == 0) Then
          Print *, "Error, confusion over data output record location."
          error = 1
          Return
        Else
          Print *, "WARNING, overwriting data in output file, record ", trec, " to ", trec + n_times - 1
        End If
      End If
    End If
 
  End If
 
  count1 (1) = n_times
  start1 (1) = trec
  Call check (nf90_put_var(ncid, time_varid, Times, start=start1, count=count1), "put times error", error)
  If (error /= 0) Return
 
  count3 = (/ iny, inx, n_times /)
  start3 = (/ 1, 1, trec /)
  Call check (nf90_put_var(ncid, pcp_varid, pcp, start=start3, count=count3), "put pcp error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, pop_varid, pop, start=start3, count=count3), "put pop error", error)
  If (error /= 0) Return
 
  Call check (nf90_put_var(ncid, pcp_error_varid, pcperror, start=start3, count=count3), "put pcp_error error", error)
  If (error /= 0) Return
 
  Call check (nf90_close(ncid), "closing file error", error)
 
 
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
End Subroutine save_precip
