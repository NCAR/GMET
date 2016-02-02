subroutine save_vars(pcp, tmean, trange, nx, ny, grdlat, grdlon, grdalt, Times, file, error)
  use netcdf
  use nrtype
  implicit none

  real(SP), intent(in) :: pcp(:,:,:), tmean(:,:,:), trange(:,:,:)
  integer(I4B), intent(in) :: nx, ny
  real(DP), intent(in) :: grdlat(:), grdlon(:), grdalt(:)
  real(DP), intent(in) :: Times(:)
  character (len = 500), intent(in) :: file
  integer, intent(out) :: error


  ! Dimension names
  character (len = *), parameter :: Y_NAME = "y"
  character (len = *), parameter :: X_NAME = "x"
  character (len = *), parameter :: TIME_NAME = "time"

  ! Variable Names
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: ALT_NAME = "elevation"
  character (len = *), parameter :: PCP_NAME = "pcp"
  character (len = *), parameter :: POP_NAME = "t_mean"
  character (len = *), parameter :: PCP_ERROR_NAME = "t_range"

  character (len = *), parameter :: LONG_NAME = "long_name"
  character (len = *), parameter :: PCP_LONG_NAME = "estimated precip in mm/day"
  character (len = *), parameter :: POP_LONG_NAME = "estimated daily mean temperature"
  character (len = *), parameter :: PCP_ERROR_LONG_NAME = "estimated diurnal range"

  ! Units
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PCP_UNITS = "mm"
  character (len = *), parameter :: POP_UNITS = "deg_C"
  character (len = *), parameter :: PCP_ERROR_UNITS = "deg_C"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: ALT_UNITS = "meters"
  character (len = *), parameter :: TIME_UNITS = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len = *), parameter :: FILL = "_FillValue"

  real(DP), allocatable :: file_times(:)

  integer :: n_chars, n_times, inx, iny
  integer :: ncid, x_dimid, y_dimid, time_dimid
  integer :: lat_varid, lon_varid, alt_varid, time_varid, pcp_varid, pop_varid, pcp_error_varid
  integer :: count1(1), start1(1), count2(2), start2(2), count3(3), start3(3), dimids2(2), dimids3(3)
  integer :: trec, nrecs, file_nx, file_ny, file_ntimes, i

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
     if(error /=0 ) return

     dimids3 = (/ x_dimid, y_dimid, time_dimid /)
     call check( nf90_def_var(ncid, PCP_NAME, NF90_FLOAT, dimids3, pcp_varid), &
          "pcp var def error", error)
     call check( nf90_def_var(ncid, POP_NAME, NF90_FLOAT, dimids3, pop_varid), &
          "pop var def error", error)
     call check( nf90_def_var(ncid, PCP_ERROR_NAME, NF90_FLOAT, dimids3, pcp_error_varid), &
          "pcp_error var def error", error)
     if(error /=0 ) return

    ! Add attributes. 
     call check( nf90_put_att(ncid, pcp_varid, LONG_NAME, PCP_LONG_NAME), "pcp long_name attribute error", error)
     call check( nf90_put_att(ncid, pop_varid, LONG_NAME, POP_LONG_NAME), "pcp long_name attribute error", error)
     call check( nf90_put_att(ncid, pcp_error_varid, LONG_NAME, PCP_ERROR_LONG_NAME), "pcp_error long_name attribute error", error)

     call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS), "lat units attribute error", error)
     call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS), "lon units attribute error", error)
     call check( nf90_put_att(ncid, alt_varid, UNITS, ALT_UNITS), "alt units attribute error", error)
     call check( nf90_put_att(ncid, time_varid, UNITS, TIME_UNITS), "time units attribute error", error)

     call check( nf90_put_att(ncid, pcp_varid, UNITS, PCP_UNITS), "pcp units attribute error", error)
     call check( nf90_put_att(ncid, pop_varid, UNITS, POP_UNITS), "pcp units attribute error", error)
     call check( nf90_put_att(ncid, pcp_error_varid, UNITS, PCP_ERROR_UNITS), "pcp_error units attribute error", error)
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
     call check( nf90_inq_dimid(ncid, Y_NAME, y_dimid), "y dim inq error", error)
     call check( nf90_inq_dimid(ncid, X_NAME, x_dimid), "x dim inq error", error)
     call check( nf90_inq_dimid(ncid, TIME_NAME, time_dimid), "time dim inq error", error)
     if(error /=0 ) return

     call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid), "lat var inq error", error)
     call check( nf90_inq_varid(ncid, LON_NAME, lon_varid), "lon var inq error", error)
     call check( nf90_inq_varid(ncid, TIME_NAME, time_varid), "time var inq error", error)
     call check( nf90_inq_varid(ncid, PCP_NAME, pcp_varid), "pcp var inq error", error)
     call check( nf90_inq_varid(ncid, POP_NAME, pop_varid), "pop var inq error", error)
     call check( nf90_inq_varid(ncid, PCP_ERROR_NAME, pcp_error_varid), "pcp_error var inq error", error)
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

  count3 = (/ inx, iny, n_times /)
  start3 = (/ 1, 1, trec /)
  call check( nf90_put_var(ncid, pcp_varid, pcp, start = start3, count = count3), &
       "put pcp error", error)
  if(error /=0 ) return

  call check( nf90_put_var(ncid, pop_varid, tmean, start = start3, count = count3), &
       "put tmean error", error)
  if(error /=0 ) return
  call check( nf90_put_var(ncid, pcp_error_varid, trange, start = start3, count = count3), &
       "put trange error", error)
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
end subroutine save_vars