subroutine read_nwp (currentTime,nwp_timestep_file,nPredict,n_nwp,nwp_vars,station_grid,x,z,error)
  use string_mod
  use utim
  use type
  use netcdf
  implicit none

  character (len=2000), intent (in) :: nwp_timestep_file
  character (len=100), intent (in)  :: nwp_vars(:)
  integer(i4b), intent(in)          :: nPredict     !total number of predictors
  integer(i4b), intent(in)          :: n_nwp        !number of NWP predictors
  integer(i4b), intent(in)          :: station_grid(:)        !nearest grid point for every station location
  real(dp), intent(in)          :: currentTime  !current timestep unix time
  real (dp), intent (inout) :: x(:,:), z(:,:) !station and grid predictor matrices
  integer(i4b), intent (inout) :: error !error integer

  !local variables
  integer(I4B)  :: utime     !unix time variable
  integer(I4B)  :: ncid     !
  integer(I4B)  :: varid     !
  integer (I4b) :: ngrids    !number of grid points
  integer (i4b) :: vshape (2) !shape of matricies
  integer (i4b) :: i,s          ! counter variables
  integer (i4b) :: nbase
  integer (i4b) :: nSta
  integer(I4b)  :: nlat,nlon    !lat and lon dimensions
  integer(I4b)  :: latid,lonid  !lat and lon dimension ids

  character(len=100) :: nwpTime

  real(dp), allocatable :: var(:,:)
  real(dp), allocatable :: var1d(:)

  !code starts below
  error = 0

  nBase = nPredict - n_nwp

  !the NWP file should have the same date as the current time step 
  !convert NWP valid time to unix time
  call check (nf90_open(nwp_timestep_file, nf90_nowrite, ncid), "File open error", error)
  if (error /= 0) stop

  call check (nf90_get_att(ncid,NF90_GLOBAL,"valid_time",nwpTime),"NWP file does not have valid_time global attribute",error)
  if (error /= 0) stop

  utime = date_to_unix (trim(nwpTime))

  !check time stamps of regression step and NWP file
  if(utime /= currentTime) then  !may need to give a range here for comparison...
    print *,'Current NWP time: ',utime,trim(nwpTime),' does not match regression timestep: ',currentTime
    stop
  end if

  !get grid points
  vshape = shape(z)
  ngrids = vshape(1)

  allocate(var1d(ngrids))

  !get nSta
  vshape = shape(x)
  nSta = vshape(1)

  !run through nwp_vars list

  do i = 1, n_nwp, 1
    call check (nf90_inq_varid(ncid, nwp_vars(i), varid),"Error inquiring NWP variable",error)
    if (error /= 0) stop

    !get dimensions
    call check (nf90_inq_dimid(ncid, "lat", latid), "Make sure there is a lat dimension in NWP file", error)
    if (error /= 0) stop
    call check (nf90_inquire_dimension(ncid, latid,len=nlat),"Error reading lat dimension in NWP file",error)
    if (error /= 0) stop
    call check (nf90_inq_dimid(ncid, "lon", lonid), "Make sure there is a lon dimension in NWP file", error)
    if (error /= 0) stop
    call check (nf90_inquire_dimension(ncid, lonid, len=nlon), "Error reading lon dimension in NWP file", error)
    if (error /= 0) stop

    if(i .eq. 1) then
      allocate( var(nlon,nlat) )
    end if

    call check (nf90_get_var(ncid,varid,var),"Error reading NWP variable",error)
    if (error /= 0) stop

    !reformat var and place in z matrix
    vshape = shape(var)

    if(vshape(1)*vshape(2) /= ngrids) then
      print *,"Variable points: ",vshape(1)*vshape(2)," do not match input DEM size",ngrids
      stop
    end if
 
    var1d = reshape(var, (/ ngrids /))
    z(1:ngrids,nbase+i) = var1d

    !now assign station matrix
    do s = 1, nSta, 1
      x(s,nbase+i) = var1d(station_grid(s))
    end do

  end do

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
end subroutine read_nwp

