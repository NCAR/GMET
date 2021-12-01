subroutine read_domain_grid (grid_list, lat, lon, elev, grad_n_s, grad_w_e, mask, nlon, nlat, error)
  use netcdf
  use type
  implicit none
 
  character (len=500),    intent (in)  :: grid_list
 
  real (dp), allocatable, intent (out) :: lat (:, :)
  real (dp), allocatable, intent (out) :: lon (:, :)
  real (dp), allocatable, intent (out) :: elev (:, :)
  real (dp), allocatable, intent (out) :: grad_n_s (:, :)
  real (dp), allocatable, intent (out) :: grad_w_e (:, :)
  real (dp), allocatable, intent (out) :: mask (:, :)
  integer (i4b),          intent (out) :: nlon
  integer (i4b),          intent (out) :: nlat
  integer,                intent (out) :: error
 
  ! local variables
  integer :: ncid                     ! netcdf file id
  integer :: i
  integer :: varid                    ! variable id
  !integer (i4b), dimension (nf90_max_var_dims) :: dim_ids ! dimension ids for dimensions of grid file
  integer, dimension (2) :: dim_ids   ! dimension ids for dimensions of grid
  !integer :: n_dims, nlat, nlon
  integer :: n_dims
 
  character (len=*), parameter :: lat_name = "latitude"
  character (len=*), parameter :: lon_name = "longitude"
  character (len=*), parameter :: elev_name = "elev"
  character (len=*), parameter :: grad_ns_name = "gradient_n_s"
  character (len=*), parameter :: grad_we_name = "gradient_w_e"
  character (len=*), parameter :: mask_name = "mask"
 
  ! ======= code starts below ========
 
  ! open netcdf grid file
  print*, 'reading domain file ', grid_list
  call check (nf90_open(grid_list, nf90_nowrite, ncid), "File open error", error)
  if (error /= 0) then 
    print*, 'Could not find domain grid -- check path/name in config file'
    stop
  endif
 
  ! inquire about latitude
  call check (nf90_inq_varid(ncid, lat_name, varid), "Latitude name error", error)
  print*, "Latitude varid=", varid
  if (error /= 0) then
     print*,  "Latitude name error in read_domain_grid()"
     stop
  endif

  ! get number of dimensions
  !call check (nf90_inquire_variable(ncid, varid, ndims=n_dims, dimids=dim_ids), "Dim inq error", error)
  call check (nf90_inquire_variable(ncid, varid, ndims=n_dims), "Dim inq error:", error)
  if (error /= 0 .or. n_dims /= 2) then
     print*,  "Variable inq error for latitude variable in read_domain_grid()"
     stop
  endif
 
  ! get x,y dimensions
  call check (nf90_inq_dimid(ncid, "x", dim_ids(1)), "x dimid error", error)
  call check (nf90_inq_dimid(ncid, "y", dim_ids(2)), "y dimid error", error)
  print*, "x dimid = ",dim_ids(1)
  print*, "y dimid = ",dim_ids(2)

  call check (nf90_inquire_dimension(ncid, dim_ids(1), len=nlon), "x dim error:", error)
  call check (nf90_inquire_dimension(ncid, dim_ids(2), len=nlat), "y dim error:", error)
  if (error /= 0) then
     print*,  "x or y dimension inq error in read_domain_grid()"
     stop
  endif
 
  print*, "nlat=", nlat, "nlon=", nlon
 
  ! allocate output variables
  allocate (lat(nlon, nlat))
  allocate (lon(nlon, nlat))
  allocate (elev(nlon, nlat))
  allocate (grad_n_s(nlon, nlat))
  allocate (grad_w_e(nlon, nlat))
  allocate (mask(nlon, nlat))
  !print*, "shape of lat var is:", shape(lat)
 
  ! --- read lat, lon, elev, grad_n_s, grad_w_e ---
 
  ! get latitude
  call check (nf90_get_var(ncid, varid, lat), "Latitude read error", error)
  if (error /= 0) then
    print*, error
    stop
  end if
 
  ! get longitude
  call check (nf90_inq_varid(ncid, lon_name, varid), "Longitude name error", error)
  if (error /= 0) return
  call check (nf90_get_var(ncid, varid, lon), "Longitude read error", error)
  if (error /= 0) return
 
  ! get elevation
  call check (nf90_inq_varid(ncid, elev_name, varid), "Elevation name error", error)
  if (error /= 0) return
  call check (nf90_get_var(ncid, varid, elev), "Elevation read error", error)
  if (error /= 0) return
 
  ! get n_s gradient
  call check (nf90_inq_varid(ncid, grad_ns_name, varid), "Grad_N_S name error", error)
  if (error /= 0) return
  call check (nf90_get_var(ncid, varid, grad_n_s), "Grad_N_S read error", error)
  if (error /= 0) return
 
  ! get w_e gradient
  call check (nf90_inq_varid(ncid, grad_we_name, varid), "Grad_W_E name error", error)
  if (error /= 0) return
  call check (nf90_get_var(ncid, varid, grad_w_e), "Grad_W_E read error", error)
  if (error /= 0) return
 
  ! get basin mask
  call check (nf90_inq_varid(ncid, mask_name, varid), "Mask name error", error)
  if (error /= 0) return
  call check (nf90_get_var(ncid, varid, mask), "Mask read error", error)
  if (error /= 0) return
 
  ! Close the file.
  i = nf90_close (ncid)

  ! print a check on first variable
  print*, "Values in (1,1) cell of domain grid:"
  print*, "  lat      =", lat(1,1)   
  print*, "  lon      =", lon(1,1)   
  print*, "  elev     =", elev(1,1)   
  print*, "  mask     =", mask(1,1)   
  print*, "  grad_n_s =", grad_n_s(1,1)   
  print*, "  grad_w_e =", grad_w_e(1,1)   
  print*
 
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
 
end subroutine read_domain_grid
