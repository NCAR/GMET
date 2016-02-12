subroutine read_nc_grid (file_name, lat, lon, elev, grad_n_s, grad_w_e, mask, nx, ny, error)
  use netcdf
  use type
  implicit none
 
 
  character (len=*), intent (in) :: file_name
 
  real (dp), allocatable, intent (out) :: lat (:, :)
  real (dp), allocatable, intent (out) :: lon (:, :)
  real (dp), allocatable, intent (out) :: elev (:, :)
  real (dp), allocatable, intent (out) :: grad_n_s (:, :)
  real (dp), allocatable, intent (out) :: grad_w_e (:, :)
  real (dp), allocatable, intent (out) :: mask (:, :)
 
  integer (i4b), intent (out) :: nx
  integer (i4b), intent (out) :: ny
  integer, intent (out) :: error
 
!local variables
 
  integer :: ncid !netcdf file id
  integer :: i
  integer :: varid !variable id
  integer (i4b), dimension (nf90_max_var_dims) :: dimids !dimension ids for dimensions of grid file
  integer (i4b) :: ndims, nlat, nlon
 
  character (len=*), parameter :: lat_name = "latitude"
  character (len=*), parameter :: lon_name = "longitude"
  character (len=*), parameter :: elev_name = "elev"
  character (len=*), parameter :: grad_ns_name = "gradient_n_s"
  character (len=*), parameter :: grad_we_name = "gradient_w_e"
  character (len=*), parameter :: mask_name = "mask"
 
 
!code starts below
 
 
!open netcdf grid file
  call check (nf90_open(file_name, nf90_nowrite, ncid), "File open error", error)
  if (error /= 0) return
 
  !inquire about latitude
  call check (nf90_inq_varid(ncid, lat_name, varid), "Latitude name error", error)
  if (error /= 0) return
 
  !get dimensions
  call check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), "Dimension inq error",&
 &  error)
  if (error /= 0 .or. ndims /= 2) return
 
  !get x,y dimensions
  call check (nf90_inquire_dimension(ncid, dimids(1), len=nlon), "x dim error", error)
  call check (nf90_inquire_dimension(ncid, dimids(2), len=nlat), "y dim error", error)
 
  print *, nlat, nlon
  ny = nlat
  nx = nlon
 
  !allocate output variables
  allocate (lat(nlon, nlat))
  allocate (lon(nlon, nlat))
  allocate (elev(nlon, nlat))
  allocate (grad_n_s(nlon, nlat))
  allocate (grad_w_e(nlon, nlat))
  allocate (mask(nlon, nlat))
 
  !get lat,lon,elev,grad_n_s,grad_w_e
 
  !get latitude
  call check (nf90_get_var(ncid, varid, lat), "Latitude read error", error)
  if (error /= 0) return
 
 
  !inquire about longitdue
  call check (nf90_inq_varid(ncid, lon_name, varid), "Longitude name error", error)
  if (error /= 0) return
  !get it
  call check (nf90_get_var(ncid, varid, lon), "Longitdue read error", error)
  if (error /= 0) return
 
  !inquire about elevation
  call check (nf90_inq_varid(ncid, elev_name, varid), "Elevation name error", error)
  if (error /= 0) return
  !get it
  call check (nf90_get_var(ncid, varid, elev), "Elevation read error", error)
  if (error /= 0) return
 
  !inquire about n_s gradient
  call check (nf90_inq_varid(ncid, grad_ns_name, varid), "Grad_N_S name error", error)
  if (error /= 0) return
  !get it
  call check (nf90_get_var(ncid, varid, grad_n_s), "Grad_N_S read error", error)
  if (error /= 0) return
 
  !inquire about n_s gradient
  call check (nf90_inq_varid(ncid, grad_we_name, varid), "Grad_W_E name error", error)
  if (error /= 0) return
  !get it
  call check (nf90_get_var(ncid, varid, grad_w_e), "Grad_W_E read error", error)
  if (error /= 0) return
 
  !inquire about basin mask
  call check (nf90_inq_varid(ncid, mask_name, varid), "Mask name error", error)
  if (error /= 0) return
  !get it
  call check (nf90_get_var(ncid, varid, mask), "Mask read error", error)
  if (error /= 0) return
 
 
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
 
end subroutine read_nc_grid
