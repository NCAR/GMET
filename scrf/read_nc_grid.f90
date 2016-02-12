Subroutine read_nc_grid (file_name, lat, lon, elev, grad_n_s, grad_w_e, mask, nx, ny, error)
  Use netcdf
  Use nrtype
  Implicit None
 
 
  Character (Len=*), Intent (In) :: file_name
 
  Real (DP), Allocatable, Intent (Out) :: lat (:, :)
  Real (DP), Allocatable, Intent (Out) :: lon (:, :)
  Real (DP), Allocatable, Intent (Out) :: elev (:, :)
  Real (DP), Allocatable, Intent (Out) :: grad_n_s (:, :)
  Real (DP), Allocatable, Intent (Out) :: grad_w_e (:, :)
  Real (DP), Allocatable, Intent (Out) :: mask (:, :)
 
  Integer (I4B), Intent (Out) :: nx
  Integer (I4B), Intent (Out) :: ny
  Integer, Intent (Out) :: error
 
!local variables
 
  Integer :: ncid !netcdf file id
  Integer :: i
  Integer :: varid !variable id
  Integer (I4B), Dimension (nf90_max_var_dims) :: dimIds !dimension ids for dimensions of grid file
  Integer (I4B) :: ndims, nlat, nlon
 
  Character (Len=*), Parameter :: lat_name = "latitude"
  Character (Len=*), Parameter :: lon_name = "longitude"
  Character (Len=*), Parameter :: elev_name = "elev"
  Character (Len=*), Parameter :: grad_ns_name = "gradient_n_s"
  Character (Len=*), Parameter :: grad_we_name = "gradient_w_e"
  Character (Len=*), Parameter :: mask_name = "mask"
 
 
!code starts below
 
 
!open netcdf grid file
  Call check (nf90_open(trim(file_name), nf90_nowrite, ncid), "File open error", error)
  If (error /= 0) Return
 
  !inquire about latitude
  Call check (nf90_inq_varid(ncid, lat_name, varid), "Latitude name error", error)
  If (error /= 0) Return
 
  !get dimensions
  Call check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimIds=dimIds), "Dimension inq error", error)
  If (error /= 0 .Or. ndims /= 2) Return
 
  !get x,y dimensions
  Call check (nf90_inquire_dimension(ncid, dimIds(1), len=nlon), "x dim error", error)
  Call check (nf90_inquire_dimension(ncid, dimIds(2), len=nlat), "y dim error", error)
  If (error /= 0) Return
 
  ny = nlat
  nx = nlon
 
  !allocate output variables
  Allocate (lat(nlon, nlat))
  Allocate (lon(nlon, nlat))
  Allocate (elev(nlon, nlat))
  Allocate (grad_n_s(nlon, nlat))
  Allocate (grad_w_e(nlon, nlat))
  Allocate (mask(nlon, nlat))
 
  !get lat,lon,elev,grad_n_s,grad_w_e
 
  !get latitude
  Call check (nf90_get_var(ncid, varid, lat), "Latitude read error", error)
  If (error /= 0) Return
 
 
  !inquire about longitdue
  Call check (nf90_inq_varid(ncid, lon_name, varid), "Longitude name error", error)
  If (error /= 0) Return
  !get it
  Call check (nf90_get_var(ncid, varid, lon), "Longitdue read error", error)
  If (error /= 0) Return
 
  !inquire about elevation
  Call check (nf90_inq_varid(ncid, elev_name, varid), "Elevation name error", error)
  If (error /= 0) Return
  !get it
  Call check (nf90_get_var(ncid, varid, elev), "Elevation read error", error)
  If (error /= 0) Return
 
  !inquire about n_s gradient
  Call check (nf90_inq_varid(ncid, grad_ns_name, varid), "Grad_N_S name error", error)
  If (error /= 0) Return
  !get it
  Call check (nf90_get_var(ncid, varid, grad_n_s), "Grad_N_S read error", error)
  If (error /= 0) Return
 
  !inquire about n_s gradient
  Call check (nf90_inq_varid(ncid, grad_we_name, varid), "Grad_W_E name error", error)
  If (error /= 0) Return
  !get it
  Call check (nf90_get_var(ncid, varid, grad_w_e), "Grad_W_E read error", error)
  If (error /= 0) Return
 
  !inquire about basin mask
  Call check (nf90_inq_varid(ncid, mask_name, varid), "Mask name error", error)
  If (error /= 0) Return
  !get it
  Call check (nf90_get_var(ncid, varid, mask), "Mask read error", error)
  If (error /= 0) Return
 
 
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
 
End Subroutine read_nc_grid
