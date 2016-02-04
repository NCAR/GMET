SUBROUTINE read_nc_grid (file_name, lat, lon, elev, grad_n_s, grad_w_e, mask, nx, ny, error)
  USE netcdf
  USE nrtype
  IMPLICIT NONE
 
 
  CHARACTER (LEN=*), INTENT (IN) :: file_name
 
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: lat (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: lon (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: elev (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: grad_n_s (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: grad_w_e (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: mask (:, :)
 
  INTEGER (I4B), INTENT (OUT) :: nx
  INTEGER (I4B), INTENT (OUT) :: ny
  INTEGER, INTENT (OUT) :: error
 
!local variables
 
  INTEGER :: ncid !netcdf file id
  INTEGER :: i
  INTEGER :: varid !variable id
  INTEGER (I4B), DIMENSION (nf90_max_var_dims) :: dimIds !dimension ids for dimensions of grid file
  INTEGER (I4B) :: ndims, nlat, nlon
 
  CHARACTER (LEN=*), PARAMETER :: lat_name = "latitude"
  CHARACTER (LEN=*), PARAMETER :: lon_name = "longitude"
  CHARACTER (LEN=*), PARAMETER :: elev_name = "elev"
  CHARACTER (LEN=*), PARAMETER :: grad_ns_name = "gradient_n_s"
  CHARACTER (LEN=*), PARAMETER :: grad_we_name = "gradient_w_e"
  CHARACTER (LEN=*), PARAMETER :: mask_name = "mask"
 
 
!code starts below
 
 
!open netcdf grid file
  CALL check (nf90_open(trim(file_name), nf90_nowrite, ncid), "File open error", error)
  IF (error /= 0) RETURN
 
  !inquire about latitude
  CALL check (nf90_inq_varid(ncid, lat_name, varid), "Latitude name error", error)
  IF (error /= 0) RETURN
 
  !get dimensions
  CALL check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimIds=dimIds), "Dimension inq error", error)
  IF (error /= 0 .OR. ndims /= 2) RETURN
 
  !get x,y dimensions
  CALL check (nf90_inquire_dimension(ncid, dimIds(1), len=nlon), "x dim error", error)
  CALL check (nf90_inquire_dimension(ncid, dimIds(2), len=nlat), "y dim error", error)
  IF (error /= 0) RETURN
 
  ny = nlat
  nx = nlon
 
  !allocate output variables
  ALLOCATE (lat(nlon, nlat))
  ALLOCATE (lon(nlon, nlat))
  ALLOCATE (elev(nlon, nlat))
  ALLOCATE (grad_n_s(nlon, nlat))
  ALLOCATE (grad_w_e(nlon, nlat))
  ALLOCATE (mask(nlon, nlat))
 
  !get lat,lon,elev,grad_n_s,grad_w_e
 
  !get latitude
  CALL check (nf90_get_var(ncid, varid, lat), "Latitude read error", error)
  IF (error /= 0) RETURN
 
 
  !inquire about longitdue
  CALL check (nf90_inq_varid(ncid, lon_name, varid), "Longitude name error", error)
  IF (error /= 0) RETURN
  !get it
  CALL check (nf90_get_var(ncid, varid, lon), "Longitdue read error", error)
  IF (error /= 0) RETURN
 
  !inquire about elevation
  CALL check (nf90_inq_varid(ncid, elev_name, varid), "Elevation name error", error)
  IF (error /= 0) RETURN
  !get it
  CALL check (nf90_get_var(ncid, varid, elev), "Elevation read error", error)
  IF (error /= 0) RETURN
 
  !inquire about n_s gradient
  CALL check (nf90_inq_varid(ncid, grad_ns_name, varid), "Grad_N_S name error", error)
  IF (error /= 0) RETURN
  !get it
  CALL check (nf90_get_var(ncid, varid, grad_n_s), "Grad_N_S read error", error)
  IF (error /= 0) RETURN
 
  !inquire about n_s gradient
  CALL check (nf90_inq_varid(ncid, grad_we_name, varid), "Grad_W_E name error", error)
  IF (error /= 0) RETURN
  !get it
  CALL check (nf90_get_var(ncid, varid, grad_w_e), "Grad_W_E read error", error)
  IF (error /= 0) RETURN
 
  !inquire about basin mask
  CALL check (nf90_inq_varid(ncid, mask_name, varid), "Mask name error", error)
  IF (error /= 0) RETURN
  !get it
  CALL check (nf90_get_var(ncid, varid, mask), "Mask read error", error)
  IF (error /= 0) RETURN
 
 
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
 
END SUBROUTINE read_nc_grid
