PROGRAM generate_ensembles
! ----------------------------------------------------------------------------------------
! Creator(s):
!   Andy Newman, 2013
!
! ----------------------------------------------------------------------------------------
! Purpose:
!
!   Driver for spatially correlated random field code from Martyn Clark
!   Generates ensebles of precipitation and temperature from regression step
!   For version 0 of CONUS ensemble product.  See Newman et al. 2015 J. Hydromet.
!
! ----------------------------------------------------------------------------------------
  USE netcdf !netcdf
  USE nrtype ! Numerical recipies types
  USE linkstruct !structure from topnet model for grid information
  USE gridweight !grid structure used by spcorr
  USE nr, ONLY: erf, erfcc ! Numerical Recipies error function
  USE qpe_namelist, ONLY: read_namelist !namelist module
  USE qpe_namelist, ONLY: nens, ntimes, start_time
  USE qpe_namelist, ONLY: out_name_base, qpe_nc_name, grid_name, clen
 
  IMPLICIT NONE
 
  INTERFACE
   SUBROUTINE read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
    USE nrtype
    CHARACTER (LEN=500), INTENT (IN) :: file_name
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
    INTEGER (I4B), INTENT (OUT) :: nx, ny
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_grid_list
 
   SUBROUTINE save_vars (pcp, tmean, trange, nx, ny, grdlat, grdlon, grdalt, times, file, error)
    USE netcdf
    USE nrtype
 
    REAL (SP), INTENT (IN) :: pcp (:, :, :), tmean (:, :, :), trange (:, :, :)
    INTEGER (I4B), INTENT (IN) :: nx, ny
    REAL (DP), INTENT (IN) :: grdlat (:), grdlon (:), grdalt (:)
    REAL (DP), INTENT (IN) :: times (:)
    CHARACTER (LEN=500), INTENT (IN) :: file
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE save_vars
 
   SUBROUTINE read_grid_qpe_nc_ens (file_name, var_name, var, lats, lons, auto_corr, tp_corr, times, tot_times, error)
    USE netcdf
    USE nrtype
 
    CHARACTER (LEN=*), INTENT (IN) :: file_name
    CHARACTER (LEN=*), INTENT (IN) :: var_name
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: var (:, :, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: lats (:, :), lons (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: times (:)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: auto_corr (:)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: tp_corr (:)
 
    INTEGER, INTENT (OUT) :: error
    INTEGER (I4B), INTENT (OUT) :: tot_times
   END SUBROUTINE read_grid_qpe_nc_ens
 
   SUBROUTINE normalize_x (x, mean, stdev)
    USE nrtype
    REAL (DP), INTENT (INOUT) :: x (:, :)
    REAL (DP), INTENT (OUT) :: mean
    REAL (DP), INTENT (OUT) :: stdev
   END SUBROUTINE normalize_x
 
   FUNCTION erfinv (x)
    USE nrtype
 
    REAL (SP), INTENT (IN) :: x
    REAL (SP) :: erfinv
   END FUNCTION erfinv
 
   SUBROUTINE read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
    USE netcdf
    USE nrtype
 
    CHARACTER (LEN=500), INTENT (IN) :: file_name
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), grad_e (:, :), mask (:, &
   & :)
    INTEGER (I4B), INTENT (OUT) :: nx, ny
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_nc_grid
 
 
  END INTERFACE
 
 
  ! Local variables
  INTEGER (I4B) :: i, j, k, igrd, istep, iens !counter variables
  INTEGER (I4B), DIMENSION (1:2) :: ORDER1 = (/ 2, 1 /)!order for reshape array
  INTEGER (I4B) :: ierr, jerr !error variables for various error checks
  INTEGER (I4B) :: NSPL1 ! # points (1st spatial dimension)
  INTEGER (I4B) :: NSPL2 ! # points (2nd spatial dimension)
  INTEGER (I4B) :: isp1 !first grid dimension location
  INTEGER (I4B) :: isp2 !second grid dimension location
  REAL (DP), DIMENSION (:, :), ALLOCATABLE :: RHO ! temporal correlation parameter
  REAL (DP), DIMENSION (:, :), ALLOCATABLE :: OLD_RANDOM ! previous correlated random field
  REAL (DP), DIMENSION (:, :), ALLOCATABLE :: pcp_RANDOM ! new correlated random field for pcp
  REAL (DP), DIMENSION (:, :), ALLOCATABLE :: tmean_RANDOM ! new correlated random field for tmean
  REAL (DP), DIMENSION (:, :), ALLOCATABLE :: trange_RANDOM ! new correlated random field for trange
 
  REAL (SP) :: acorr !value from scrf
  REAL (SP) :: aprob !probability from scrf
  REAL (SP) :: a_ra
  REAL (SP) :: aprob_ra
 
  REAL (DP) :: cprob !cdf value from scrf
  REAL (DP) :: amult !multiplier value to get actual precip from normalized value                                 ::
  REAL (DP) :: rn
  REAL (DP) :: ra
  REAL (DP) :: ra_err
  REAL (DP) :: cs
  REAL (DP) :: cprob_ra
 
  REAL (DP), ALLOCATABLE :: transform_exp (:)
  REAL (DP) :: transform
 
  CHARACTER (LEN=1024) :: out_name !base output name for netcdf files
  CHARACTER (LEN=128) :: suffix !suffix for ensemble member output
  CHARACTER (LEN=1024) :: var_name !name of netcdf variable grabbed from jason's netcdf file
  REAL (DP), ALLOCATABLE :: lon_out (:)! lon output to netcdf
  REAL (DP), ALLOCATABLE :: lat_out (:)! lat output to netcdf
  REAL (DP), ALLOCATABLE :: hgt_out (:)! hgt output to netcdf
  REAL (DP), ALLOCATABLE :: lat (:, :)
  REAL (DP), ALLOCATABLE :: lon (:, :)
  REAL (DP), ALLOCATABLE :: hgt (:, :)
  REAL (DP), ALLOCATABLE :: slp_e (:, :)
  REAL (DP), ALLOCATABLE :: slp_n (:, :)
  REAL (DP), ALLOCATABLE :: mask (:, :)
  REAL (DP), ALLOCATABLE :: weight (:, :)!weights from spcorr
  REAL (DP), ALLOCATABLE :: std (:, :)!std from spcorr
  REAL (DP), ALLOCATABLE :: var (:, :, :)!generic variable
  REAL (DP), ALLOCATABLE :: pcp (:, :, :)!output from qpe code, normalized precip
  REAL (DP), ALLOCATABLE :: pop (:, :, :)!output from qpe code, normalized pop
  REAL (DP), ALLOCATABLE :: pcp_error (:, :, :)!error from ols regression in qpe code
  REAL (DP), ALLOCATABLE :: tmean (:, :, :)
  REAL (DP), ALLOCATABLE :: tmean_error (:, :, :)
  REAL (DP), ALLOCATABLE :: trange (:, :, :)
  REAL (DP), ALLOCATABLE :: trange_error (:, :, :)
 
  REAL (DP), ALLOCATABLE :: pcp_2 (:, :, :)!output from qpe code, normalized precip
  REAL (DP), ALLOCATABLE :: pop_2 (:, :, :)!output from qpe code, normalized pop
  REAL (DP), ALLOCATABLE :: pcp_error_2 (:, :, :)!error from ols regression in qpe code
  REAL (DP), ALLOCATABLE :: tmean_2 (:, :, :)
  REAL (DP), ALLOCATABLE :: tmean_error_2 (:, :, :)
  REAL (DP), ALLOCATABLE :: trange_2 (:, :, :)
  REAL (DP), ALLOCATABLE :: trange_error_2 (:, :, :)
 
  REAL (DP), ALLOCATABLE :: lons (:, :)!lons array from qpe code
  REAL (DP), ALLOCATABLE :: lats (:, :)!lats array from qpe code
  REAL (DP), ALLOCATABLE :: times (:)!time vector from qpe code
  REAL (DP), ALLOCATABLE :: auto_corr (:)!lag-1 autocorrelation vector from qpe code
  REAL (DP), ALLOCATABLE :: tpc_corr (:)!temp-precip correlation vector from qpe code
  REAL (DP), ALLOCATABLE :: y_mean (:, :, :)!mean of transformed non-zero precip (at each timestep)
  REAL (DP), ALLOCATABLE :: y_std (:, :, :)!std dev of transformed non-zero precip (at each timestep)
  REAL (DP), ALLOCATABLE :: y_std_all (:, :, :)!std dev of transformed non-zero precip (at each timestep)
  REAL (DP), ALLOCATABLE :: y_min (:, :, :)!min of normalized transformed non-zero precip (at each timestep)
  REAL (DP), ALLOCATABLE :: y_max (:, :, :)!max of normalized transformed non-zero precip (at each timestep)
  REAL (SP), ALLOCATABLE :: pcp_out (:, :, :)!
  REAL (SP), ALLOCATABLE :: tmean_out (:, :, :)!
  REAL (SP), ALLOCATABLE :: trange_out (:, :, :)!
  INTEGER (I4B) :: nx, ny !grid size
  INTEGER (I4B) :: spl1_start, spl2_start !starting point of x,y grid
  INTEGER (I4B) :: spl1_count, spl2_count !length of x,y grid
  INTEGER (I4B) :: tot_times
 
  INTEGER :: ncid, varid, error
 
 
  TYPE (COORDS), POINTER :: grid !coordinate structure for grid
 
  TYPE (SPLNUM), DIMENSION (:, :), POINTER :: sp_pcp, sp_temp !structures of spatially correlated random field weights
 
  !read namelist in
  CALL read_namelist
 
  !set output file name from namelist
  out_name = out_name_base
 
  error = 0
  ierr = 0
  jerr = 0
 
  !read in netcdf grid file
  CALL read_nc_grid (grid_name, lat, lon, hgt, slp_n, slp_e, mask, nx, ny, error)
 
 
  IF (error .NE. 0) CALL exit_scrf (1, 'problem in read_nc_grid ')
 
  ALLOCATE (lat_out(nx*ny), lon_out(nx*ny), hgt_out(nx*ny), STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating for 1-d output variables')
 
 
 !allocate a few other variables
  ALLOCATE (pcp(nx, ny, ntimes))
  ALLOCATE (pop(nx, ny, ntimes))
  ALLOCATE (pcp_error(nx, ny, ntimes))
  ALLOCATE (tmean(nx, ny, ntimes))
  ALLOCATE (tmean_error(nx, ny, ntimes))
  ALLOCATE (trange(nx, ny, ntimes))
  ALLOCATE (trange_error(nx, ny, ntimes))
 
  ALLOCATE (pcp_2(nx, ny, ntimes))
  ALLOCATE (pop_2(nx, ny, ntimes))
  ALLOCATE (pcp_error_2(nx, ny, ntimes))
  ALLOCATE (tmean_2(nx, ny, ntimes))
  ALLOCATE (tmean_error_2(nx, ny, ntimes))
  ALLOCATE (trange_2(nx, ny, ntimes))
  ALLOCATE (trange_error_2(nx, ny, ntimes))
 
  ALLOCATE (y_mean(nx, ny, ntimes))
  ALLOCATE (y_std(nx, ny, ntimes))
  ALLOCATE (y_std_all(nx, ny, ntimes))
  ALLOCATE (y_min(nx, ny, ntimes))
  ALLOCATE (y_max(nx, ny, ntimes))
 
  ALLOCATE (times(ntimes))
  ALLOCATE (auto_corr(ntimes))
  ALLOCATE (tpc_corr(ntimes))
 
  PRINT *, 'Reading in Regression data, this will take a bit...'
 
  error = nf90_open (trim(qpe_nc_name), nf90_nowrite, ncid)
  IF (ierr /= 0) RETURN
 
  var_name = 'time'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, times, start= (/ start_time /), count= (/ ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'auto_corr'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, auto_corr, start= (/ start_time /), count= (/ ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'tp_corr'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, tpc_corr, start= (/ start_time /), count= (/ ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'pcp'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, pcp, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'pcp_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, pcp_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'pop'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, pop, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'pop_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, pop_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'pcp_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, pcp_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'pcp_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, pcp_error_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'tmean'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, tmean, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'tmean_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, tmean_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'tmean_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, tmean_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'tmean_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, tmean_error_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'trange'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, trange, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'trange_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, trange_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'trange_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, trange_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'trange_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, trange_error_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  var_name = 'ymean'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, y_mean, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'ystd'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, y_std, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'ystd_all'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, y_std_all, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'ymin'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, y_min, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
  var_name = 'ymax'
  error = nf90_inq_varid (ncid, var_name, varid)
  IF (ierr /= 0) RETURN
  error = nf90_get_var (ncid, varid, y_max, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  IF (ierr /= 0) RETURN
 
 
  error = nf90_close (ncid)
  IF (ierr /= 0) RETURN
 
!setup a few variables for spcorr structure
  NSPL1 = nx
  NSPL2 = ny
 
  spl1_start = 1
  spl2_start = 1
  spl1_count = nx
  spl2_count = ny
 
!allocate space for scrfs
  ALLOCATE (sp_pcp(NSPL1, NSPL2), STAT=ierr)
  IF (ierr .NE. 0) THEN
   CALL exit_scrf (1, 'problem deallocating space for sp_pcp ')
  END IF
 
  ALLOCATE (sp_temp(NSPL1, NSPL2), STAT=ierr)
  IF (ierr .NE. 0) THEN
   CALL exit_scrf (1, 'problem deallocating space for sp_temp ')
  END IF
 
 
  IF (allocated(RHO)) THEN
   DEALLOCATE (RHO, STAT=ierr)
   IF (ierr .NE. 0) CALL exit_scrf (1, 'problem deallocating space for rho ')
  END IF
  ALLOCATE (RHO(NSPL1, NSPL2), STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating space for rho ')
 
  IF (allocated(OLD_RANDOM)) THEN
   DEALLOCATE (OLD_RANDOM, STAT=ierr)
   IF (ierr .NE. 0) CALL exit_scrf (1, 'problem deallocating space for old_random ')
  END IF
  ALLOCATE (OLD_RANDOM(NSPL1, NSPL2), STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating space for old_random ')
 
  IF (allocated(pcp_RANDOM)) THEN
   DEALLOCATE (pcp_RANDOM, STAT=ierr)
   IF (ierr .NE. 0) CALL exit_scrf (1, 'problem deallocating space for pcp_random ')
  END IF
  ALLOCATE (pcp_RANDOM(NSPL1, NSPL2), STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating space for pcp_random ')
 
  IF (allocated(tmean_RANDOM)) THEN
   DEALLOCATE (tmean_RANDOM, STAT=ierr)
   IF (ierr .NE. 0) CALL exit_scrf (1, 'problem deallocating space for tmean_random ')
  END IF
  ALLOCATE (tmean_RANDOM(NSPL1, NSPL2), STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating space for tmean_random ')
 
  IF (allocated(trange_RANDOM)) THEN
   DEALLOCATE (trange_RANDOM, STAT=ierr)
   IF (ierr .NE. 0) CALL exit_scrf (1, 'problem deallocating space for trange_random ')
  END IF
  ALLOCATE (trange_RANDOM(NSPL1, NSPL2), STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating space for trange_random ')
 
 
  NULLIFY (grid)
  ALLOCATE (grid, STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating structure grid')
 
!place info into grid structure
  grid%idx%spl1_start = spl1_start
  grid%idx%spl2_start = spl2_start
  grid%idx%spl1_count = spl1_count
  grid%idx%spl2_count = spl2_count
 
  ! --------------------------------------------------------------------------------------
  ! allocate space for spatial arrays in grid structure
  ALLOCATE (grid%lat(spl1_count, spl2_count), grid%lon(spl1_count, spl2_count), grid%ELV(spl1_count, spl2_count), &
 & STAT=jerr)
  IF (ierr .NE. 0 .OR. jerr .NE. 0) CALL exit_scrf (1, ' problem allocating space for lat-lon-elv coordinates ')
 
 
  ALLOCATE (pcp_out(nx, ny, ntimes), tmean_out(nx, ny, ntimes), trange_out(nx, ny, ntimes), STAT=ierr)
  IF (ierr .NE. 0) CALL exit_scrf (1, 'problem allocating for 2-d output variables')
 
 
 
  pcp_out = 0.0
  tmean_out = 0.0
  trange_out = 0.0
 
  lon_out = pack (lon, .TRUE.)
  lat_out = pack (lat, .TRUE.)
  hgt_out = pack (hgt, .TRUE.)
 
  grid%lat = lat
  grid%lon = lon
  grid%ELV = hgt
 
 
  PRINT *, 'Generating weights for spatially correlated random field (SCRF)...'
 
  CALL spcorr_grd (NSPL1, NSPL2, grid)
 
  sp_pcp = spcorr
 
  CALL field_rand (NSPL1, NSPL2, pcp_RANDOM)
 
 
!setup sp_corr structure for temperature with larger correlation length
  clen = 800.0 !rough estimate based on observations
  CALL spcorr_grd (NSPL1, NSPL2, grid)
  sp_temp = spcorr
 
  CALL field_rand (NSPL1, NSPL2, tmean_RANDOM)
  CALL field_rand (NSPL1, NSPL2, trange_RANDOM)
 
  PRINT *, 'Done generating weights'
  PRINT *, 'Generating ensembles...'
 
  !transform power, shouldn't be hard-coded, but it is right now...
  transform = 4.0d0
 
 
    ! loop through the ensemble members
  DO iens = 1, nens
 
      !Loop through time
   DO istep = 1, ntimes
 
    DO igrd = 1, NSPL1 * NSPL2

     ! identify the (i,j) position of the igrd-th point
     isp1 = IORDER (igrd)
     isp2 = JORDER (igrd)
 
     !only compute values for valid grid points
     IF (grid%ELV(isp1, isp2) .GT.-300.0) THEN
 
   !find cumulative probability
      acorr = real (pcp_RANDOM(isp1, isp2), kind(SP)) / Sqrt (2._SP)
      aprob = erfcc (acorr)
      cprob = (2.d0-real(aprob, kind(DP))) / 2.d0
 
     ! check thresholds of slope fields to see which regression to use
     !For precipitation only
      IF (Abs(slp_n(isp1, isp2)) .LE. 3.6 .AND. Abs(slp_e(isp1, isp2)) .LE. 3.6) THEN
       pop (isp1, isp2, istep) = pop_2 (isp1, isp2, istep)
       pcp (isp1, isp2, istep) = pcp_2 (isp1, isp2, istep)
       pcp_error (isp1, isp2, istep) = pcp_error_2 (isp1, isp2, istep)
      END IF
 
   !for temperature don't use regression that included slope, only lat,lon,elevation based regression
      tmean (isp1, isp2, istep) = tmean_2 (isp1, isp2, istep)
      tmean_error (isp1, isp2, istep) = tmean_error_2 (isp1, isp2, istep)
      trange (isp1, isp2, istep) = trange_2 (isp1, isp2, istep)
      trange_error (isp1, isp2, istep) = trange_error_2 (isp1, isp2, istep)

      IF (cprob .LT. (1.0_DP-real(pop(isp1, isp2, istep), kind(DP)))) THEN !Don't generate precip
 
       pcp_out (isp1, isp2, istep) = 0.0d0
 
 
      ELSE !generate a precip amount
 
   !scale cumulative probability by regression pop
       cs = (cprob-(1.0_DP-real(pop(isp1, isp2, istep), kind(DP)))) / real (pop(isp1, isp2, istep), kind(DP))
 
   !convert cs to a z-score from standard normal
   !use erfinv
       IF (cs .LE. 3e-5) THEN
        rn = - 3.99
       ELSE IF (cs .GE. 0.99997) THEN
        rn = 3.99
       ELSE
        rn = Sqrt (2._SP) * erfinv ((2._SP*real(cs, kind(SP)))-1.0_SP)
       END IF
 
 
       ra = (real(pcp(isp1, isp2, istep), kind(DP))*real(y_std(isp1, isp2, istep), kind(DP))) + real (y_mean(isp1, &
      & isp2, istep), kind(DP)) + real (y_std(isp1, isp2, istep), kind(DP)) * rn * real (pcp_error(isp1, isp2, istep), &
      & kind(DP))
 
       IF (ra .GT. 0.0) THEN
        ra = ra ** transform
       ELSE
        ra = 0.01
       END IF
 
 
      !limit max value to y_max + pcp_error (max station value plus some portion of error)
       IF (ra .GT. (real(y_max(isp1, isp2, istep), kind(DP))+0.2*real(pcp_error(isp1, isp2, istep), &
      & kind(DP)))**transform) THEN
        ra = (real(y_max(isp1, isp2, istep), kind(DP))+0.2*real(pcp_error(isp1, isp2, istep), kind(DP))) ** transform
       END IF
 
 
       pcp_out (isp1, isp2, istep) = real (ra, kind(SP))
 
 
      END IF !end if statement for precip generation
 
 
     !tmean
      ra = real (tmean(isp1, isp2, istep), kind(DP)) + real (tmean_RANDOM(isp1, isp2), kind(DP)) * real &
     & (tmean_error(isp1, isp2, istep)/3.0, kind(DP))
      tmean_out (isp1, isp2, istep) = real (ra, kind(SP))
 
     !trange
      ra = real (trange(isp1, isp2, istep), kind(DP)) + real (trange_RANDOM(isp1, isp2), kind(DP)) * real &
     & (trange_error(isp1, isp2, istep)/3.0, kind(DP))
      trange_out (isp1, isp2, istep) = real (ra, kind(SP))
 
   !using +/- 3 std dev of uncertainty for temp gives unrealistic min and max exnsemble membner temps and diurnal ranges
   !ad hoc fix is to limit temp ensemble to roughly +/- 1 uncertainty range.  needs to be looked at further in future releases but
   !this at least gives reasonable temp results and covers the daymet, nldas, maurer spread for basins we've looked at
 
 
   !check for unrealistic and non-physical trange values
      IF (trange_out(isp1, isp2, istep) .GT. 40.0) THEN
       trange_out (isp1, isp2, istep) = 40.0
      ELSE IF (trange_out(isp1, isp2, istep) .LT. 2.0) THEN
       trange_out (isp1, isp2, istep) = 2.0
      END IF
 
   !check for unrealistic tmean values
      IF (tmean_out(isp1, isp2, istep) .GT. 35.0) THEN
       tmean_out (isp1, isp2, istep) = 35.0
      ELSE IF (tmean_out(isp1, isp2, istep) .LT.-35.0) THEN
       tmean_out (isp1, isp2, istep) = - 35.0
      END IF
 
 
     END IF !end valid elevation if check
 
    END DO !end loop for grid pts
 
    !Generate new SCRFs
    !generate new random numbers for tmean
    spcorr = sp_temp
    OLD_RANDOM = tmean_RANDOM
    CALL field_rand (NSPL1, NSPL2, tmean_RANDOM)
 
      !want to condition random numbers in the following way:
      !use the temp auto correlation to condition the tmean and trange
    tmean_RANDOM = OLD_RANDOM * auto_corr (1) + Sqrt (1-auto_corr(1)*auto_corr(1)) * tmean_RANDOM
 
      !generate new random numbers for trange
    OLD_RANDOM = trange_RANDOM
    CALL field_rand (NSPL1, NSPL2, trange_RANDOM)
    trange_RANDOM = OLD_RANDOM * auto_corr (1) + Sqrt (1-auto_corr(1)*auto_corr(1)) * trange_RANDOM
 
 
      !then use t-p correlation and trange_random to condition the precip random numbers
      !generate new random numbers for precip
    spcorr = sp_pcp
    CALL field_rand (NSPL1, NSPL2, pcp_RANDOM)
    pcp_RANDOM = trange_RANDOM * tpc_corr (1) + Sqrt (1-tpc_corr(1)*tpc_corr(1)) * pcp_RANDOM
 
   END DO !end time step loop
 
 
   PRINT *, 'Done with ensemble member: ', iens
   WRITE (suffix, '(I3.3)') iens
 
 
 
    !setup output name
   out_name = trim (out_name) // '_' // trim (suffix)
   PRINT *, trim (out_name)
 
 
 
    !save to netcdf file
   CALL save_vars (pcp_out, tmean_out, trange_out, nx, ny, lat_out, lon_out, hgt_out, &
  & times(start_time:start_time+ntimes-1), out_name, ierr)
 
 
   IF (ierr /= 0) RETURN
 
    !reset out_name
   out_name = out_name_base
 
  END DO !end ensemble member loop
 
 
END PROGRAM generate_ensembles
