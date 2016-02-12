Program generate_ensembles
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
  Use netcdf !netcdf
  Use nrtype ! Numerical recipies types
  Use linkstruct !structure from topnet model for grid information
  Use gridweight !grid structure used by spcorr
  Use nr, Only: erf, erfcc ! Numerical Recipies error function
  Use qpe_namelist, Only: read_namelist !namelist module
  Use qpe_namelist, Only: nens, ntimes, start_time
  Use qpe_namelist, Only: out_name_base, qpe_nc_name, grid_name, clen
 
  Implicit None
 
  Interface
    Subroutine read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
      Use nrtype
      Character (Len=500), Intent (In) :: file_name
      Real (DP), Allocatable, Intent (Out) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
      Integer (I4B), Intent (Out) :: nx, ny
      Integer, Intent (Out) :: error
    End Subroutine read_grid_list
 
    Subroutine save_vars (pcp, tmean, trange, nx, ny, grdlat, grdlon, grdalt, times, file, error)
      Use netcdf
      Use nrtype
 
      Real (SP), Intent (In) :: pcp (:, :, :), tmean (:, :, :), trange (:, :, :)
      Integer (I4B), Intent (In) :: nx, ny
      Real (DP), Intent (In) :: grdlat (:), grdlon (:), grdalt (:)
      Real (DP), Intent (In) :: times (:)
      Character (Len=500), Intent (In) :: file
      Integer, Intent (Out) :: error
    End Subroutine save_vars
 
    Subroutine read_grid_qpe_nc_ens (file_name, var_name, var, lats, lons, auto_corr, tp_corr, times, tot_times, error)
      Use netcdf
      Use nrtype
 
      Character (Len=*), Intent (In) :: file_name
      Character (Len=*), Intent (In) :: var_name
      Real (DP), Allocatable, Intent (Out) :: var (:, :, :)
      Real (DP), Allocatable, Intent (Out) :: lats (:, :), lons (:, :)
      Real (DP), Allocatable, Intent (Out) :: times (:)
      Real (DP), Allocatable, Intent (Out) :: auto_corr (:)
      Real (DP), Allocatable, Intent (Out) :: tp_corr (:)
 
      Integer, Intent (Out) :: error
      Integer (I4B), Intent (Out) :: tot_times
    End Subroutine read_grid_qpe_nc_ens
 
    Subroutine normalize_x (x, mean, stdev)
      Use nrtype
      Real (DP), Intent (Inout) :: x (:, :)
      Real (DP), Intent (Out) :: mean
      Real (DP), Intent (Out) :: stdev
    End Subroutine normalize_x
 
    Function erfinv (x)
      Use nrtype
 
      Real (SP), Intent (In) :: x
      Real (SP) :: erfinv
    End Function erfinv
 
    Subroutine read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
      Use netcdf
      Use nrtype
 
      Character (Len=500), Intent (In) :: file_name
      Real (DP), Allocatable, Intent (Out) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), grad_e (:, :), mask &
     & (:, :)
      Integer (I4B), Intent (Out) :: nx, ny
      Integer, Intent (Out) :: error
    End Subroutine read_nc_grid
 
 
  End Interface
 
 
  ! Local variables
  Integer (I4B) :: i, j, k, igrd, istep, iens !counter variables
  Integer (I4B), Dimension (1:2) :: ORDER1 = (/ 2, 1 /)!order for reshape array
  Integer (I4B) :: ierr, jerr !error variables for various error checks
  Integer (I4B) :: NSPL1 ! # points (1st spatial dimension)
  Integer (I4B) :: NSPL2 ! # points (2nd spatial dimension)
  Integer (I4B) :: isp1 !first grid dimension location
  Integer (I4B) :: isp2 !second grid dimension location
  Real (DP), Dimension (:, :), Allocatable :: RHO ! temporal correlation parameter
  Real (DP), Dimension (:, :), Allocatable :: OLD_RANDOM ! previous correlated random field
  Real (DP), Dimension (:, :), Allocatable :: pcp_RANDOM ! new correlated random field for pcp
  Real (DP), Dimension (:, :), Allocatable :: tmean_RANDOM ! new correlated random field for tmean
  Real (DP), Dimension (:, :), Allocatable :: trange_RANDOM ! new correlated random field for trange
 
  Real (SP) :: acorr !value from scrf
  Real (SP) :: aprob !probability from scrf
  Real (SP) :: a_ra
  Real (SP) :: aprob_ra
 
  Real (DP) :: cprob !cdf value from scrf
  Real (DP) :: amult !multiplier value to get actual precip from normalized value                                 ::
  Real (DP) :: rn
  Real (DP) :: ra
  Real (DP) :: ra_err
  Real (DP) :: cs
  Real (DP) :: cprob_ra
 
  Real (DP), Allocatable :: transform_exp (:)
  Real (DP) :: transform
 
  Character (Len=1024) :: out_name !base output name for netcdf files
  Character (Len=128) :: suffix !suffix for ensemble member output
  Character (Len=1024) :: var_name !name of netcdf variable grabbed from jason's netcdf file
  Real (DP), Allocatable :: lon_out (:)! lon output to netcdf
  Real (DP), Allocatable :: lat_out (:)! lat output to netcdf
  Real (DP), Allocatable :: hgt_out (:)! hgt output to netcdf
  Real (DP), Allocatable :: lat (:, :)
  Real (DP), Allocatable :: lon (:, :)
  Real (DP), Allocatable :: hgt (:, :)
  Real (DP), Allocatable :: slp_e (:, :)
  Real (DP), Allocatable :: slp_n (:, :)
  Real (DP), Allocatable :: mask (:, :)
  Real (DP), Allocatable :: weight (:, :)!weights from spcorr
  Real (DP), Allocatable :: std (:, :)!std from spcorr
  Real (DP), Allocatable :: var (:, :, :)!generic variable
  Real (DP), Allocatable :: pcp (:, :, :)!output from qpe code, normalized precip
  Real (DP), Allocatable :: pop (:, :, :)!output from qpe code, normalized pop
  Real (DP), Allocatable :: pcp_error (:, :, :)!error from ols regression in qpe code
  Real (DP), Allocatable :: tmean (:, :, :)
  Real (DP), Allocatable :: tmean_error (:, :, :)
  Real (DP), Allocatable :: trange (:, :, :)
  Real (DP), Allocatable :: trange_error (:, :, :)
 
  Real (DP), Allocatable :: pcp_2 (:, :, :)!output from qpe code, normalized precip
  Real (DP), Allocatable :: pop_2 (:, :, :)!output from qpe code, normalized pop
  Real (DP), Allocatable :: pcp_error_2 (:, :, :)!error from ols regression in qpe code
  Real (DP), Allocatable :: tmean_2 (:, :, :)
  Real (DP), Allocatable :: tmean_error_2 (:, :, :)
  Real (DP), Allocatable :: trange_2 (:, :, :)
  Real (DP), Allocatable :: trange_error_2 (:, :, :)
 
  Real (DP), Allocatable :: lons (:, :)!lons array from qpe code
  Real (DP), Allocatable :: lats (:, :)!lats array from qpe code
  Real (DP), Allocatable :: times (:)!time vector from qpe code
  Real (DP), Allocatable :: auto_corr (:)!lag-1 autocorrelation vector from qpe code
  Real (DP), Allocatable :: tpc_corr (:)!temp-precip correlation vector from qpe code
  Real (DP), Allocatable :: y_mean (:, :, :)!mean of transformed non-zero precip (at each timestep)
  Real (DP), Allocatable :: y_std (:, :, :)!std dev of transformed non-zero precip (at each timestep)
  Real (DP), Allocatable :: y_std_all (:, :, :)!std dev of transformed non-zero precip (at each timestep)
  Real (DP), Allocatable :: y_min (:, :, :)!min of normalized transformed non-zero precip (at each timestep)
  Real (DP), Allocatable :: y_max (:, :, :)!max of normalized transformed non-zero precip (at each timestep)
  Real (SP), Allocatable :: pcp_out (:, :, :)!
  Real (SP), Allocatable :: tmean_out (:, :, :)!
  Real (SP), Allocatable :: trange_out (:, :, :)!
  Integer (I4B) :: nx, ny !grid size
  Integer (I4B) :: spl1_start, spl2_start !starting point of x,y grid
  Integer (I4B) :: spl1_count, spl2_count !length of x,y grid
  Integer (I4B) :: tot_times
 
  Integer :: ncid, varid, error
 
 
  Type (COORDS), Pointer :: grid !coordinate structure for grid
 
  Type (SPLNUM), Dimension (:, :), Pointer :: sp_pcp, sp_temp !structures of spatially correlated random field weights
 
  !read namelist in
  Call read_namelist
 
  !set output file name from namelist
  out_name = out_name_base
 
  error = 0
  ierr = 0
  jerr = 0
 
  !read in netcdf grid file
  Call read_nc_grid (grid_name, lat, lon, hgt, slp_n, slp_e, mask, nx, ny, error)
 
 
  If (error .Ne. 0) Call exit_scrf (1, 'problem in read_nc_grid ')
 
  Allocate (lat_out(nx*ny), lon_out(nx*ny), hgt_out(nx*ny), Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating for 1-d output variables')
 
 
 !allocate a few other variables
  Allocate (pcp(nx, ny, ntimes))
  Allocate (pop(nx, ny, ntimes))
  Allocate (pcp_error(nx, ny, ntimes))
  Allocate (tmean(nx, ny, ntimes))
  Allocate (tmean_error(nx, ny, ntimes))
  Allocate (trange(nx, ny, ntimes))
  Allocate (trange_error(nx, ny, ntimes))
 
  Allocate (pcp_2(nx, ny, ntimes))
  Allocate (pop_2(nx, ny, ntimes))
  Allocate (pcp_error_2(nx, ny, ntimes))
  Allocate (tmean_2(nx, ny, ntimes))
  Allocate (tmean_error_2(nx, ny, ntimes))
  Allocate (trange_2(nx, ny, ntimes))
  Allocate (trange_error_2(nx, ny, ntimes))
 
  Allocate (y_mean(nx, ny, ntimes))
  Allocate (y_std(nx, ny, ntimes))
  Allocate (y_std_all(nx, ny, ntimes))
  Allocate (y_min(nx, ny, ntimes))
  Allocate (y_max(nx, ny, ntimes))
 
  Allocate (times(ntimes))
  Allocate (auto_corr(ntimes))
  Allocate (tpc_corr(ntimes))
 
  Print *, 'Reading in Regression data, this will take a bit...'
 
  error = nf90_open (trim(qpe_nc_name), nf90_nowrite, ncid)
  If (ierr /= 0) Return
 
  var_name = 'time'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, times, start= (/ start_time /), count= (/ ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'auto_corr'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, auto_corr, start= (/ start_time /), count= (/ ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'tp_corr'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, tpc_corr, start= (/ start_time /), count= (/ ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'pcp'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, pcp, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'pcp_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, pcp_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'pop'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, pop, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'pop_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, pop_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'pcp_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, pcp_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'pcp_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, pcp_error_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'tmean'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, tmean, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'tmean_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, tmean_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'tmean_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, tmean_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'tmean_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, tmean_error_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'trange'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, trange, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'trange_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, trange_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'trange_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, trange_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'trange_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, trange_error_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  var_name = 'ymean'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, y_mean, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'ystd'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, y_std, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'ystd_all'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, y_std_all, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'ymin'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, y_min, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
  var_name = 'ymax'
  error = nf90_inq_varid (ncid, var_name, varid)
  If (ierr /= 0) Return
  error = nf90_get_var (ncid, varid, y_max, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes /))
  If (ierr /= 0) Return
 
 
  error = nf90_close (ncid)
  If (ierr /= 0) Return
 
!setup a few variables for spcorr structure
  NSPL1 = nx
  NSPL2 = ny
 
  spl1_start = 1
  spl2_start = 1
  spl1_count = nx
  spl2_count = ny
 
!allocate space for scrfs
  Allocate (sp_pcp(NSPL1, NSPL2), Stat=ierr)
  If (ierr .Ne. 0) Then
    Call exit_scrf (1, 'problem deallocating space for sp_pcp ')
  End If
 
  Allocate (sp_temp(NSPL1, NSPL2), Stat=ierr)
  If (ierr .Ne. 0) Then
    Call exit_scrf (1, 'problem deallocating space for sp_temp ')
  End If
 
 
  If (allocated(RHO)) Then
    Deallocate (RHO, Stat=ierr)
    If (ierr .Ne. 0) Call exit_scrf (1, 'problem deallocating space for rho ')
  End If
  Allocate (RHO(NSPL1, NSPL2), Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating space for rho ')
 
  If (allocated(OLD_RANDOM)) Then
    Deallocate (OLD_RANDOM, Stat=ierr)
    If (ierr .Ne. 0) Call exit_scrf (1, 'problem deallocating space for old_random ')
  End If
  Allocate (OLD_RANDOM(NSPL1, NSPL2), Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating space for old_random ')
 
  If (allocated(pcp_RANDOM)) Then
    Deallocate (pcp_RANDOM, Stat=ierr)
    If (ierr .Ne. 0) Call exit_scrf (1, 'problem deallocating space for pcp_random ')
  End If
  Allocate (pcp_RANDOM(NSPL1, NSPL2), Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating space for pcp_random ')
 
  If (allocated(tmean_RANDOM)) Then
    Deallocate (tmean_RANDOM, Stat=ierr)
    If (ierr .Ne. 0) Call exit_scrf (1, 'problem deallocating space for tmean_random ')
  End If
  Allocate (tmean_RANDOM(NSPL1, NSPL2), Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating space for tmean_random ')
 
  If (allocated(trange_RANDOM)) Then
    Deallocate (trange_RANDOM, Stat=ierr)
    If (ierr .Ne. 0) Call exit_scrf (1, 'problem deallocating space for trange_random ')
  End If
  Allocate (trange_RANDOM(NSPL1, NSPL2), Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating space for trange_random ')
 
 
  Nullify (grid)
  Allocate (grid, Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating structure grid')
 
!place info into grid structure
  grid%idx%spl1_start = spl1_start
  grid%idx%spl2_start = spl2_start
  grid%idx%spl1_count = spl1_count
  grid%idx%spl2_count = spl2_count
 
  ! --------------------------------------------------------------------------------------
  ! allocate space for spatial arrays in grid structure
  Allocate (grid%lat(spl1_count, spl2_count), grid%lon(spl1_count, spl2_count), grid%ELV(spl1_count, spl2_count), &
 & Stat=jerr)
  If (ierr .Ne. 0 .Or. jerr .Ne. 0) Call exit_scrf (1, ' problem allocating space for lat-lon-elv coordinates ')
 
 
  Allocate (pcp_out(nx, ny, ntimes), tmean_out(nx, ny, ntimes), trange_out(nx, ny, ntimes), Stat=ierr)
  If (ierr .Ne. 0) Call exit_scrf (1, 'problem allocating for 2-d output variables')
 
 
 
  pcp_out = 0.0
  tmean_out = 0.0
  trange_out = 0.0
 
  lon_out = pack (lon, .True.)
  lat_out = pack (lat, .True.)
  hgt_out = pack (hgt, .True.)
 
  grid%lat = lat
  grid%lon = lon
  grid%ELV = hgt
 
 
  Print *, 'Generating weights for spatially correlated random field (SCRF)...'
 
  Call spcorr_grd (NSPL1, NSPL2, grid)
 
  sp_pcp = spcorr
 
  Call field_rand (NSPL1, NSPL2, pcp_RANDOM)
 
 
!setup sp_corr structure for temperature with larger correlation length
  clen = 800.0 !rough estimate based on observations
  Call spcorr_grd (NSPL1, NSPL2, grid)
  sp_temp = spcorr
 
  Call field_rand (NSPL1, NSPL2, tmean_RANDOM)
  Call field_rand (NSPL1, NSPL2, trange_RANDOM)
 
  Print *, 'Done generating weights'
  Print *, 'Generating ensembles...'
 
  !transform power, shouldn't be hard-coded, but it is right now...
  transform = 4.0d0
 
 
    ! loop through the ensemble members
  Do iens = 1, nens
 
      !Loop through time
    Do istep = 1, ntimes
 
      Do igrd = 1, NSPL1 * NSPL2
 
     ! identify the (i,j) position of the igrd-th point
        isp1 = IORDER (igrd)
        isp2 = JORDER (igrd)
 
     !only compute values for valid grid points
        If (grid%ELV(isp1, isp2) .Gt.-300.0) Then
 
   !find cumulative probability
          acorr = real (pcp_RANDOM(isp1, isp2), kind(SP)) / Sqrt (2._SP)
          aprob = erfcc (acorr)
          cprob = (2.d0-real(aprob, kind(DP))) / 2.d0
 
     ! check thresholds of slope fields to see which regression to use
     !For precipitation only
          If (Abs(slp_n(isp1, isp2)) .Le. 3.6 .And. Abs(slp_e(isp1, isp2)) .Le. 3.6) Then
            pop (isp1, isp2, istep) = pop_2 (isp1, isp2, istep)
            pcp (isp1, isp2, istep) = pcp_2 (isp1, isp2, istep)
            pcp_error (isp1, isp2, istep) = pcp_error_2 (isp1, isp2, istep)
          End If
 
   !for temperature don't use regression that included slope, only lat,lon,elevation based regression
          tmean (isp1, isp2, istep) = tmean_2 (isp1, isp2, istep)
          tmean_error (isp1, isp2, istep) = tmean_error_2 (isp1, isp2, istep)
          trange (isp1, isp2, istep) = trange_2 (isp1, isp2, istep)
          trange_error (isp1, isp2, istep) = trange_error_2 (isp1, isp2, istep)
 
          If (cprob .Lt. (1.0_DP-real(pop(isp1, isp2, istep), kind(DP)))) Then !Don't generate precip
 
            pcp_out (isp1, isp2, istep) = 0.0d0
 
 
          Else !generate a precip amount
 
   !scale cumulative probability by regression pop
            cs = (cprob-(1.0_DP-real(pop(isp1, isp2, istep), kind(DP)))) / real (pop(isp1, isp2, istep), kind(DP))
 
   !convert cs to a z-score from standard normal
   !use erfinv
            If (cs .Le. 3e-5) Then
              rn = - 3.99
            Else If (cs .Ge. 0.99997) Then
              rn = 3.99
            Else
              rn = Sqrt (2._SP) * erfinv ((2._SP*real(cs, kind(SP)))-1.0_SP)
            End If
 
 
            ra = (real(pcp(isp1, isp2, istep), kind(DP))*real(y_std(isp1, isp2, istep), kind(DP))) + real (y_mean(isp1, &
           & isp2, istep), kind(DP)) + real (y_std(isp1, isp2, istep), kind(DP)) * rn * real (pcp_error(isp1, isp2, &
           & istep), kind(DP))
 
            If (ra .Gt. 0.0) Then
              ra = ra ** transform
            Else
              ra = 0.01
            End If
 
 
      !limit max value to y_max + pcp_error (max station value plus some portion of error)
            If (ra .Gt. (real(y_max(isp1, isp2, istep), kind(DP))+0.2*real(pcp_error(isp1, isp2, istep), &
           & kind(DP)))**transform) Then
              ra = (real(y_max(isp1, isp2, istep), kind(DP))+0.2*real(pcp_error(isp1, isp2, istep), kind(DP))) ** &
             & transform
            End If
 
 
            pcp_out (isp1, isp2, istep) = real (ra, kind(SP))
 
 
          End If !end if statement for precip generation
 
 
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
          If (trange_out(isp1, isp2, istep) .Gt. 40.0) Then
            trange_out (isp1, isp2, istep) = 40.0
          Else If (trange_out(isp1, isp2, istep) .Lt. 2.0) Then
            trange_out (isp1, isp2, istep) = 2.0
          End If
 
   !check for unrealistic tmean values
          If (tmean_out(isp1, isp2, istep) .Gt. 35.0) Then
            tmean_out (isp1, isp2, istep) = 35.0
          Else If (tmean_out(isp1, isp2, istep) .Lt.-35.0) Then
            tmean_out (isp1, isp2, istep) = - 35.0
          End If
 
 
        End If !end valid elevation if check
 
      End Do !end loop for grid pts
 
    !Generate new SCRFs
    !generate new random numbers for tmean
      spcorr = sp_temp
      OLD_RANDOM = tmean_RANDOM
      Call field_rand (NSPL1, NSPL2, tmean_RANDOM)
 
      !want to condition random numbers in the following way:
      !use the temp auto correlation to condition the tmean and trange
      tmean_RANDOM = OLD_RANDOM * auto_corr (1) + Sqrt (1-auto_corr(1)*auto_corr(1)) * tmean_RANDOM
 
      !generate new random numbers for trange
      OLD_RANDOM = trange_RANDOM
      Call field_rand (NSPL1, NSPL2, trange_RANDOM)
      trange_RANDOM = OLD_RANDOM * auto_corr (1) + Sqrt (1-auto_corr(1)*auto_corr(1)) * trange_RANDOM
 
 
      !then use t-p correlation and trange_random to condition the precip random numbers
      !generate new random numbers for precip
      spcorr = sp_pcp
      Call field_rand (NSPL1, NSPL2, pcp_RANDOM)
      pcp_RANDOM = trange_RANDOM * tpc_corr (1) + Sqrt (1-tpc_corr(1)*tpc_corr(1)) * pcp_RANDOM
 
    End Do !end time step loop
 
 
    Print *, 'Done with ensemble member: ', iens
    Write (suffix, '(I3.3)') iens
 
 
 
    !setup output name
    out_name = trim (out_name) // '_' // trim (suffix)
    Print *, trim (out_name)
 
 
 
    !save to netcdf file
    Call save_vars (pcp_out, tmean_out, trange_out, nx, ny, lat_out, lon_out, hgt_out, &
   & times(start_time:start_time+ntimes-1), out_name, ierr)
 
 
    If (ierr /= 0) Return
 
    !reset out_name
    out_name = out_name_base
 
  End Do !end ensemble member loop
 
 
End Program generate_ensembles
