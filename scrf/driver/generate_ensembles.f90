program generate_ensembles
! -----------------------------------------------------------------------------
! Creator(s):
!   Andy Newman, 2013
! Modified:
!   Andy Wood, 2016 -- adding namelist file as argument
!                   -- no longer hardwired; clean formatting
!                   -- adding documentation
!                   -- altered namelist args to specify ens range
! -----------------------------------------------------------------------------
! Purpose:
!   Driver for spatially correlated random field code from Martyn Clark
!   Generates ensebles of precipitation and temperature from regression step
!   For version 0 of CONUS ensemble product.  See Newman et al. 2015 JHM
! -----------------------------------------------------------------------------
 
  use netcdf !netcdf
  use nrtype ! Numerical recipies types
  use linkstruct !structure from topnet model for grid information
  use gridweight !grid structure used by spcorr
  use nr, only: erf, erfcc ! Numerical Recipies error function
  use namelist_module, only: read_namelist !namelist module
  ! use namelist_module, only: nens, ntimes, start_time AW edited
  use namelist_module, only: start_ens, stop_ens, ntimes, start_time
  use namelist_module, only: out_forc_name_base, in_regr_name, grid_name, clen
 
  implicit none
 
  interface
    subroutine read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
      use nrtype
      character (len=500), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_grid_list
 
    subroutine save_vars (pcp, tmean, trange, nx, ny, grdlat, grdlon, grdalt, times, file, error)
      use netcdf
      use nrtype
 
      real (sp), intent (in) :: pcp (:, :, :), tmean (:, :, :), trange (:, :, :)
      integer (i4b), intent (in) :: nx, ny
      real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
      real (dp), intent (in) :: times (:)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_vars
 
    subroutine read_grid_qpe_nc_ens (file_name, var_name, var, lats, lons, auto_corr, tp_corr, &
   & times, tot_times, error)
      use netcdf
      use nrtype
 
      character (len=*), intent (in) :: file_name
      character (len=*), intent (in) :: var_name
      real (dp), allocatable, intent (out) :: var (:, :, :)
      real (dp), allocatable, intent (out) :: lats (:, :), lons (:, :)
      real (dp), allocatable, intent (out) :: times (:)
      real (dp), allocatable, intent (out) :: auto_corr (:)
      real (dp), allocatable, intent (out) :: tp_corr (:)
 
      integer, intent (out) :: error
      integer (i4b), intent (out) :: tot_times
    end subroutine read_grid_qpe_nc_ens
 
    subroutine normalize_x (x, mean, stdev)
      use nrtype
      real (dp), intent (inout) :: x (:, :)
      real (dp), intent (out) :: mean
      real (dp), intent (out) :: stdev
    end subroutine normalize_x
 
    function erfinv (x)
      use nrtype
 
      real (sp), intent (in) :: x
      real (sp) :: erfinv
    end function erfinv
 
    subroutine read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
      use netcdf
      use nrtype
 
      character (len=*), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), &
     & grad_e (:, :), mask (:, :)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_nc_grid
  
  end interface
  ! ================== END of INTERFACES ===============
 
  ! Local variables
  integer (i4b) :: i, j, k, igrd, istep, iens !counter variables
  integer (i4b) :: nens  ! AW number of ensemble members to generate
  integer (i4b), dimension (1:2) :: order1 = (/ 2, 1 /)!order for reshape array
  integer (i4b) :: ierr, jerr !error variables for various error checks
  integer (i4b) :: nspl1 ! # points (1st spatial dimension)
  integer (i4b) :: nspl2 ! # points (2nd spatial dimension)
  integer (i4b) :: isp1  ! first grid dimension location
  integer (i4b) :: isp2  ! second grid dimension location
  real (dp), dimension (:, :), allocatable :: rho ! temporal correlation parameter
  real (dp), dimension (:, :), allocatable :: old_random ! previous correlated random field
  real (dp), dimension (:, :), allocatable :: pcp_random ! new correlated random field for pcp
  real (dp), dimension (:, :), allocatable :: tmean_random ! new correlated rand field, tmean
  real (dp), dimension (:, :), allocatable :: trange_random ! new correlated rand field, trange
 
  real (sp) :: acorr !value from scrf
  real (sp) :: aprob !probability from scrf
  real (sp) :: a_ra
  real (sp) :: aprob_ra
 
  real (dp) :: cprob !cdf value from scrf
  real (dp) :: amult !multiplier value to get actual precip from normalized value                                 ::
  real (dp) :: rn
  real (dp) :: ra
  real (dp) :: ra_err
  real (dp) :: cs
  real (dp) :: cprob_ra
 
  real (dp), allocatable :: transform_exp (:)
  real (dp) :: transform
 
  integer :: f ! AWW for command line argument read
  character (len=200) :: namelist_filename !AWW now an argument to the program
  character (len=1024) :: arg ! AWW command line arg for configuration file
  character (len=1024) :: out_name !base output name for netcdf files
  character (len=128) :: suffix !suffix for ensemble member output
  character (len=1024) :: var_name !name of netcdf variable grabbed from jason's netcdf file
  real (dp), allocatable :: lon_out (:)! lon output to netcdf
  real (dp), allocatable :: lat_out (:)! lat output to netcdf
  real (dp), allocatable :: hgt_out (:)! hgt output to netcdf
  real (dp), allocatable :: lat (:, :)
  real (dp), allocatable :: lon (:, :)
  real (dp), allocatable :: hgt (:, :)
  real (dp), allocatable :: slp_e (:, :)
  real (dp), allocatable :: slp_n (:, :)
  real (dp), allocatable :: mask (:, :)
  real (dp), allocatable :: weight (:, :)!weights from spcorr
  real (dp), allocatable :: std (:, :)!std from spcorr
  real (dp), allocatable :: var (:, :, :)!generic variable
  real (dp), allocatable :: pcp (:, :, :)!output from qpe code, normalized precip
  real (dp), allocatable :: pop (:, :, :)!output from qpe code, normalized pop
  real (dp), allocatable :: pcp_error (:, :, :)!error from ols regression in qpe code
  real (dp), allocatable :: tmean (:, :, :)
  real (dp), allocatable :: tmean_error (:, :, :)
  real (dp), allocatable :: trange (:, :, :)
  real (dp), allocatable :: trange_error (:, :, :)
 
  real (dp), allocatable :: pcp_2 (:, :, :)!output from qpe code, normalized precip
  real (dp), allocatable :: pop_2 (:, :, :)!output from qpe code, normalized pop
  real (dp), allocatable :: pcp_error_2 (:, :, :)!error from ols regression in qpe code
  real (dp), allocatable :: tmean_2 (:, :, :)
  real (dp), allocatable :: tmean_error_2 (:, :, :)
  real (dp), allocatable :: trange_2 (:, :, :)
  real (dp), allocatable :: trange_error_2 (:, :, :)
 
  real (dp), allocatable :: lons (:, :)!lons array from qpe code
  real (dp), allocatable :: lats (:, :)!lats array from qpe code
  real (dp), allocatable :: times (:)!time vector from qpe code
  real (dp), allocatable :: auto_corr (:)!lag-1 autocorrelation vector from qpe code
  real (dp), allocatable :: tpc_corr (:)!temp-precip correlation vector from qpe code
  ! Hongli remove
  !real (dp), allocatable :: y_mean (:, :, :) !mean of transformed non-0 pcp (at each timestep)
  !real (dp), allocatable :: y_std (:, :, :) !stdev of transformed non-0 pcpp (at each t-step)
  !real (dp), allocatable :: y_std_all (:, :, :) !stddev of transf. non-0 pcp (at each tstep)
  !real (dp), allocatable :: y_min (:, :, :) !min of normalized transf. non-0 pcp (each tstep)
  !real (dp), allocatable :: y_max (:, :, :) !max of normalized transf. non-0 pcp (each tstep)
  real (dp), allocatable :: obs_max_pcp (:, :, :) !max of normalized transf. non-0 pcp (each tstep)
  real (sp), allocatable :: pcp_out (:, :, :)!
  real (sp), allocatable :: tmean_out (:, :, :)!
  real (sp), allocatable :: trange_out (:, :, :)!
  integer (i4b) :: nx, ny !grid size
  integer (i4b) :: spl1_start, spl2_start !starting point of x,y grid
  integer (i4b) :: spl1_count, spl2_count !length of x,y grid
  integer (i4b) :: tot_times
  integer       :: ncid, varid, error
 
  type (coords), pointer :: grid !coordinate structure for grid
  type (splnum), dimension (:, :), pointer :: sp_pcp, sp_temp ! structures of spatially correlated random field weights

  ! ========== code starts below ==============================
 
  ! AWW: get namelist filename from command line (no longer hardwired)
  f = 0
  do
    call get_command_argument (f, arg)
    if (f .eq. 1) namelist_filename = arg
    if (len_trim(arg) == 0) exit
    f = f + 1
  end do

  ! set transform power, shouldn't be hard-coded, but it is for now...
  transform = 4.0d0
 
  ! read namelist in
  call read_namelist (namelist_filename)
 
  ! set output file name from namelist
  out_name = out_forc_name_base
  error = 0
  ierr = 0
  jerr = 0

  ! AW set number of ensembles to generate
  nens = stop_ens - start_ens + 1
  if(nens <= 0) call exit_scrf (1, 'number of ensembles to generate is 0 or less')
  print*, 'Generating ',nens,' ensembles from ',start_ens,' to ',stop_ens
  nens = stop_ens - start_ens + 1
  if(stop_ens <= start_ens) call exit_scrf (1, 'stop_ens equal or before start_ens')
  print*, 'Generating ',nens,' ensembles from ',start_ens,' to ',stop_ens
 
  !read in netcdf grid file
  call read_nc_grid (grid_name, lat, lon, hgt, slp_n, slp_e, mask, nx, ny, error)
 
  if (error .ne. 0) call exit_scrf (1, 'problem in read_nc_grid ')
 
  allocate (lat_out(nx*ny), lon_out(nx*ny), hgt_out(nx*ny), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating for 1-d output variables')
 
 !allocate a few other variables
  allocate (pcp(nx, ny, ntimes))
  allocate (pop(nx, ny, ntimes))
  allocate (pcp_error(nx, ny, ntimes))
  allocate (tmean(nx, ny, ntimes))
  allocate (tmean_error(nx, ny, ntimes))
  allocate (trange(nx, ny, ntimes))
  allocate (trange_error(nx, ny, ntimes))
 
  allocate (pcp_2(nx, ny, ntimes))
  allocate (pop_2(nx, ny, ntimes))
  allocate (pcp_error_2(nx, ny, ntimes))
  allocate (tmean_2(nx, ny, ntimes))
  allocate (tmean_error_2(nx, ny, ntimes))
  allocate (trange_2(nx, ny, ntimes))
  allocate (trange_error_2(nx, ny, ntimes))
 
  ! Hongli remove
  !allocate (y_mean(nx, ny, ntimes))
  !allocate (y_std(nx, ny, ntimes))
  !allocate (y_std_all(nx, ny, ntimes))
  !allocate (y_min(nx, ny, ntimes))
  !allocate (y_max(nx, ny, ntimes))
  allocate (obs_max_pcp(nx, ny, ntimes))
 
  allocate (times(ntimes))
  allocate (auto_corr(ntimes))
  allocate (tpc_corr(ntimes))
 
  print *, 'Reading in Regression data, this will take a bit...'
 
  error = nf90_open (trim(in_regr_name), nf90_nowrite, ncid)
  if (ierr /= 0) stop
 
  var_name = 'time'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, times, start= (/ start_time /), count= (/ ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'auto_corr'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, auto_corr, start= (/ start_time /), count= (/ ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'tp_corr'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, tpc_corr, start= (/ start_time /), count= (/ ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'pcp'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pcp, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
  
  var_name = 'pcp_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pcp_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop

  var_name = 'pop'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pop, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
 
  var_name = 'pop_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pop_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
 
  var_name = 'pcp_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pcp_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'pcp_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pcp_error_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'tmean'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, tmean, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
 
  var_name = 'tmean_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, tmean_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'tmean_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, tmean_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'tmean_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, tmean_error_2, start= (/ 1, 1, start_time /), count= (/ nx, &
 & ny, ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'trange'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, trange, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'trange_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, trange_2, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'trange_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, trange_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'trange_error_2'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, trange_error_2, start= (/ 1, 1, start_time /), count= (/ nx, &
 & ny, ntimes /))
  if (ierr /= 0) stop
 
  ! Hongli remove
  !var_name = 'ymean'
  !error = nf90_inq_varid (ncid, var_name, varid)
  !if (ierr /= 0) stop
  !error = nf90_get_var (ncid, varid, y_mean, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 !& ntimes /))
 ! if (ierr /= 0) stop
 
  !var_name = 'ystd'
  !error = nf90_inq_varid (ncid, var_name, varid)
  !if (ierr /= 0) stop
  !error = nf90_get_var (ncid, varid, y_std, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 !& /))
  !if (ierr /= 0) stop
 
  !var_name = 'ystd_all'
  !error = nf90_inq_varid (ncid, var_name, varid)
  !if (ierr /= 0) stop
  !error = nf90_get_var (ncid, varid, y_std_all, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 !& ntimes /))
  !if (ierr /= 0) stop
 
  !var_name = 'ymin'
  !error = nf90_inq_varid (ncid, var_name, varid)
  !if (ierr /= 0) stop
  !error = nf90_get_var (ncid, varid, y_min, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 !& /))
  !if (ierr /= 0) stop
 
  !var_name = 'ymax'
  !error = nf90_inq_varid (ncid, var_name, varid)
  !if (ierr /= 0) stop
  !error = nf90_get_var (ncid, varid, y_max, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 !& /))
  !if (ierr /= 0) stop

  var_name = 'obs_max_pcp'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, obs_max_pcp, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
 
  error = nf90_close (ncid)
  if (ierr /= 0) stop
 
  ! set up a few variables for spcorr structure
  nspl1 = nx
  nspl2 = ny
  spl1_start = 1
  spl2_start = 1
  spl1_count = nx
  spl2_count = ny
 
  ! allocate space for scrfs
  allocate (sp_pcp(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) then
    call exit_scrf (1, 'problem deallocating space for sp_pcp ')
  end if
 
  allocate (sp_temp(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) then
    call exit_scrf (1, 'problem deallocating space for sp_temp ')
  end if
 
  if (allocated(rho)) then
    deallocate (rho, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for rho ')
  end if
  allocate (rho(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for rho ')
 
  if (allocated(old_random)) then
    deallocate (old_random, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for old_random ')
  end if
  allocate (old_random(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for old_random ')
 
  if (allocated(pcp_random)) then
    deallocate (pcp_random, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for pcp_random ')
  end if
  allocate (pcp_random(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for pcp_random ')
 
  if (allocated(tmean_random)) then
    deallocate (tmean_random, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for tmean_random ')
  end if
  allocate (tmean_random(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for tmean_random ')
 
  if (allocated(trange_random)) then
    deallocate (trange_random, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for trange_random ')
  end if
  allocate (trange_random(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for trange_random ')
 
 
  nullify (grid)
  allocate (grid, stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating structure grid')
 
  ! place info into grid structure
  grid%idx%spl1_start = spl1_start
  grid%idx%spl2_start = spl2_start
  grid%idx%spl1_count = spl1_count
  grid%idx%spl2_count = spl2_count
 
  ! --------------------------------------------------------------------------------------
  ! allocate space for spatial arrays in grid structure
  allocate (grid%lat(spl1_count, spl2_count), grid%lon(spl1_count, spl2_count), &
 & grid%elv(spl1_count, spl2_count), stat=jerr)
  if (ierr .ne. 0 .or. jerr .ne. 0) call exit_scrf (1, ' problem allocating space for&
 & lat-lon-elev coordinates ')
 
  allocate (pcp_out(nx, ny, ntimes), tmean_out(nx, ny, ntimes), trange_out(nx, ny, ntimes), &
 & stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating for 2-d output variables')
  
  pcp_out = 0.0
  tmean_out = 0.0
  trange_out = 0.0
 
  lon_out = pack (lon, .true.)
  lat_out = pack (lat, .true.)
  hgt_out = pack (hgt, .true.)
 
  grid%lat = lat
  grid%lon = lon
  grid%elv = hgt
 
  print *, 'Generating weights for spatially correlated random field (SCRF)...'
 
  call spcorr_grd (nspl1, nspl2, grid)
  sp_pcp = spcorr
 
  call field_rand (nspl1, nspl2, pcp_random)
 
  ! setup sp_corr structure for temperature with larger correlation length
  clen = 800.0 !rough estimate based on observations
  call spcorr_grd (nspl1, nspl2, grid)
  sp_temp = spcorr
 
  call field_rand (nspl1, nspl2, tmean_random)
  call field_rand (nspl1, nspl2, trange_random)
 
  print *, 'Done generating weights...'
  print *, 'Generating ensembles...'
 
  ! ============ loop through the ensemble members ============
  ! do iens = 1, nens
  do iens = start_ens, stop_ens
 
    ! Loop through time
    do istep = 1, ntimes
 
      do igrd = 1, nspl1 * nspl2
 
        ! identify the (i,j) position of the igrd-th point
        isp1 = iorder (igrd)
        isp2 = jorder (igrd)
 
        ! only compute values for valid grid points
        if (grid%elv(isp1, isp2) .gt.-300.0) then
 
          ! find cumulative probability
          acorr = real (pcp_random(isp1, isp2), kind(sp)) / sqrt (2._sp)
          aprob = erfcc (acorr)
          cprob = (2.d0-real(aprob, kind(dp))) / 2.d0
 
          ! check thresholds of slope fields to see which regression to use

          ! For precipitation only
          !if (abs(slp_n(isp1, isp2)) .le. 3.6 .and. abs(slp_e(isp1, isp2)) .le. 3.6) then 
          ! hongli change slope threshold from 3.6 to 0.3 for nldas application ONLY
          if (abs(slp_n(isp1, isp2)) .le. 0.3 .and. abs(slp_e(isp1, isp2)) .le. 0.3) then
            pop (isp1, isp2, istep) = pop_2 (isp1, isp2, istep)
            pcp (isp1, isp2, istep) = pcp_2 (isp1, isp2, istep)
            pcp_error (isp1, isp2, istep) = pcp_error_2 (isp1, isp2, istep)
          end if
 
          ! For temperature don't use regression that included slope --
          !   only lat,lon,elevation based regression
          tmean (isp1, isp2, istep) = tmean_2 (isp1, isp2, istep)
          tmean_error (isp1, isp2, istep) = tmean_error_2 (isp1, isp2, istep)
          trange (isp1, isp2, istep) = trange_2 (isp1, isp2, istep)
          trange_error (isp1, isp2, istep) = trange_error_2 (isp1, isp2, istep)
 
          if (cprob .lt. (1.0_dp-real(pop(isp1, isp2, istep), kind(dp)))) then 
            ! Don't generate precip
            pcp_out (isp1, isp2, istep) = 0.0d0
 
          else 
            ! generate a precip amount
 
            ! scale cumulative probability by regression pop
            cs = (cprob-(1.0_dp-real(pop(isp1, isp2, istep), kind(dp)))) / real (pop(isp1, isp2, &
           & istep), kind(dp))
 
            ! convert cs to a z-score from standard normal
            ! use erfinv
            if (cs .le. 3e-5) then
              rn = - 3.99
            else if (cs .ge. 0.99997) then
              rn = 3.99
            else
              rn = sqrt (2._sp) * erfinv ((2._sp*real(cs, kind(sp)))-1.0_sp)
            end if
 
 
           ! ra = (real(pcp(isp1, isp2, istep), kind(dp))*real(y_std(isp1, isp2, istep), kind(dp))) &
           !& + real (y_mean(isp1, isp2, istep), kind(dp)) + real (y_std(isp1, isp2, istep), &
           !& kind(dp)) * rn * real (pcp_error(isp1, isp2, istep), kind(dp))
           ! Hongli modify
           ra = real(pcp(isp1, isp2, istep), kind(dp)) + rn * real(pcp_error(isp1, isp2, istep), kind(dp))           
 
            if (ra .gt. 0.0) then
              !ra = ra ** transform
              ! Hongli changed for box-cox back-transformation              
              ra = (ra*(1.0d0/transform)+1) ** transform 
            else
              ra = 0.01
            end if
 
 
            !! limit max value to y_max + pcp_error (max station value + some portion of error)
            !if (ra .gt. (real(y_max(isp1, isp2, istep), kind(dp))+0.2*real(pcp_error(isp1, isp2, &
           !& istep), kind(dp)))**transform) then
           !   ra = (real(y_max(isp1, isp2, istep), kind(dp))+0.2*real(pcp_error(isp1, isp2, istep), &
           !  //////(dp))) ** transform
           ! end if
            
            ! limit max value to obs_max_pcp + pcp_error (max station value + some portion of error)
            ! Hongli changed for box-cox back-transformation
            if (ra .gt. ((real(obs_max_pcp(isp1, isp2, istep), kind(dp))+0.2*real(pcp_error(isp1, isp2, &
           & istep), kind(dp)))*(1.0d0/transform)+1)**transform) then
              ra = ((real(obs_max_pcp(isp1, isp2, istep), kind(dp))+0.2*real(pcp_error(isp1, isp2, istep), &
             & kind(dp)))*(1.0d0/transform)+1) ** transform
            end if
 
            pcp_out (isp1, isp2, istep) = real (ra, kind(sp))
 
 
          end if 
          ! end IF statement for precip generation
 
 
          ! tmean
          !ra = real (tmean(isp1, isp2, istep), kind(dp)) + real (tmean_random(isp1, isp2), &
         !& kind(dp)) * real (tmean_error(isp1, isp2, istep)/3.0, kind(dp))! Hongli removed division of 3 as suggested by Andy Newman
         
         ra = real (tmean(isp1, isp2, istep), kind(dp)) + real (tmean_random(isp1, isp2), &
         & kind(dp)) * real (tmean_error(isp1, isp2, istep), kind(dp))
          tmean_out (isp1, isp2, istep) = real (ra, kind(sp)) 
 
          ! trange
          !ra = real (trange(isp1, isp2, istep), kind(dp)) + real (trange_random(isp1, isp2), &
         !& kind(dp)) * real (trange_error(isp1, isp2, istep)/3.0, kind(dp))! Hongli removed division of 3 as suggested by Andy Newman
         ra = real (trange(isp1, isp2, istep), kind(dp)) + real (trange_random(isp1, isp2), &
         & kind(dp)) * real (trange_error(isp1, isp2, istep), kind(dp))
          trange_out (isp1, isp2, istep) = real (ra, kind(sp))
 
          ! using +/- 3 std dev of uncertainty for temp gives unrealistic min and max ensemble 
          ! member temps and diurnal ranges
          ! ad hoc fix is to limit temp ensemble to roughly +/- 1 uncertainty range.  
          ! needs to be looked at further in future releases but
          ! this at least gives reasonable temp results and covers the daymet, 
          ! nldas, maurer spread for basins we've looked at
 
          ! check for unrealistic and non-physical trange values
          if (trange_out(isp1, isp2, istep) .gt. 40.0) then
            trange_out (isp1, isp2, istep) = 40.0
          else if (trange_out(isp1, isp2, istep) .lt. 2.0) then
            trange_out (isp1, isp2, istep) = 2.0
          end if
 
          ! check for unrealistic tmean values
          if (tmean_out(isp1, isp2, istep) .gt. 35.0) then
            tmean_out (isp1, isp2, istep) = 35.0
          else if (tmean_out(isp1, isp2, istep) .lt.-35.0) then
            tmean_out (isp1, isp2, istep) = - 35.0
          end if
 
 
        end if !end valid elevation if check
 
      end do !end loop for grid pts
 
      ! Generate new SCRFs
      ! generate new random numbers for tmean
      spcorr = sp_temp
      old_random = tmean_random
      call field_rand (nspl1, nspl2, tmean_random)
 
      ! want to condition random numbers in the following way:
      ! use the temp auto correlation to condition the tmean and trange
      tmean_random = old_random * auto_corr (1) + sqrt (1-auto_corr(1)*auto_corr(1)) * tmean_random
      ! generate new random numbers for trange
      old_random = trange_random
      call field_rand (nspl1, nspl2, trange_random)
      trange_random = old_random * auto_corr (1) + sqrt (1-auto_corr(1)*auto_corr(1)) * &
     & trange_random
 
      ! then use t-p correlation and trange_random to condition the precip random numbers
      ! generate new random numbers for precip
      spcorr = sp_pcp
      call field_rand (nspl1, nspl2, pcp_random)
      pcp_random = trange_random * tpc_corr (1) + sqrt (1-tpc_corr(1)*tpc_corr(1)) * pcp_random
 
    end do !end time step loop
 
    ! ============ now WRITE out the data file ============
    print *, 'Done with ensemble member: ', iens
    write (suffix, '(I3.3)') iens
    ! print *, 'Done with ensemble member: ', iens + start_ens - 1
    ! write (suffix, '(I3.3)') iens + start_ens - 1
 
    ! setup output name
    out_name = trim (out_name) // '.' // trim (suffix) // '.nc'
    print *, trim (out_name)
 
    ! save to netcdf file
    call save_vars (pcp_out, tmean_out, trange_out, nx, ny, lat_out, lon_out, hgt_out, &
   & times(start_time:start_time+ntimes-1), out_name, ierr)
  
    if (ierr /= 0) stop
 
    ! reset out_name
    out_name = out_forc_name_base
 
  end do !end ensemble member loop
 
end program generate_ensembles
