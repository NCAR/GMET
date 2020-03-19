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
  use utim   !time utility routines
  use nrtype ! Numerical recipies types
  use linkstruct !structure from topnet model for grid information
  use gridweight !grid structure used by spcorr
  use nr, only: erf, erfcc ! Numerical Recipies error function
  use namelist_module, only: read_namelist !namelist module
  ! use namelist_module, only: nens, ntimes, start_time AW edited
  use namelist_module, only: start_ens, stop_ens, ntimes, start_time
  use namelist_module, only: out_forc_name_base, in_regr_name, grid_name, clen
  use namelist_module, only: time_mode
  use namelist_module, only: climo_path

 
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

    subroutine read_climo_precip (climo_path,current_month,climo_field,climo_out,error)
      use netcdf
      use nrtype
  
      character (len = *), intent(in)     :: climo_path      !path to climo netcdf files
      integer, intent(in)                 :: current_month   !integer month
      character (len = *), intent(in)     :: climo_field     !name of climo variable to read
      real(SP), intent(out)               :: climo_out(:,:,:)    !climo variable grid
      integer, intent(out)                :: error           !error code
    end subroutine read_climo_precip

    subroutine read_climo_temp(climo_path,current_month,climo_field,climo_out,error)
      use netcdf
      use nrtype

      character (len = *), intent(in)     :: climo_path      !path to climo netcdf files
      integer, intent(in)                 :: current_month   !integer month
      character (len = *), intent(in)     :: climo_field     !name of climo variable to read
      real(SP), intent(out)               :: climo_out(:,:)    !climo variable grid 
      integer, intent(out)                :: error           !error code
    end subroutine read_climo_temp

    subroutine read_climo_uncertainty(uncert_path,current_month,uncert_field,uncert_out,error)
      use netcdf
      use nrtype

      character (len = *), intent(in)     :: uncert_path      !path to uncert netcdf files
      integer, intent(in)                 :: current_month   !integer month
      character (len = *), intent(in)     :: uncert_field     !name of uncert variable to read
      real(SP), intent(out)               :: uncert_out(:,:)    !uncertainty variable grid
      integer, intent(out)                :: error           !error code
    end subroutine read_climo_uncertainty

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
  real (dp) :: amult !multiplier value to get actual precip from normalized value
  real (dp) :: rn
  real (dp) :: ra
  real (dp) :: ra_err
  real (dp) :: cs
  real (dp) :: cprob_ra
  integer(I4B)   :: cs_percentile !for climo precip distribution
 
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
  real(DP)               :: obs_max !precipitation limit using observations (used as cap for ensemble members)

  real (dp), allocatable :: lons (:, :)!lons array from qpe code
  real (dp), allocatable :: lats (:, :)!lats array from qpe code
  real (dp), allocatable :: times (:)!time vector from qpe code
  real (dp), allocatable :: auto_corr (:)!lag-1 autocorrelation vector from qpe code
  real (dp), allocatable :: tpc_corr (:)!temp-precip correlation vector from qpe code
  real (dp), allocatable :: obs_max_pcp (:, :, :) !max of non-0 pcp (each tstep)
  real (sp), allocatable :: pcp_out (:, :, :)!
  real (sp), allocatable :: tmean_out (:, :, :)!
  real (sp), allocatable :: trange_out (:, :, :)!
  integer (i4b) :: nx, ny !grid size
  integer (i4b) :: spl1_start, spl2_start !starting point of x,y grid
  integer (i4b) :: spl1_count, spl2_count !length of x,y grid
  integer (i4b) :: tot_times
  integer (i4b) :: ncid, dimid, varid, error
  integer (i4b) :: nTimesRegression 

  !climo grid variables
  real(SP),allocatable  :: climo_tmin(:,:,:)        !monthly climo tmax grids
  real(SP),allocatable  :: climo_tmax(:,:,:)        !monthly climo tmin grids
  real(SP),allocatable  :: climo_precip(:,:,:)        !climo precip grid for current day
  real(SP),allocatable  :: climo_tmean(:,:)         !climo tmean grid for current day
  real(SP),allocatable  :: climo_trange(:,:)        !climo trange grid for current day
  real(SP),allocatable  :: uncert_tmean(:,:)         !uncertainty tmean grid for current month
  real(SP),allocatable  :: uncert_trange(:,:)        !uncertainty trange grid for current month
  logical               :: first_climo = .FALSE.    !logical for first read
  integer               :: prev_month = -999        !previous month
  integer               :: current_month            !current month
  integer               :: prev_delta               !number of days from previous month for temporal interpolation
  integer               :: next_delta               !number of days to next month for temporal interpolation
  integer               :: year, day, hour, minute, second  !dates
  integer,dimension(12) :: month_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  character(len=2)      :: mnth_str    !string to contain current month we're in
  character(len=1024)   :: climo_file
  real(DP)              :: combined_error           !total error of daily anomaly uncertainty and climo uncertainty
  real(SP)              :: max_pcp                  !maximum allowable precip for a grid cell



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
 
 !allocate regression input variables variables
  allocate (pcp(nx, ny, ntimes))
  allocate (pop(nx, ny, ntimes))
  allocate (pcp_error(nx, ny, ntimes))
  allocate (tmean(nx, ny, ntimes))
  allocate (tmean_error(nx, ny, ntimes))
  allocate (trange(nx, ny, ntimes))
  allocate (trange_error(nx, ny, ntimes))
 
  allocate (obs_max_pcp(nx, ny, ntimes))
 
  allocate (times(ntimes))
  allocate (auto_corr(ntimes))
  allocate (tpc_corr(ntimes))
 
  print *, 'Reading in Regression data, this will take a bit...'
 
  error = nf90_open (trim(in_regr_name), nf90_nowrite, ncid)
  if (ierr /= 0) stop

  error = nf90_inq_dimid(ncid,'time',dimid)
  error = nf90_inquire_dimension(ncid,dimid,len=nTimesRegression)
  if(ierr /= 0)then; print *,'Error inquiring time dimension'; stop; endif

  if(nTimesRegression<ntimes)then
    print *,'Regression timesteps: ',nTimesRegression, ' less than namelist timesteps: ',ntimes
    stop
  endif
  if((start_time+ntimes)-1 > nTimesRegression)then
    print *,'Start timestep+number of timesteps: ',start_time+ntimes, ' greater than number of timesteps in regression output: ',nTimesRegression
    stop
  endif
 
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
  
  var_name = 'pop'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pop, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
 
  var_name = 'pcp_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, pcp_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'tmean'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, tmean, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
 
  var_name = 'tmean_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, tmean_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'trange'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, trange, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'trange_error'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, trange_error, start= (/ 1, 1, start_time /), count= (/ nx, ny, &
 & ntimes /))
  if (ierr /= 0) stop
 
  var_name = 'ymax'
  error = nf90_inq_varid (ncid, var_name, varid)
  if (ierr /= 0) stop
  error = nf90_get_var (ncid, varid, obs_max_pcp, start= (/ 1, 1, start_time /), count= (/ nx, ny, ntimes &
 & /))
  if (ierr /= 0) stop
 
  error = nf90_close (ncid)
  if (ierr /= 0) stop
 
  !sanity check on the observed maximum value in transformed anomaly space
  where(obs_max_pcp .gt. 5.)
    obs_max_pcp = 5.0
  end where

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
 

  if(trim(time_mode) .eq. 'daily' .or. trim(time_mode) .eq. 'DAILY') then
    ! set transform power, shouldn't be hard-coded, but it is for now...
    transform = 4.0d0
  elseif(trim(time_mode) .eq. 'climo' .or. trim(time_mode) .eq. 'CLIMO') then
    transform = 4.0d0
  elseif(trim(time_mode) .eq. 'daily_anom' .or. trim(time_mode) .eq. 'DAILY_ANOM') then
    transform = 3.0d0

    allocate(climo_precip(100,nx,ny))  !climo precip grid
    allocate(climo_tmean(nx,ny))   !climo tmean grid
    allocate(climo_trange(nx,ny))  !climo trange grid
    allocate(uncert_tmean(nx,ny))   !uncert tmean grid
    allocate(uncert_trange(nx,ny))  !uncert trange grid
  else
    print *, 'Incorrect time mode: ',trim(time_mode)
    print *, 'Current options are: daily, climo, daily_anom'
    stop
  end if

  ! ============ loop through the ensemble members ============
  ! do iens = 1, nens
  do iens = start_ens, stop_ens
 
    ! Loop through time
    do istep = 1, ntimes
 
      do igrd = 1, nspl1 * nspl2

        !check time mode and read climatological fields if necessary
        !climo fields only needed for daily anomaly time mode
        if(trim(time_mode) .eq. 'DAILY_ANOM' .or. trim(time_mode) .eq. 'daily_anom') then
          !read in climo grid if needed
          call unix_to_date(times(istep),year,current_month,day,hour,minute,second)

          if(current_month .ne. prev_month) then
            prev_month = current_month

            write(mnth_str,"(I0.2)") current_month

            !precipitation climo
            call read_climo_precip(climo_path,current_month,'pcp_sorted',climo_precip,error)

            !tmean climo
            call read_climo_temp(climo_path,current_month,'t_mean',climo_tmean,error)

            !trange climo
            call read_climo_temp(climo_path,current_month,'t_range',climo_trange,error)

            !read temperature uncertainty fields
            !tmean uncertainty
            call read_climo_uncertainty(climo_path,current_month,'t_mean_stddev',uncert_tmean,error)

            !trange uncertainty
            call read_climo_uncertainty(climo_path,current_month,'t_range_stddev',uncert_trange,error)

!            !convert climo precip to daily rates
!            climo_precip = climo_precip / real(month_days(current_month),kind(sp))
          end if !end of climo read check

        end if   !end of time_mode check
        

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
 
            cs_percentile = ceiling(cs*100.)
            if(cs_percentile .eq. 0) then
              cs_percentile = 1
            endif

            !time mode defines how ensemble values are generated
            if(trim(time_mode) .eq. 'daily_anom' .or. trim(time_mode) .eq. 'DAILY_ANOM') then
              ra = real(pcp(isp1,isp2,istep),kind(dp))+rn*real(pcp_error(isp1,isp2,istep),kind(dp))
              ra = ((ra*(1.0/transform))+1.0_dp)**transform

              pcp_out(isp1,isp2,istep) = ra*climo_precip(cs_percentile,isp1,isp2)

              max_pcp = obs_max_pcp(isp1,isp2,istep)+0.2*cs
              max_pcp = (max_pcp*(1.0/transform)+1.0_dp)**transform
              max_pcp = max_pcp*climo_precip(50,isp1,isp2)

              if(pcp_out(isp1,isp2,istep) .gt. max_pcp) then
                pcp_out(isp1,isp2,istep) = max_pcp
              endif

              if(pcp_out(isp1,isp2,istep) .lt. 0.1) then
                pcp_out(isp1,isp2,istep) = 0.1
              endif

            elseif(trim(time_mode) .eq. 'daily' .or. trim(time_mode) .eq. 'DAILY') then
              
              if(pcp_error(isp1,isp2,istep)<0.1) then
                pcp_error(isp1,isp2,istep) = 0.1
              endif
              ra =  real(pcp(isp1,isp2,istep),kind(dp))+rn*real(pcp_error(isp1,isp2,istep),kind(dp))
              ra = ((ra*(1.0/transform))+1.0_dp)**transform


              !limit max value to obs_max_pcp + pcp_error (max station value plus some portion of error)
              obs_max = 1.5*((((obs_max_pcp(isp1,isp2,istep)+0.2*cs)*(1.0/transform))+1.0)**transform)
              if(ra .gt. obs_max) then
                ra = obs_max
              endif

              pcp_out(isp1,isp2,istep) = ra

              if(pcp_out(isp1,isp2,istep) .lt. 0.1) then
                pcp_out(isp1,isp2,istep) = 0.1
              endif

            elseif(trim(time_mode) .eq. 'climo' .or. trim(time_mode) .eq. 'CLIMO') then
              ra = real(pcp(isp1,isp2,istep),kind(dp)) + rn*real(pcp_error(isp1,isp2,istep),kind(dp))
              ra = ((ra*(1.0_dp/transform))+1.0_dp)**transform

              if(ra .lt. 0.01) then
                pcp_out(isp1,isp2,istep) = 0.01
              else
                pcp_out(isp1,isp2,istep) = real(ra,kind(sp))
              endif
            else
              print *, 'Incorrect time mode: ',trim(time_mode)
              print *, 'Current options are: daily, climo, daily_anom'
              stop
            end if  !end if for precipitation time_mode check
          end if ! end IF statement for precip generation
 
          !time mode defines how ensemble values are generated
          if(trim(time_mode) .eq. 'daily_anom' .or. trim(time_mode) .eq. 'DAILY_ANOM') then
            !TMEAN
            combined_error = sqrt( real(tmean_error(isp1,isp2,istep),kind(dp))**2 + real(uncert_tmean(isp1,isp2),kind(dp))**2 )

            ra = ( real(tmean(isp1,isp2,istep),kind(dp))+real(climo_tmean(isp1,isp2),kind(dp)) ) + real(tmean_random(isp1,isp2),kind(dp))*combined_error
            tmean_out(isp1,isp2,istep) = ra

            !TRANGE
            combined_error = sqrt( real(trange_error(isp1,isp2,istep),kind(dp))**2 + real(uncert_trange(isp1,isp2),kind(dp))**2 )
 
            ra = ( real(trange(isp1,isp2,istep),kind(dp))+real(climo_trange(isp1,isp2),kind(dp)) ) + real(trange_random(isp1,isp2),kind(dp))*combined_error
            trange_out(isp1,isp2,istep) = ra

          elseif(trim(time_mode) .eq. 'daily' .or. trim(time_mode) .eq. 'DAILY') then
            !TMEAN
            ra = real(tmean(isp1,isp2,istep),kind(dp)) + real(tmean_random(isp1,isp2),kind(dp))*real(tmean_error(isp1,isp2,istep)/3.0,kind(dp))
            tmean_out(isp1,isp2,istep) = real(ra,kind(sp))
            !error term for tmean seems unrealistically large in many to all situations...  Limit to roughly +/- 1 std of error rather than ~+/- 3 std

            !trange
            ra = real(trange(isp1,isp2,istep),kind(dp)) + real(trange_random(isp1,isp2),kind(dp))*real(trange_error(isp1,isp2,istep)/3.0,kind(dp))
            trange_out(isp1,isp2,istep) = real(ra,kind(sp))
            !error term for trange seems unrealistically large in many to all situations...  Limit to roughly +/- 1 std of error rather than ~+/- 3 std
          elseif(trim(time_mode) .eq. 'climo' .or. trim(time_mode) .eq. 'CLIMO') then
            !TMEAN
            ra = real(tmean(isp1,isp2,istep),kind(dp)) + real(tmean_random(isp1,isp2),kind(dp))*real(tmean_error(isp1,isp2,istep),kind(dp))
            tmean_out(isp1,isp2,istep) = ra

            !TRANGE
            ra = real(trange(isp1,isp2,istep),kind(dp)) + real(trange_random(isp1,isp2),kind(dp))*real(trange_error(isp1,isp2,istep),kind(dp))
            trange_out(isp1,isp2,istep) = ra
          else
            print *, 'Incorrect time mode: ',trim(time_mode)
            print *, 'Current options are: daily, climo, daily_anom'
            stop
          end if  !end of time_mode check for temperature

          ! check for unrealistic and non-physical trange values
          if (trange_out(isp1, isp2, istep) .gt. 40.0) then
            trange_out (isp1, isp2, istep) = 40.0
          else if (trange_out(isp1, isp2, istep) .lt. 1.0) then
            trange_out (isp1, isp2, istep) = 1.0
          end if
 
          ! check for unrealistic tmean values
          if (tmean_out(isp1, isp2, istep) .gt. 50.0) then
            tmean_out (isp1, isp2, istep) = 50.0
          else if (tmean_out(isp1, isp2, istep) .lt.-50.0) then
            tmean_out (isp1, isp2, istep) = -50.0
          end if
 
        else  !if elevation not valid, fill with missing
          pcp_out(isp1,isp2,istep) = -999.0
          tmean_out(isp1,isp2,istep) = -999.0
          trange_out(isp1,isp2,istep) = -999.0
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
