program generate_ensembles
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
  use netcdf							  !netcdf
  use nrtype							  ! Numerical recipies types
  use linkstruct						  !structure from topnet model for grid information
  use gridweight					          !grid structure used by spcorr
  use nr, only: erf,erfcc                                         ! Numerical Recipies error function
  use qpe_namelist, only: read_namelist				  !namelist module
  use qpe_namelist, only: nens, ntimes, start_time
  use qpe_namelist, only: out_name_base, qpe_nc_name, grid_name,clen

  implicit none

  interface
    subroutine read_grid_list(file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
      use nrtype
      character(len=500), intent(in) :: file_name
      real(DP), allocatable, intent(out) :: lats(:), lons(:), alts(:), slp_n(:), slp_e(:)
      integer(I4B), intent(out) :: nx, ny
      integer, intent(out) :: error
    end subroutine read_grid_list

    subroutine save_vars(pcp, tmean, trange, nx, ny, grdlat, grdlon, grdalt, times, file, error)
      use netcdf
      use nrtype

      real(SP), intent(in) :: pcp(:,:,:), tmean(:,:,:), trange(:,:,:)
!      real(SP), intent(in) :: pcp(:,:), tmean(:,:), trange(:,:)
      integer(I4B), intent(in) :: nx, ny
      real(DP), intent(in) :: grdlat(:), grdlon(:), grdalt(:)
      real(DP), intent(in) :: Times(:)
      character (len = 500), intent(in) :: file
      integer, intent(out) :: error
    end subroutine save_vars

    subroutine read_grid_qpe_nc_ens(file_name,var_name,var,lats,lons,auto_corr,tp_corr, & 
                                times,tot_times,error)
      use netcdf
      use nrtype

      character (len = *), intent(in) :: file_name
      character (len = *), intent(in) :: var_name
      real(DP), allocatable, intent(out) :: var(:, :, :)
      real(DP), allocatable, intent(out) :: lats(:,:), lons(:,:)
      real(DP), allocatable, intent(out) :: times(:)
      real(DP), allocatable, intent(out) :: auto_corr(:)
      real(DP), allocatable, intent(out) :: tp_corr(:)

      integer, intent(out) :: error
      integer(I4B),intent(out) :: tot_times
    end subroutine read_grid_qpe_nc_ens

    subroutine normalize_x(x,mean,stdev)
      use nrtype
      real(DP), intent(inout) :: x(:,:)
      real(DP), intent(out) :: mean
      real(DP), intent(out) :: stdev
    end subroutine normalize_x

    function erfinv(x)
      use nrtype

      real(sp), intent(in) :: x
      real(sp)             :: erfinv
    end function erfinv

    subroutine read_nc_grid(file_name,lat,lon,elev,grad_n,grad_e,mask,nx,ny,error)
      use netcdf
      use nrtype

      character(len=500), intent(in)		:: file_name
      real(DP), allocatable, intent(out)	:: lat(:,:),lon(:,:),elev(:,:),grad_n(:,:),grad_e(:,:),mask(:,:)
      integer(I4B), intent(out)			:: nx,ny
      integer, intent(out)			:: error
    end subroutine read_nc_grid


  end interface


  ! Local variables
  INTEGER(I4B)                                :: i, j, k, igrd,istep,iens  !counter variables
  INTEGER(I4B), DIMENSION (1:2)               :: ORDER1 = (/ 2, 1 /) !order for reshape array
  integer(I4B)                                :: ierr,jerr        !error variables for various error checks
  INTEGER(I4B)                                :: NSPL1       ! # points (1st spatial dimension)
  INTEGER(I4B)                                :: NSPL2       ! # points (2nd spatial dimension)
  integer(i4B)                                :: isp1	     !first grid dimension location
  integer(i4B)                                :: isp2	     !second grid dimension location
  REAL(DP), DIMENSION(:,:), ALLOCATABLE       :: RHO        ! temporal correlation parameter
  REAL(DP), DIMENSION(:,:), ALLOCATABLE       :: OLD_RANDOM ! previous correlated random field
  REAL(DP), DIMENSION(:,:), ALLOCATABLE       :: pcp_RANDOM ! new correlated random field for pcp
  REAL(DP), DIMENSION(:,:), ALLOCATABLE       :: tmean_RANDOM ! new correlated random field for tmean
  REAL(DP), DIMENSION(:,:), ALLOCATABLE       :: trange_RANDOM ! new correlated random field for trange
  
  real(sp)                                    :: acorr !value from scrf
  real(sp)				      :: aprob !probability from scrf
  real(sp)			              :: a_ra
  real(sp)			              :: aprob_ra

  real(dp)				      :: cprob !cdf value from scrf
  real(dp)                                    :: amult !multiplier value to get actual precip from normalized value                                 ::
  real(dp)			              :: rn
  real(dp)			              :: ra
  real(dp)			              :: ra_err
  real(dp)			              :: cs
  real(dp)			              :: cprob_ra

  real(DP), allocatable			      :: transform_exp(:)
  real(DP)					      :: transform

  character (len=1024)                        :: out_name       !base output name for netcdf files
  character (len=128)			      :: suffix    !suffix for ensemble member output
  character (len=1024)                        :: var_name !name of netcdf variable grabbed from jason's netcdf file
  real(DP), allocatable                       :: lon_out(:) ! lon output to netcdf
  real(DP), allocatable                       :: lat_out(:) ! lat output to netcdf
  real(DP), allocatable                       :: hgt_out(:) ! hgt output to netcdf
  real(DP), allocatable				:: lat(:,:)
  real(DP), allocatable				:: lon(:,:)
  real(DP), allocatable				:: hgt(:,:)
  real(DP), allocatable				:: slp_e(:,:)
  real(DP), allocatable				:: slp_n(:,:)
  real(DP), allocatable				:: mask(:,:)
  real(DP), allocatable                       :: weight(:,:)	  !weights from spcorr
  real(DP), allocatable                       :: std(:,:)	  !std from spcorr
  real(DP), allocatable			      :: var(:,:,:)	!generic variable
  real(DP), allocatable			      :: pcp(:,:,:)	!output from qpe code, normalized precip
  real(DP), allocatable			      :: pop(:,:,:)	!output from qpe code, normalized pop
  real(DP), allocatable			      :: pcp_error(:,:,:)	!error from ols regression in qpe code
  real(DP), allocatable			      :: tmean(:,:,:)
  real(DP), allocatable			      :: tmean_error(:,:,:)
  real(DP), allocatable			      :: trange(:,:,:)
  real(DP), allocatable			      :: trange_error(:,:,:)

  real(DP), allocatable			      :: pcp_2(:,:,:)	!output from qpe code, normalized precip
  real(DP), allocatable			      :: pop_2(:,:,:)	!output from qpe code, normalized pop
  real(DP), allocatable			      :: pcp_error_2(:,:,:)	!error from ols regression in qpe code
  real(DP), allocatable			      :: tmean_2(:,:,:)
  real(DP), allocatable			      :: tmean_error_2(:,:,:)
  real(DP), allocatable			      :: trange_2(:,:,:)
  real(DP), allocatable			      :: trange_error_2(:,:,:)

  real(DP), allocatable                       :: lons(:,:)       !lons array from qpe code
  real(DP), allocatable                       :: lats(:,:)	  !lats array from qpe code
  real(DP), allocatable                       :: times(:)	  !time vector from qpe code
  real(DP), allocatable                       :: auto_corr(:)	  !lag-1 autocorrelation vector from qpe code
  real(DP), allocatable                       :: tpc_corr(:)	  !temp-precip correlation vector from qpe code
  real(DP), allocatable			      :: y_mean(:,:,:)       !mean of transformed non-zero precip (at each timestep)
  real(DP), allocatable			      :: y_std(:,:,:)       !std dev of transformed non-zero precip (at each timestep)
  real(DP), allocatable			      :: y_std_all(:,:,:)       !std dev of transformed non-zero precip (at each timestep)
  real(DP), allocatable			      :: y_min(:,:,:)       !min of normalized transformed non-zero precip (at each timestep)
  real(DP), allocatable			      :: y_max(:,:,:)       !max of normalized transformed non-zero precip (at each timestep)
  real(SP), allocatable                       :: pcp_out(:,:,:)	  !
  real(SP), allocatable                       :: tmean_out(:,:,:)	  !
  real(SP), allocatable                       :: trange_out(:,:,:)	  !
  integer(I4B)                                :: nx, ny    !grid size
  integer(I4B)				      :: spl1_start,spl2_start !starting point of x,y grid
  integer(I4B)				      :: spl1_count,spl2_count !length of x,y grid
  integer(I4B)				      :: tot_times

  integer				      :: ncid,varid,error


  TYPE(COORDS),POINTER                        :: grid	!coordinate structure for grid

  TYPE(SPLNUM),DIMENSION(:,:),POINTER	      :: sp_pcp,sp_temp !structures of spatially correlated random field weights

  !read namelist in
  call read_namelist

  !set output file name from namelist
  out_name = out_name_base

  error = 0
  ierr = 0
  jerr = 0

  !read in netcdf grid file
  call read_nc_grid(grid_name,lat,lon,hgt,slp_n,slp_e,mask,nx,ny,error)


  if(error .ne. 0) &
    call exit_scrf(1,'problem in read_nc_grid ')

  ALLOCATE(lat_out(nx*ny),lon_out(nx*ny),hgt_out(nx*ny),stat=ierr)
  if(ierr .ne. 0) &
    call exit_scrf(1,'problem allocating for 1-d output variables')


 !allocate a few other variables
  allocate(pcp(nx,ny,ntimes))
  allocate(pop(nx,ny,ntimes))
  allocate(pcp_error(nx,ny,ntimes))
  allocate(tmean(nx,ny,ntimes))
  allocate(tmean_error(nx,ny,ntimes))
  allocate(trange(nx,ny,ntimes))
  allocate(trange_error(nx,ny,ntimes))

  allocate(pcp_2(nx,ny,ntimes))
  allocate(pop_2(nx,ny,ntimes))
  allocate(pcp_error_2(nx,ny,ntimes))
  allocate(tmean_2(nx,ny,ntimes))
  allocate(tmean_error_2(nx,ny,ntimes))
  allocate(trange_2(nx,ny,ntimes))
  allocate(trange_error_2(nx,ny,ntimes))

  allocate(y_mean(nx,ny,ntimes))
  allocate(y_std(nx,ny,ntimes))
  allocate(y_std_all(nx,ny,ntimes))
  allocate(y_min(nx,ny,ntimes))
  allocate(y_max(nx,ny,ntimes))

  allocate(times(ntimes))
  allocate(auto_corr(ntimes))
  allocate(tpc_corr(ntimes))

  print *,'Reading in Regression data, this will take a bit...'

  error = nf90_open(trim(qpe_nc_name), nf90_nowrite, ncid)
  if(ierr /= 0) return

  var_name = 'time'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,times,start=(/start_time/),count=(/ntimes/))
  if(ierr /= 0) return

  var_name = 'auto_corr'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,auto_corr,start=(/start_time/),count=(/ntimes/))
  if(ierr /= 0) return

  var_name = 'tp_corr'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,tpc_corr,start=(/start_time/),count=(/ntimes/))
  if(ierr /= 0) return

  var_name = 'pcp'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,pcp,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'pcp_2'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,pcp_2,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'pop'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,pop,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'pop_2'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,pop_2,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'pcp_error'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,pcp_error,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'pcp_error_2'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,pcp_error_2,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'tmean'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,tmean,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'tmean_2'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,tmean_2,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'tmean_error'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,tmean_error,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'tmean_error_2'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,tmean_error_2,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'trange'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,trange,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'trange_2'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,trange_2,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'trange_error'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,trange_error,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'trange_error_2'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,trange_error_2,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  var_name = 'ymean'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,y_mean,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'ystd'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,y_std,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'ystd_all'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,y_std_all,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'ymin'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,y_min,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return

  var_name = 'ymax'
  error = nf90_inq_varid(ncid, var_name, varid)
  if(ierr /= 0) return
  error = nf90_get_var(ncid,varid,y_max,start=(/1,1,start_time/),count=(/nx,ny,ntimes/))
  if(ierr /= 0) return


  error = nf90_close(ncid)


!setup a few variables for spcorr structure
  nspl1 = nx
  nspl2 = ny

  spl1_start = 1
  spl2_start = 1
  spl1_count = nx
  spl2_count = ny

!allocate space for scrfs
  ALLOCATE(sp_pcp(NSPL1,NSPL2), STAT=IERR)
  if(ierr .ne. 0) then
      call exit_scrf(1,'problem deallocating space for sp_pcp ')
  endif

  ALLOCATE(sp_temp(NSPL1,NSPL2), STAT=IERR)
  if(ierr .ne. 0) then
      call exit_scrf(1,'problem deallocating space for sp_temp ')
  endif


  if (allocated(rho)) then
    deallocate(rho, stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem deallocating space for rho ')
    endif
    allocate(rho(nspl1,nspl2),stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem allocating space for rho ')

  if (allocated(old_random)) then
    deallocate(old_random, stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem deallocating space for old_random ')
    endif
    allocate(old_random(nspl1,nspl2),stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem allocating space for old_random ')

  if (allocated(pcp_random)) then
    deallocate(pcp_random, stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem deallocating space for pcp_random ')
    endif
    allocate(pcp_random(nspl1,nspl2),stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem allocating space for pcp_random ')

  if (allocated(tmean_random)) then
    deallocate(tmean_random, stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem deallocating space for tmean_random ')
    endif
    allocate(tmean_random(nspl1,nspl2),stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem allocating space for tmean_random ')

  if (allocated(trange_random)) then
    deallocate(trange_random, stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem deallocating space for trange_random ')
    endif
    allocate(trange_random(nspl1,nspl2),stat=ierr)
    if(ierr .ne. 0) &
      call exit_scrf(1,'problem allocating space for trange_random ')


  NULLIFY(grid)
  ALLOCATE(grid,STAT=IERR)
  IF (IERR.NE.0) &
   CALL EXIT_SCRF (1,'problem allocating structure grid')

!place info into grid structure
  grid%idx%spl1_start = spl1_start
  grid%idx%spl2_start = spl2_start
  grid%idx%spl1_count = spl1_count
  grid%idx%spl2_count = spl2_count

  ! --------------------------------------------------------------------------------------
  ! allocate space for spatial arrays in grid structure
  ALLOCATE(grid%LAT(SPL1_COUNT,SPL2_COUNT),grid%LON(SPL1_COUNT,SPL2_COUNT),&
           grid%ELV(SPL1_COUNT,SPL2_COUNT), STAT=JERR)
  IF (IERR.NE.0 .OR. JERR.NE.0) &
    CALL EXIT_SCRF(1,' problem allocating space for lat-lon-elv coordinates ')


ALLOCATE(pcp_out(nx,ny,ntimes),tmean_out(nx,ny,ntimes),trange_out(nx,ny,ntimes),stat=ierr)
  if(ierr .ne. 0) &
    call exit_scrf(1,'problem allocating for 2-d output variables')



  pcp_out = 0.0
  tmean_out = 0.0
  trange_out = 0.0

  lon_out = pack(lon,.true.)
  lat_out = pack(lat,.true.) 
  hgt_out = pack(hgt,.true.)

  grid%LAT=lat
  grid%LON=lon
  grid%ELV=hgt


  print *,'Generating weights for spatially correlated random field (SCRF)...'

  call spcorr_grd(nspl1, nspl2, grid)

  sp_pcp = spcorr

  call field_rand(nspl1, nspl2, pcp_random)


!setup sp_corr structure for temperature with larger correlation length
  clen = 800.0   !rough estimate based on observations
  call spcorr_grd(nspl1,nspl2,grid)
  sp_temp = spcorr

  call field_rand(nspl1, nspl2, tmean_random)
  call field_rand(nspl1, nspl2, trange_random)

  print *,'Done generating weights'
  print *,'Generating ensembles...'
 
  !transform power, shouldn't be hard-coded, but it is right now...
  transform =  4.0d0


    ! loop through the ensemble members
    DO IENS=1,nens

      !Loop through time
      DO ISTEP=1,ntimes

      DO IGRD=1,NSPL1*NSPL2


	! identify the (i,j) position of the igrd-th point
	ISP1 = IORDER(IGRD)
	ISP2 = JORDER(IGRD)

        !only compute values for valid grid points
	if(grid%elv(isp1,isp2) .gt. -300.0) then
 
	  !find cumulative probability
	  ACORR = real(pcp_random(ISP1,ISP2),kind(sp))/SQRT(2._sp)
	  APROB = ERFCC(ACORR)
	  CPROB = ( 2.d0 - real(APROB,kind(dp)) ) / 2.d0

     ! check thresholds of slope fields to see which regression to use 
     !For precipitation only
	if(abs(slp_n(isp1,isp2)) .le. 3.6 .and. abs(slp_e(isp1,isp2)) .le. 3.6) then
	  pop(isp1,isp2,istep) = pop_2(isp1,isp2,istep)
	  pcp(isp1,isp2,istep) = pcp_2(isp1,isp2,istep)
	  pcp_error(isp1,isp2,istep) = pcp_error_2(isp1,isp2,istep)
	endif

	!for temperature don't use regression that included slope, only lat,lon,elevation based regression
	tmean(isp1,isp2,istep) = tmean_2(isp1,isp2,istep)
	tmean_error(isp1,isp2,istep) = tmean_error_2(isp1,isp2,istep)
	trange(isp1,isp2,istep) = trange_2(isp1,isp2,istep)
	trange_error(isp1,isp2,istep) = trange_error_2(isp1,isp2,istep)

	if(cprob .lt. (1.0_dp-real(pop(isp1,isp2,istep),kind(dp)))) then 	!Don't generate precip

	  pcp_out(isp1,isp2,istep) = 0.0d0


	else  !generate a precip amount

	  !scale cumulative probability by regression pop
	  cs = (cprob - (1.0_dp-real(pop(isp1,isp2,istep),kind(dp)) ))/real(pop(isp1,isp2,istep),kind(dp))


	  !convert cs to a z-score from standard normal
	  !use erfinv
	  if(cs .le. 3e-5) then
	    rn = -3.99
	  elseif(cs .ge. 0.99997) then
	    rn = 3.99
	  else
	    rn = sqrt(2._sp) * erfinv((2._sp * real(cs,kind(sp))) - 1.0_sp)
	  endif


	  ra = ( real(pcp(isp1,isp2,istep),kind(dp))*real(y_std(isp1,isp2,istep),kind(dp)) ) + real(y_mean(isp1,isp2,istep),kind(dp)) & 
		+ real(y_std(isp1,isp2,istep),kind(dp))*rn*real(pcp_error(isp1,isp2,istep),kind(dp))

	  if(ra .gt. 0.0) then
	    ra = ra**transform
!	      ra = ra
	  else
!	      ra = real(y_min(isp1,isp2,istep),kind(dp))**transform
	    ra = 0.01
	  endif
      

!limit max value to y_max + pcp_error (max station value plus some portion of error)
	  if(ra .gt. (real(y_max(isp1,isp2,istep),kind(dp))+0.2*real(pcp_error(isp1,isp2,istep),kind(dp)))**transform) then
	    ra = (real(y_max(isp1,isp2,istep),kind(dp))+0.2*real(pcp_error(isp1,isp2,istep),kind(dp)))**transform
	  endif


	  pcp_out(isp1,isp2,istep) = real(ra,kind(sp))


	  endif   !end if statement for precip generation


          !tmean
	  ra = real(tmean(isp1,isp2,istep),kind(dp)) + real(tmean_random(isp1,isp2),kind(dp))*real(tmean_error(isp1,isp2,istep)/3.0,kind(dp))
	  tmean_out(isp1,isp2,istep) = real(ra,kind(sp))

	  !trange
	  ra = real(trange(isp1,isp2,istep),kind(dp)) + real(trange_random(isp1,isp2),kind(dp))*real(trange_error(isp1,isp2,istep)/3.0,kind(dp))
	  trange_out(isp1,isp2,istep) = real(ra,kind(sp))

	  !using +/- 3 std dev of uncertainty for temp gives unrealistic min and max exnsemble membner temps and diurnal ranges
	  !ad hoc fix is to limit temp ensemble to roughly +/- 1 uncertainty range.  needs to be looked at further in future releases but 
	  !this at least gives reasonable temp results and covers the daymet, nldas, maurer spread for basins we've looked at


	  !check for unrealistic and non-physical trange values
	  if(trange_out(isp1,isp2,istep) .gt. 40.0) then
	    trange_out(isp1,isp2,istep) = 40.0
	  elseif(trange_out(isp1,isp2,istep) .lt. 2.0) then
	    trange_out(isp1,isp2,istep) = 2.0
	  endif

	  !check for unrealistic tmean values
	  if(tmean_out(isp1,isp2,istep) .gt. 35.0) then
	    tmean_out(isp1,isp2,istep) = 35.0
	  elseif(tmean_out(isp1,isp2,istep) .lt. -35.0) then
	    tmean_out(isp1,isp2,istep) = -35.0
	  endif


	endif  !end valid elevation if check

      enddo	!end loop for grid pts
!      enddo

  !Generate new SCRFs
      !generate new random numbers for tmean
      spcorr = sp_temp
      old_random = tmean_random
      call field_rand(nspl1, nspl2, tmean_random)

      !want to condition random numbers in the following way:
      !use the temp auto correlation to condition the tmean and trange
      tmean_random = old_random*auto_corr(1) + sqrt(1-auto_corr(1)*auto_corr(1))*tmean_random

      !generate new random numbers for trange
      old_random = trange_random
      call field_rand(nspl1, nspl2, trange_random)
      trange_random = old_random*auto_corr(1) + sqrt(1-auto_corr(1)*auto_corr(1))*trange_random


      !then use t-p correlation and trange_random to condition the precip random numbers
      !generate new random numbers for precip
      spcorr = sp_pcp
      call field_rand(nspl1, nspl2, pcp_random)
      pcp_random = trange_random*tpc_corr(1) + sqrt(1-tpc_corr(1)*tpc_corr(1))*pcp_random

    enddo	!end time step loop


    print *,'Done with ensemble member: ',iens
    write(suffix,'(I3.3)') iens

   

    !setup output name
    out_name = trim(out_name)//'_'//trim(suffix)
    print *,trim(out_name)



    !save to netcdf file
    call save_vars(pcp_out, tmean_out, trange_out, nx, ny, lat_out, lon_out, hgt_out, &
                   Times(start_time:start_time+ntimes-1), out_name, ierr)


    if(ierr /= 0) return

    !reset out_name
    out_name = out_name_base

  enddo	  !end ensemble member loop


end program generate_ensembles