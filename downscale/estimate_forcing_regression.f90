! AWW-2016Jan, modifications to handle time subsetting and reduce mem alloc, and clean up
!   renamed from estimate_precip; add also 'directory' var, changed some var names

subroutine estimate_forcing_regression (nPredict, gen_sta_weights, sta_weight_name, nwp_input_list, &
  & n_nwp, nwp_vars, x, z, ngrid, maxdistance, times, st_rec, end_rec, &
  & stnid, stnvar, directory, pcp, pop, pcperr, tmean, tmean_err, &
  & trange, trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, error, &
  & pcp_2, pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2)

  ! ==============================================================================================
  ! This routine is called during MODE 2 usage:  creates gridded ensembles from station/point data
  ! ==============================================================================================

  use string_mod
  use utim
  use type
  implicit none

  ! ===== start interfaces =======
  interface
    subroutine read_transform_exp (ntimes, file_name, texp)
      use type
      integer (i4b), intent (in) :: ntimes
      character (len=*), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: texp (:)
    end subroutine read_transform_exp

    subroutine read_station (stnvar, stnid, directory, st_rec, end_rec, vals, tair_vals, &
   & vals_miss, vals_miss_t, error)
      use type
      character (len=100), intent (in) :: stnvar
      character (len=100), intent (in) :: stnid
      character (len=500), intent (in) :: directory
      integer (i4b), intent (in) :: st_rec, end_rec
      real (dp), allocatable, intent (out) :: vals (:), tair_vals (:, :)
      logical, allocatable, intent (out) :: vals_miss (:), vals_miss_t (:)
      integer, intent (out) :: error
    end subroutine read_station

    subroutine compute_station_weights(sta_weight_name,ngrid,nstns,X,Z,search_distance, &
                                   sta_limit,sta_data,tair_data, &
                                   close_meta,close_meta_t,close_loc,close_loc_t, &
                                   close_count,close_count_t,close_weights,close_weights_t,error)
      use type
      !inputs
      character(len=500),intent(in) :: sta_weight_name   !name of station weight binary file
      integer(I4B), intent(in)      :: ngrid               !number of grid points
      integer(I4B), intent(in)      :: nstns               !number of stations
      real(DP), intent(in)          :: Z(:,:)              !grid metadata array
      real(DP), intent(in)          :: X(:,:)              !station metadata array
      real(DP), intent(in)          :: search_distance     !default station search distance
      integer(I4B), intent(in)      :: sta_limit           !maximum number of stations for a grid point
      real(DP), intent(in)          :: sta_data(:,:)       !station data values for precipitation
      real(DP), intent(in)          :: tair_data(:,:,:)    !station air temperature data
      !in/out
      real(DP), intent(inout)     :: close_meta(:,:,:)
      real(DP), intent(inout)     :: close_meta_t(:,:,:)
      integer(I4B), intent(inout) :: close_loc(:,:)
      integer(I4B), intent(inout) :: close_loc_t(:,:)
      integer(I4B), intent(inout) :: close_count(:)
      integer(I4B), intent(inout) :: close_count_t(:)
      real(DP), intent(inout)     :: close_weights(:,:)
      real(DP), intent(inout)     :: close_weights_t(:,:) 
      integer(I4B),intent(inout)  :: error
    end subroutine compute_station_weights

    subroutine write_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !input
                                close_weights_t,close_count,close_count_t,error) !input
      use type
      !inputs
      character(len=500), intent(in)    :: sta_weight_name
      real(DP), intent(in)      :: close_meta(:,:,:)
      real(DP), intent(in)      :: close_meta_t(:,:,:)
      integer(I4B), intent(in)  :: close_loc(:,:)
      integer(I4B), intent(in)  :: close_loc_t(:,:)
      real(DP), intent(in)      :: close_weights(:,:)
      real(DP), intent(in)      :: close_weights_t(:,:)
      integer(I4B), intent(in)  :: close_count(:)
      integer(I4B), intent(in)  :: close_count_t(:)
      integer(I4B), intent(inout):: error
    end subroutine write_station_weights

    subroutine read_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !output
                                close_weights_t,close_count,close_count_t,error) !output
      use type
      !input
      character(len=500), intent(in)    :: sta_weight_name
      !output
      real(DP), intent(out)      :: close_meta(:,:,:)
      real(DP), intent(out)      :: close_meta_t(:,:,:)
      integer(I4B), intent(out)  :: close_loc(:,:)
      integer(I4B), intent(out)  :: close_loc_t(:,:)
      real(DP), intent(out)      :: close_weights(:,:)
      real(DP), intent(out)      :: close_weights_t(:,:)
      integer(I4B), intent(out)  :: close_count(:)
      integer(I4B), intent(out)  :: close_count_t(:)
      integer(I4B), intent(inout):: error
    end subroutine read_station_weights


    subroutine read_nwp(currentTime,nwp_timestep_file,nPredict,n_nwp,nwp_vars,station_grid,x,z,error)
      use string_mod
      use utim
      use type
      character (len=2000), intent (in) :: nwp_timestep_file
      character (len=100), intent (in)  :: nwp_vars(:)
      integer(i4b), intent(in)          :: nPredict     !total number of predictors
      integer(i4b), intent(in)          :: n_nwp        !number of NWP predictors
      integer(i4b), intent(in)          :: station_grid(:)        !nearest grid point for every station location
      real(dp), intent(in)          :: currentTime  !current timestep unix time
      real (dp), intent (inout) :: x(:,:), z(:,:) !station and grid predictor matrices
      integer, intent (inout) :: error !error integer
    end subroutine read_nwp

    subroutine station_grid_correspondence(close_weights,close_weights_t,close_loc,close_loc_t,nSta,nearestGridpoint)
      use type

      real(dp), intent(in)        :: close_weights(:,:)
      real(dp), intent(in)        :: close_weights_t(:,:)
      integer(I4B), intent(in)    :: close_loc(:,:)
      integer(I4B), intent(in)    :: close_loc_t(:,:)
      integer(I4B), intent(in)    :: nSta
      integer(I4B), intent(out)   :: nearestGridpoint(:)
    end subroutine station_grid_correspondence

    subroutine normalize_x (x)
      use type
      real (dp), intent (inout) :: x (:, :)
    end subroutine normalize_x

    ! added AJN Sept 2013
    subroutine normalize_xv (x, weight, mean, stdev, stdev_all, smin, smax, yp)
      use type
      real (dp), intent (inout) :: x (:)
      real (dp), intent (in) :: weight (:)
      real (dp), intent (out) :: mean
      real (dp), intent (out) :: stdev
      real (dp), intent (out) :: stdev_all
      real (dp), intent (out) :: smin
      real (dp), intent (out) :: smax
      integer (i4b), intent (out) :: yp (:)
    end subroutine normalize_xv

    subroutine normalize_y (texp, y)
      use type
      real (dp), intent (in) :: texp !transform exponent
      real (dp), intent (inout) :: y (:)
    end subroutine normalize_y

    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares

    subroutine logistic_regression (x, y, tx, yp, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      integer (i4b), intent (in) :: yp (:)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine logistic_regression

    ! added AJN Sept 2013
    subroutine generic_corr (prcp_data, tair_data, lag, window, auto_corr, t_p_corr)
      use type
      real (dp), intent (in) :: prcp_data (:)
      real (dp), intent (in) :: tair_data (:, :)
      integer (i4b), intent (in) :: lag
      integer (i4b), intent (in) :: window
      real (dp), intent (out) :: auto_corr
      real (dp), intent (out) :: t_p_corr
    end subroutine generic_corr

    subroutine calc_distance (lat1, lon1, lat2, lon2, dist)
      use type
      implicit none
      real (dp), intent (in) :: lat1, lon1, lat2, lon2
      real (dp), intent (out) :: dist
    end subroutine calc_distance

    subroutine heapsort (n, ra, rn)
      use type
      implicit none
      integer (i4b), intent (in) :: n
      integer (i4b), dimension (:), intent (inout) :: rn
      real (dp), dimension (:), intent (inout) :: ra
    end subroutine heapsort

  end interface
  ! =========== end interfaces, start code =============

  real (dp), intent (inout) :: x (:, :), z (:, :)  ! station and grid point description arrays
  real (dp), intent (in) :: maxdistance         ! max distance for weight function
  integer (i4b), intent (in) :: ngrid           ! number of grid points
  real (dp), intent (in) :: times (:)!time step array

  ! AWW added next
  integer (i4b), intent (in) :: st_rec, end_rec

  character (len=100), intent (in)  :: stnid (:)!station id array
  character (len=100), intent (in)  :: stnvar !control file variables
  character (len=500), intent (in)  :: directory
  character (len = 500),intent(in)  :: gen_sta_weights     ! flag for generating station weight file
  character (len = 500),intent(in)  :: sta_weight_name     ! station weight file name
  character (len = 2000),intent(in) :: nwp_input_list      !file containing list of NWP predictor input files
  character (len=100), intent(in)   :: nwp_vars(:)    !list of nwp predictor variables

  character (len=2000)  ::  nwp_timestep_file        !name of a NWP predictor file

  integer(I4B), intent(in)          :: nPredict            ! number of total predictors
  integer(I4B), intent(in)          :: n_nwp               ! number of NWP predictors

  real (sp), allocatable, intent (out) :: pcp (:, :), pop (:, :), pcperr (:, :)!output variables for precipitation
  real (sp), allocatable, intent (out) :: tmean (:, :), tmean_err (:, :)!OLS tmean estimate and error
  real (sp), allocatable, intent (out) :: trange (:, :), trange_err (:, :)!OLS trange estimate and error

  real (sp), allocatable, intent (out) :: tmean_2 (:, :), tmean_err_2 (:, :)!OLS tmean estimate and error
  real (sp), allocatable, intent (out) :: trange_2 (:, :), trange_err_2 (:, :)!OLS trange estimate and error
  real (sp), allocatable, intent (out) :: pcp_2 (:, :), pop_2 (:, :), pcperr_2 (:, :)

  integer, intent (out) :: error ! integer error flag
  real (dp), intent (out) :: mean_autocorr (:)!mean autocorrelation from all stations over entire time period
  real (dp), intent (out) :: mean_tp_corr (:)!mean correlation for mean temp and precip

  ! vary at each grid point and time step
  real (dp), intent (out) :: y_mean (:, :), y_std (:, :)!std and mean of time step precipitation
  real (dp), intent (out) :: y_std_all (:, :)!std of time step precip including stations with zero precip
  real (dp), intent (out) :: y_min (:, :), y_max (:, :)!min & max  of normalized time step precipitation

  real (dp), allocatable :: y (:), b (:)
  !real (dp), allocatable :: twx (:, :), tx (:, :) ! not used

  real (dp), allocatable :: y_red (:)! reduced matrix for ...
  real (dp), allocatable :: x_red (:, :)! reduced matrix for ...
  real (dp), allocatable :: x_red_t (:, :)! reduced matrix for ...

  ! condensed these variables into just 4 that get re-used
  real (dp), allocatable :: twx_red (:, :), tx_red (:, :)!reduced matricies (orig)
  real (dp), allocatable :: twx_red_2 (:, :), tx_red_2 (:, :)!reduced matricies (dims 2)

  real (dp), allocatable :: w_base (:, :)!initial distance weight matrix
  real (dp), allocatable :: w_pcp_1d (:), w_temp_1d (:)
  integer (i4b), allocatable :: w_pcp_1d_loc (:), w_temp_1d_loc (:)
  !real(DP), allocatable :: w_pcp(:,:), w_temp(:,:) !distance weight matrices for a specific grid point
  real (dp), allocatable :: w_pcp_red (:, :), w_temp_red (:, :)!reduced distance weigth matricies

  real (dp), allocatable :: y_tmean (:), y_trange (:)!transformed station data arrays
  real (dp), allocatable :: y_tmean_red (:), y_trange_red (:)!transformed station data arrays
  real (dp), allocatable :: stn_prcp (:), prcp_data (:, :), tair_data (:, :, :), stn_tair (:, :)! orig stn data arrays
  real (dp), allocatable :: auto_corr (:)   ! lag-1 autocorrelation for stations over entire time period used
  real (dp), allocatable :: t_p_corr (:)    ! correlation between temp and precip
  integer (i4b), allocatable :: yp (:)      ! binary for logistic regression
  integer (i4b), allocatable :: yp_red (:)  ! reduced binary for logistic regression

  logical, allocatable :: stn_miss (:), stn_miss_t (:) ! missing value logical arrays

  real (dp) :: errsum, wgtsum, sta_temp
  real (dp) :: auto_corr_sum, tp_corr_sum
  real (dp) :: step_mean, step_std, step_std_all, step_min, step_max ! timestep statistics
  real (dp) :: ss_tot, ss_res ! r-squared and variance correction

  integer (i4b) :: xsize !size of second dimension of input X array
  integer (i4b) :: ntimes, nstns
  integer (i4b) :: t, i, g, ndata, nodata
  integer (i4b) :: ndata_t, nodata_t
  integer (i4b) :: lag, window
  integer (i4b) :: auto_cnt, tp_cnt

  integer(I4B)  :: nBase  !number of geophysical predictors
  integer(I4B),allocatable  :: noSlopePredicts(:)

  ! variables for tracking closest N stations for precipitation
  integer (i4b), parameter :: sta_limit = 30
  integer (i4b), allocatable :: close_loc (:, :)
  integer (i4b), allocatable :: close_count (:)
  real (dp), allocatable :: close_weights (:, :)
  real (dp), allocatable :: close_meta (:, :, :)
  real (dp) :: max_distance
  real (dp), parameter :: search_distance = 1000.0

  ! variables for tracking closest N stations for temperature
  integer (i4b), allocatable :: close_loc_t (:, :)
  integer (i4b), allocatable :: close_count_t (:)
  real (dp), allocatable :: close_weights_t (:, :)
  real (dp), allocatable :: close_meta_t (:, :, :)
  real (dp) :: max_distance_t
  real (dp) :: tmp_weight

  integer(I4B), allocatable :: nearestGridpoint(:)

  integer (i4b) :: slope_flag_pcp
  integer (i4b) :: slope_flag_temp

  ! variables to check for singular matrix
  real (dp), allocatable :: tmp (:, :)
  real (dp), allocatable :: vv (:)

  ! variables for timing code AJN
  integer (i4b) :: t1, t2, count_rate
  integer (i4b) :: tg1, tg2

  ! variables for keeping track of max_distance modifications
  integer (i4b), allocatable :: expand_flag (:), expand_flag_t (:)
  real (dp), allocatable :: expand_dist (:), expand_dist_t (:)

  !==============================================================!
  !                     code starts below here                   !
  !==============================================================!

  nstns = size (stnid)
  ntimes = size (times)
  xsize = size (x, 2)

  ! allocate variables
  allocate (y(nstns))
  allocate (y_tmean(nstns), y_trange(nstns))
  allocate (prcp_data(nstns, ntimes))
  allocate (tair_data(2, nstns, ntimes))
  allocate (w_pcp_red(sta_limit, sta_limit))
  allocate (w_temp_red(sta_limit, sta_limit))
  allocate (y_red(sta_limit))
  allocate (y_tmean_red(sta_limit), y_trange_red(sta_limit))
  allocate (x_red(sta_limit, xsize))
  allocate (x_red_t(sta_limit, xsize))
  allocate (w_pcp_1d(sta_limit))
  allocate (w_temp_1d(sta_limit))
  allocate (w_pcp_1d_loc(sta_limit))
  allocate (w_temp_1d_loc(sta_limit))
  allocate (tmp(nPredict, nPredict))
  allocate (vv(nPredict))
  allocate (pcp(ngrid, ntimes))
  allocate (pop(ngrid, ntimes))
  allocate (pcperr(ngrid, ntimes))
  allocate (tmean(ngrid, ntimes))
  allocate (tmean_err(ngrid, ntimes))
  allocate (trange(ngrid, ntimes))
  allocate (trange_err(ngrid, ntimes))
  allocate (pcp_2(ngrid, ntimes))
  allocate (pop_2(ngrid, ntimes))
  allocate (pcperr_2(ngrid, ntimes))
  allocate (tmean_2(ngrid, ntimes))
  allocate (tmean_err_2(ngrid, ntimes))
  allocate (trange_2(ngrid, ntimes))
  allocate (trange_err_2(ngrid, ntimes))
  allocate (auto_corr(nstns))
  allocate (t_p_corr(nstns))
  allocate (yp(nstns))
  allocate (yp_red(sta_limit))

  ! station limit arrays (precip)
  allocate (close_weights(ngrid, sta_limit))
  allocate (close_loc(ngrid, sta_limit))
  allocate (close_meta(5, ngrid, sta_limit))
  allocate (close_count(ngrid))

  ! station limit arrays (temp)
  allocate (close_weights_t(ngrid, sta_limit))
  allocate (close_loc_t(ngrid, sta_limit))
  allocate (close_meta_t(5, ngrid, sta_limit))
  allocate (close_count_t(ngrid))

  ! base weight array
  allocate (w_base(ngrid, nstns))

  !nearest grid index to stations
  allocate (nearestGridpoint(nstns))
 
  !predictor index array
  allocate (noSlopePredicts(nPredict-2))

  ! max_dist tracking variables
  allocate (expand_dist(ngrid), expand_flag(ngrid))
  allocate (expand_dist_t(ngrid), expand_flag_t(ngrid))

  ! initializations
  pcp = 0.0d0
  pop = 0.0d0
  pcperr = 0.0d0
  auto_corr = 0.0d0

  tmean = 0.0d0
  trange = 0.0d0
  tmean_err = 0.0d0
  trange_err = 0.0d0

  auto_corr_sum = 0.0d0
  auto_cnt = 0
  tp_corr_sum = 0.0d0
  tp_cnt = 0

  w_base = 0.0d0

  expand_dist = 0.0d0
  expand_flag = 0
  expand_dist_t = 0.0d0
  expand_flag_t = 0


  ! ================= LOOP OVER STATIONS ============
  ! this part calls subroutines that calc various  correlations
  ! can do autocorrelations and correlation between temperature and precipitation
  ! uses an n-day moving average (window) to remove "monthly" cycle from temp
  ! and computes autocorrelation on the anomalies

  do i = 1, nstns, 1

    call read_station (stnvar, stnid(i), directory, st_rec, end_rec, stn_prcp, stn_tair, &
   & stn_miss, stn_miss_t, error)

    print*, "first value check: "
    print*, "   prcp(1): ",stn_prcp(1)
    print*, "   tavg(1): ",stn_tair(1,1)
    print*, "   trng(1): ",stn_tair(2,1)
    print*

    prcp_data (i, :) = stn_prcp
    tair_data (1, i, :) = stn_tair (1, :)
    tair_data (2, i, :) = stn_tair (2, :)

    ! check data if needed
    ! print *, '-----------'
    ! print *, 'precip',prcp_data(i,:),'MISS',stn_miss
    ! print *, '-----------'
    ! print *,'tmin',tair_data(1,i,:),'MISS',stn_miss_t
    ! print *, '-----------'

    ! compute mean autocorrelation for all stations and all times
    lag = 1
    window = 31 ! AWW:  hardwired parameter; should bring out
    call generic_corr (prcp_data(i, :), tair_data(:, i, :), lag, window, auto_corr(i), t_p_corr(i))
    ! print *,auto_corr(i)

    ! check for values outside of -1 to 1
    ! stations with incomplete data are set to -999
    print *, 'station id & auto_corr'
    print *, stnid(i), auto_corr(i)
    if (auto_corr(i) .ge.-1.0 .and. auto_corr(i) .le. 1.0) then
      auto_corr_sum = auto_corr_sum + auto_corr (i)
      auto_cnt = auto_cnt + 1
    end if
    if (t_p_corr(i) .ge.-1.0 .and. t_p_corr(i) .le. 1.0) then
      tp_corr_sum = tp_corr_sum + t_p_corr (i)
      tp_cnt = tp_cnt + 1
    end if

    deallocate (stn_miss_t)  ! must be allocated within read_station
    deallocate (stn_miss)
    deallocate (stn_prcp)
    deallocate (stn_tair)

  end do
  ! =========== end station read loop ============

  error = 0 ! AWW:  why is this set?  not used again in subroutine

  ! AWW: some checks
  print *, 'auto_cnt, tp_cnt=', auto_cnt, tp_cnt
  if (auto_cnt == 0 .or. tp_cnt == 0) then
    print *, 'ERROR:  autocorr or crosscorr (TxP) could not be calculated due to lack of matching p&
   &airs'
    stop
  end if
  mean_autocorr = auto_corr_sum / real (auto_cnt, kind(dp))
  mean_tp_corr = tp_corr_sum / real (tp_cnt, kind(dp))

  print *, ' '
  print *, '===================================================='
  print *, 'Temp lag-1 autocorrelation: ', mean_autocorr (1)
  print *, 'Temp-precip correlation: ', mean_tp_corr (1)
  print *, '===================================================='
  print *, ' '

  call system_clock (t1, count_rate)

  ! ========= LOOP OVER GRID CELLS ==================
  !Create station-grid cell weight matrices before time stepping
  print *, 'Generating base weight matrix and finding nearest stations for each gridpoint'

  if(gen_sta_weights .eq. "TRUE" .or. gen_sta_weights .eq. "true") then
    call compute_station_weights(sta_weight_name,ngrid,nstns,X,Z,search_distance, & !input
                                 sta_limit,prcp_data,tair_data, & !input
                                 close_meta,close_meta_t,close_loc,close_loc_t, &  !output
                                 close_count,close_count_t,close_weights,close_weights_t,error) !output

    if(error /= 0) then
       return
    endif
    call write_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !input
                                close_weights_t,close_count,close_count_t,error) !input
    if(error /= 0) then
       return
    endif
  else
    call read_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !output
                                close_weights_t,close_count,close_count_t,error) !output
    if(error /= 0) then
       return
    endif
  endif

  call system_clock (t2, count_rate)
  print *, 'Elapsed time for weight generation: ', real (t2-t1) / real (count_rate)

  ! AWW-Feb2016:  just allocate grids once time, and re-use in code below
  allocate (twx_red(nPredict, sta_limit))    ! these have dim1 = 6
  allocate (tx_red(nPredict, sta_limit))
  allocate (twx_red_2(nPredict-2, sta_limit))  ! these are for no slope calcs, have dim1 = 4
  allocate (tx_red_2(nPredict-2, sta_limit))


  !open NWP predcitor file list
  open(unit=55,file=trim(nwp_input_list),form='formatted',status='old',iostat=error)
  if( error /= 0) then
    print *,'Error opening NWP predictor file list: ',trim(nwp_input_list)
    stop
  end if

  !create index array of closest grid point to every station
  call station_grid_correspondence(close_weights,close_weights_t,close_loc,close_loc_t,nstns,nearestGridpoint)

  nBase = nPredict-n_nwp
  noSlopePredicts(1:4) = (/1,2,3,4/)
!print *,'here ',nBase,nPredict
  do i = 1,n_nwp,1
    noSlopePredicts(nBase-2+i) = nBase+i
!print *,'here ',nBase,nPredict,i,noSlopePredicts
  end do

  ! =========== now LOOP through all TIME steps and populate grids ===============
  do t = 1, ntimes, 1

    call system_clock (tg1, count_rate)
    print *, "TIME STEP = ", times (t), " (", t, "/", ntimes, ")"

    ! --- assign vectors of station values for prcp, temp, for current time step
    do i = 1, nstns, 1
      y (i) = prcp_data (i, t)
      y_tmean (i) = (tair_data(1, i, t)+tair_data(2, i, t)) / 2.0d0
      y_trange (i) = abs (tair_data(2, i, t)-tair_data(1, i, t))
    end do

    ! do power transformation on precip vector (AWW: consider alternate transforms)
    call normalize_y (4.0d0, y)    ! SHOULD NOT BE HARDWIRED


    !read from nwp predictor file list
    read(unit=55,fmt='(A)',iostat=error) nwp_timestep_file
    print *,trim(nwp_timestep_file)
    if(error /= 0) then
      print *, 'Error reading from NWP predictor file list: ',trim(nwp_input_list)
      stop
    end if

    !read NWP predictor file for current timestep
    call read_nwp(times(t),nwp_timestep_file,nPredict,n_nwp,nwp_vars,nearestGridpoint,x,z,error)
    if(error /= 0) then
      print *, 'Error in read_nwp'
      stop
    end if


    ! -------- loop through all grid cells for a given time step --------
    do g = 1, ngrid, 1
      ! call system_clock(tg1,count_rate)

      deallocate (twx_red)
      deallocate (tx_red)
      allocate (twx_red(nPredict, sta_limit))! these have dim1 = 6
      allocate (tx_red(nPredict, sta_limit))

      ! IF the elevation is valid for this grid cell
      ! (this starts a long section working first on precip, then temp)
      if (z(g, 4) .gt.-200) then
        ! call system_clock(t1,count_rate)

        ! want to reset weights for closest sta_limit stations...
        ! recalc calc_distance_weight function for selected stations
        ! set max_distance equal to the farthest station distance

        ! ---- first, PRECIP ----

        ! set data count integers and initialize reduced arrays to zero
        ndata = 0
        nodata = 0
        w_pcp_red = 0.0
        y_red = 0.0
        x_red = 0.0
        yp_red = 0

        max_distance = 0.0d0
        do i = 1, (close_count(g)-1)
          if (close_meta(5, g, i) .gt. max_distance) then
            max_distance = close_meta (5, g, i)
          end if
        end do

        if (max_distance .le. maxdistance) then
          expand_dist (g) = max_distance
          expand_flag (g) = 0
          max_distance = maxdistance
        else
          max_distance = max_distance + 1.0d0
          expand_flag (g) = 1
          expand_dist (g) = max_distance
        end if

        ! reduced matrices for precip
        slope_flag_pcp = 0
        do i = 1, (close_count(g)-1)
          call calc_distance_weight (max_distance, close_meta(1, g, i), close_meta(2, g, i), &
         & close_meta(3, g, i), close_meta(4, g, i), tmp_weight)

          w_pcp_red (i, i) = tmp_weight
          w_pcp_1d (i) = tmp_weight
          w_pcp_1d_loc (i) = close_loc (g, i)
          y_red (i) = y (close_loc(g, i))
          x_red (i, :) = x (close_loc(g, i), :)

          if (prcp_data(close_loc(g, i), t) .gt. 0.0) then
            ndata = ndata + 1    ! count data points with non-zero precip
            yp_red (i) = 1
          else
            nodata = nodata + 1  ! count data points with zero precip
          end if
        end do

        call normalize_xv (y_red, w_pcp_1d, step_mean, step_std, step_std_all, step_min, step_max, &
       & yp_red)

        y_mean (g, t) = step_mean
        y_std (g, t) = step_std
        y_std_all (g, t) = step_std_all
        y_min (g, t) = step_min
        y_max (g, t) = step_max

        ! ---- second, TEMPERATURES ----

        ! start with initializations
        ndata_t = 0
        nodata_t = 0
        w_temp_red = 0.0
        y_tmean_red = 0.0
        y_trange_red = 0.0
        x_red_t = 0.0

        ! max_distance_t = maxval(close_meta_t(5,g,:))
        max_distance_t = 0.0d0
        do i = 1, (close_count_t(g)-1)
          if (close_meta_t(5, g, i) .gt. max_distance_t) then
            max_distance_t = close_meta_t (5, g, i)
          end if
        end do

        if (max_distance_t .le. maxdistance) then
          max_distance_t = maxdistance
        else
          max_distance_t = max_distance_t + 1.0d0
          expand_flag_t (g) = 1
          expand_dist_t (g) = max_distance_t
        end if

        ! reduced matrices for temperature
        slope_flag_temp = 0
        do i = 1, (close_count_t(g)-1)
          call calc_distance_weight (max_distance_t, close_meta_t(1, g, i), close_meta_t(2, g, i), &
         & close_meta_t(3, g, i), close_meta_t(4, g, i), tmp_weight)

          w_temp_red (i, i) = tmp_weight
          w_temp_1d (i) = tmp_weight
          y_tmean_red (i) = y_tmean (close_loc_t(g, i))
          y_trange_red (i) = y_trange (close_loc_t(g, i))
          x_red_t (i, :) = x (close_loc_t(g, i), :)

          if (y_tmean(close_loc_t(g, i)) .gt.-100.0) then
            ndata_t = ndata_t + 1    ! count data points with valid temperature
          else
            nodata_t = nodata_t + 1  ! count data points with invalid temperature
          end if
        end do

        ! ---- checks on station availability for precip and temp

        if (ndata == 0 .and. nodata == 0) then
          print *, "WARNING:  No stations with data within max distance of grid cell!"
          ! this should not happen if station data are filled
          pop (g, t) = 0.0
          pcp (g, t) = 0.0
          pcperr (g, t) = 0.0
          pop_2 (g, t) = 0.0
          pcp_2 (g, t) = 0.0
          pcperr_2 (g, t) = 0.0
        end if

        ! added AJN Sept 2013
        if (ndata_t == 0 .and. nodata_t == 0) then
          if (t .gt. 1) then
            tmean (g, t) = tmean (g, t-1)
            trange (g, t) = trange (g, t-1)
            tmean_err (g, t) = tmean_err (g, t-1)
            trange_err (g, t) = trange_err (g, t-1)

            tmean_2 (g, t) = tmean_2 (g, t-1)
            trange_2 (g, t) = trange_2 (g, t-1)
            tmean_err_2 (g, t) = tmean_err_2 (g, t-1)
            trange_err_2 (g, t) = trange_err_2 (g, t-1)
          else
            tmean (g, t) = -999
            trange (g, t) = -999
            tmean_err (g, t) = 0.0
            trange_err (g, t) = 0.0

            tmean_2 (g, t) = -999
            trange_2 (g, t) = -999
            tmean_err_2 (g, t) = 0.0
            trange_err_2 (g, t) = 0.0
          end if
        end if
        ! print *,ndata

        ! ========= Precip & temp are processed sequentially, again ===========
        ! this is the start of the PRECIP processing block ---

        if (ndata >= 1) then  ! at least one station close by has pcp > 0

          ! original call
          ! TWX = matmul(TX, w_pcp)

          ! tmp needs to be matmul(TX,X) where TX = TWX_red and X = X_red
          twx_red = matmul (transpose(x_red), w_pcp_red)
          tmp = matmul (twx_red, x_red)
          vv = maxval (abs(tmp), dim=2)

          if (any(vv == 0.0)) then
            slope_flag_pcp = 0
          else
            slope_flag_pcp = 1
          end if

          ! -------------- 1. CALCULATING POP -----------------
          if (nodata == 0) then
            ! print *, "All stations have precip>0, POP = 1.0"
            pop (g, t) = 1.0
            pop_2 (g, t) = 1.0
          else
            ! some stations don't have precip > 0
            if (slope_flag_pcp .eq. 0) then
              pop (g, t) = -999.   ! when not using slope regressions
            else
              ! --- regression with slope ---
              tx_red = transpose (x_red)
              twx_red = matmul (tx_red, w_pcp_red)
              call logistic_regression (x_red, y_red, twx_red, yp_red, b) ! AJN

              !pop (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, :), b))), kind(sp))
              if(-dot_product(Z(g,:),B) < 25.) then
                pop (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, :), b))), kind(sp))
              else
                POP(g,t) = 0.0
              endif

              deallocate (b)
            end if

            ! --- regression without slope ---
            ! AWW note that these now use the 2nd set of T* variables (different dimension)

            deallocate(tx_red_2)   ! just testing
            deallocate(twx_red_2)
            allocate(tx_red_2(nPredict-2, sta_limit))
            allocate(twx_red_2(nPredict-2, sta_limit))

!            tx_red_2(:,1:4) = x_red(:, 1:4)
!            tx_red_2(:,nBase+1:nPredict) = x_red(:,nBase+1:nPredict)
            tx_red_2 = x_red(:,noSlopePredicts)
            tx_red_2 = transpose(tx_red_2)
            twx_red_2 = matmul(tx_red_2, w_pcp_red)
!print *,'npredict ',nPredict,nBase,n_nwp,' indices ',noSlopePredicts
!print *,x_red(:, noSlopePredicts)
!print *,size(x_red(:, noSlopePredicts)),size(y_red),size(twx_red_2),size(yp_red)
            call logistic_regression (x_red(:, noSlopePredicts), y_red, twx_red_2, yp_red, b)!AJN
!print *,'b ',b
            !pop_2 (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, 1:4), b))), kind(sp))
            if(-dot_product(Z(g,noSlopePredicts),B) < 25.) then
!print *,'z ',z(g,noSlopePredicts)
              pop_2 (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, noSlopePredicts), b))), kind(sp))
            else
              POP(g,t) = 0.0
            endif

            deallocate (b)! B must get allocated in logistic reg.; could this also be allocated just once?
          end if
          ! print *, "POP: ", POP(g,t)

          ! -------------- 2. NOW CALCULATING PCP -----------------

          deallocate(twx_red)
          deallocate(tx_red)   ! just testing
          allocate(twx_red(nPredict, sta_limit))
          allocate(tx_red(nPredict, sta_limit))

          if(slope_flag_pcp .eq. 0) then
            pcp (g, t) = -999.
          else
            ! regression with slope
            tx_red = transpose (x_red)
            twx_red = matmul (tx_red, w_pcp_red)

            call least_squares (x_red, y_red, twx_red, b)
            pcp (g, t) = real (dot_product(z(g, :), b), kind(sp))
            deallocate (b)  !AWW-seems to be missing

            wgtsum = 0.0
            errsum = 0.0
            ss_tot = 0.0
            ss_res = 0.0
            do i = 1, (close_count(g)-1), 1
              wgtsum = wgtsum + w_pcp_red (i, i)
              errsum = errsum + (w_pcp_red(i, i)*(pcp(g, t)-y_red(i))**2)
            end do

            pcperr (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))
          end if

          ! regression without slope

          deallocate(tx_red_2)   ! just testing
          deallocate(twx_red_2)
          allocate(tx_red_2(nPredict-2, sta_limit))
          allocate(twx_red_2(nPredict-2, sta_limit))

          ! AWW note that these use the 2nd set of T* variables (different dimension)
          tx_red_2 = transpose (x_red(:, noSlopePredicts))
          twx_red_2 = matmul (tx_red_2, w_pcp_red)
          call least_squares (x_red(:, noSlopePredicts), y_red, twx_red_2, b)

          pcp_2 (g, t) = real (dot_product(z(g, noSlopePredicts), b), kind(sp))

          wgtsum = 0.0
          errsum = 0.0
          ss_tot = 0.0
          ss_res = 0.0
          do i = 1, (close_count(g)-1), 1
            wgtsum = wgtsum + w_pcp_red (i, i)
            errsum = errsum + (w_pcp_red(i, i)*(pcp_2(g, t)-y_red(i))**2)
          end do

          pcperr_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))

          deallocate (b)
          ! print *,'done precip'

        else
          ! this means ndata = 0 for this grid cell and timestep
          ! print *, "INFO:  No stations nearby have pcp > 0, so precip for this cell being set to zero"
          pop (g, t) = 0.0
          pcp (g, t) = 0.0
          pcperr (g, t) = 0.0
          pop_2 (g, t) = 0.0
          pcp_2 (g, t) = 0.0
          pcperr_2 (g, t) = 0.0

        end if ! done with precip if (ndata>=1) block

        ! added AJN Sept 2013
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Temperature OLS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (ndata_t .ge. 1) then ! AJN

          ! regression with slope
          ! AWW note that these use the 1st set of T* variables (6 dim)
          tx_red = transpose (x_red_t)
          twx_red = matmul (tx_red, w_temp_red)
          call least_squares (x_red_t, y_tmean_red, twx_red, b)

          tmean (g, t) = real (dot_product(z(g, :), b), kind(sp))

          errsum = 0.0
          wgtsum = 0.0
          do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
            errsum = errsum + (w_temp_red(i, i)*(tmean(g, t)-y_tmean_red(i))**2)
          end do
          tmean_err (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))
          deallocate (b)

          ! regression without slope

          deallocate(tx_red_2)   ! just testing
          deallocate(twx_red_2)
          allocate(tx_red_2(nPredict-2, sta_limit))
          allocate(twx_red_2(nPredict-2, sta_limit))

          ! AWW note that these use the 2nd set of T* variables
          tx_red_2 = transpose (x_red_t(:, noSlopePredicts))
          twx_red_2 = matmul (tx_red_2, w_temp_red)
          call least_squares (x_red_t(:, noSlopePredicts), y_tmean_red, twx_red_2, b)

          tmean_2 (g, t) = real (dot_product(z(g, noSlopePredicts), b), kind(sp))

          errsum = 0.0
          wgtsum = 0.0
          do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
            errsum = errsum + (w_temp_red(i, i)*(tmean_2(g, t)-y_tmean_red(i))**2)
          end do
          tmean_err_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))
          deallocate (b)

          ! ===== NOW do TRANGE ============

          ! regression with slope

          deallocate(tx_red)   ! just testing
          deallocate(twx_red)
          allocate(tx_red(nPredict, sta_limit))
          allocate(twx_red(nPredict, sta_limit))

          ! AWW note that these use the 1st set of T* variables
          tx_red = transpose (x_red_t)
          twx_red = matmul (tx_red, w_temp_red)
          call least_squares (x_red_t, y_trange_red, twx_red, b)

          trange (g, t) = real (dot_product(z(g, :), b), kind(sp))

          errsum = 0.0
          wgtsum = 0.0
          do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
            errsum = errsum + (w_temp_red(i, i)*(trange(g, t)-y_trange_red(i))**2)
          end do
          trange_err (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))
          deallocate (b)

          ! --- regression without slope ---

          deallocate(tx_red_2)   ! just testing
          deallocate(twx_red_2)
          allocate(tx_red_2(nPredict-2, sta_limit))
          allocate(twx_red_2(nPredict-2, sta_limit))

          ! AWW note that these use the 2nd set of T* variables
          tx_red_2 = transpose (x_red_t(:, noSlopePredicts))
          twx_red_2 = matmul (tx_red_2, w_temp_red)
          call least_squares (x_red_t(:, noSlopePredicts), y_trange_red, twx_red_2, b)

          trange_2 (g, t) = real (dot_product(z(g, noSlopePredicts), b), kind(sp))

          errsum = 0.0
          wgtsum = 0.0
          do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
            sta_temp = real (dot_product(x_red_t(i, noSlopePredicts), b), kind(sp))
            errsum = errsum + (w_temp_red(i, i)*(trange_2(g, t)-y_trange_red(i))**2)
          end do
          trange_err_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))
          deallocate (b)  !AWW-seems to be missing

        else ! alternative to having (ndata_t <= 1)

          ! if not enough stations with data
          ! just use value from previous grid point for now AJN
          print *, 'WARNING:  not enough data stations for current point for temperature; using las&
         &t grid point'

          if (g .gt. 1) then
            trange (g, t) = trange (g-1, t)
            trange_err (g, t) = trange_err (g-1, t)
            tmean (g, t) = tmean (g-1, t)
            tmean_err (g, t) = tmean_err (g-1, t)

            trange_2 (g, t) = trange_2 (g-1, t)
            trange_err_2 (g, t) = trange_err_2 (g-1, t)
            tmean_2 (g, t) = tmean_2 (g-1, t)
            tmean_err_2 (g, t) = tmean_err_2 (g-1, t)
          else
            trange (g, t) = trange (g, t-1)
            trange_err (g, t) = trange_err (g-1, t-1)
            tmean (g, t) = tmean (g, t-1)
            tmean_err (g, t) = tmean_err (g, t-1)

            trange_2 (g, t) = trange_2 (g, t-1)
            trange_err_2 (g, t) = trange_err_2 (g-1, t-1)
            tmean_2 (g, t) = tmean_2 (g, t-1)
            tmean_err_2 (g, t) = tmean_err_2 (g, t-1)
          end if
        end if ! end data check if statement for temperature

      end if ! end check for valid elevation

    end do ! end grid loop

    call system_clock (tg2, count_rate)
    print *, 'Elapsed time for one time step: ', real (tg2-tg1) / real (count_rate)

  end do ! end time record loop

  ! AWW -- just deallocate once at end of subroutine
  deallocate (twx_red)
  deallocate (tx_red)
  deallocate (twx_red_2)
  deallocate (tx_red_2)

end subroutine estimate_forcing_regression
