! Main routine for processing station and grid predictor data and calculating the spatial regression

subroutine estimate_forcing_regression (nTotPredictors, gen_sta_weights, sta_weight_name, nwp_input_list, &
  & nDynPRedictors, nwp_vars, nwp_prcp_var, x, z, ngrid, maxdistance, times, st_rec, end_rec, &
  & stnid, stnvar, directory, sta_limit, kfold_trials, pcp, pop, pcperr, obs_max_pcp, tmean, &
  & tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, error, &
  & pcp_2, pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2, use_stn_weights)

  ! ==============================================================================================
  ! This routine is called during MODE 2 usage:  creates gridded ensembles from station/point data
  ! ==============================================================================================

  use string_mod
  use utim
  use combination_routines
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
      ! inputs
      character(len=500),intent(in) :: sta_weight_name     ! name of station weight binary file
      integer(I4B), intent(in)      :: ngrid               ! number of grid points
      integer(I4B), intent(in)      :: nstns               ! number of stations
      real(DP), intent(in)          :: Z(:,:)              ! grid metadata array
      real(DP), intent(in)          :: X(:,:)              ! station metadata array
      real(DP), intent(in)          :: search_distance     ! default station search distance
      integer(I4B), intent(in)      :: sta_limit           ! maximum number of stations used in regression
      real(DP), intent(in)          :: sta_data(:,:)       ! station data values for precipitation
      real(DP), intent(in)          :: tair_data(:,:,:)    ! station air temperature data
      !in/out
      real(DP), intent(inout)     :: close_meta(:,:,:)     ! 
      real(DP), intent(inout)     :: close_meta_t(:,:,:)
      integer(I4B), intent(inout) :: close_loc(:,:)        ! indices of nearest neighbors for pcp, dim (ngrid, sta_limit)
      integer(I4B), intent(inout) :: close_loc_t(:,:)      ! indices of nearest neighbors for pcp, dim (ngrid, sta_limit)
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
      ! inputs
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

    subroutine read_nwp(currentTime,nwp_timestep_file,nTotPredictors,nDynPRedictors,nwp_vars,station_grid,x,z,error)
      use string_mod
      use utim
      use type
      character (len=2000), intent (in) :: nwp_timestep_file
      character (len=100), intent (in)  :: nwp_vars(:)
      integer(i4b), intent(in)          :: nTotPredictors     !total number of predictors
      integer(i4b), intent(in)          :: nDynPRedictors        !number of NWP predictors
      integer(i4b), intent(in)          :: station_grid(:)        !nearest grid point for every station location
      real(dp), intent(in)          :: currentTime  !current timestep unix time
      real (dp), intent (inout) :: x(:,:), z(:,:) !station and grid predictor matrices
      integer, intent (inout) :: error !error integer
    end subroutine read_nwp

    subroutine station_grid_correspondence(X,Z,close_weights,close_weights_t,close_loc,close_loc_t,&
                                           & nSta,nearestGridpoint)
      use type
      real(dp), intent(in)        :: X(:,:)
      real(dp), intent(in)        :: Z(:,:)
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

    !subroutine max_x (x, smax)
    !  use type
    !  real (dp), intent (in) :: x (:)
    !  real (dp), intent (out) :: smax
    !end subroutine max_x
    
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

    subroutine kfold_crossval(X, Y, W, kfold_trials, kfold_nsamp, n_train, max_n_test, xval_combinations, varUncert)
      use type
      implicit none
      !inputs/outputs 
      real(dp), intent(in)      :: X(:,:)                  ! full station attribute array
      real(dp), intent(in)      :: Y(:)                    ! full station variable value vector  
      real(dp), intent(in)      :: W(:,:)                  ! full station diagonal weight matrix
      integer(I4B), intent(in)  :: kfold_trials            ! number of kfold xval trials
      integer(I4B), intent(in)  :: kfold_nsamp             ! number of samples to draw from
      integer(I4B), intent(in)  :: n_train                 ! number of stations in training sample
      integer(i4b), intent(in)  :: max_n_test              ! maximum test sample size over all trials
      integer(I4B), intent(in)  :: xval_combinations(:,:)  ! array of sampling combinations (integer array indicies)
      real(sp), intent(out)     :: varUncert               ! output: uncertainty estimate from kfold trials
    end subroutine kfold_crossval
    
    !subroutine comb()? --- in a module!
  end interface
  ! =========== end interfaces, start code =============

  real (dp), intent (inout)  :: x (:, :), z (:, :)  ! station and grid point description arrays
  real (dp), intent (in)     :: maxdistance         ! max distance for weight function
  integer (i4b), intent (in) :: ngrid               ! number of grid points
  real (dp), intent (in)     :: times (:)           ! time step array

  ! AWW added next
  integer (i4b), intent (in) :: st_rec, end_rec

  character (len=100), intent (in)  :: stnid (:)           ! station id array
  character (len=100), intent (in)  :: stnvar              ! control file variables
  character (len=500), intent (in)  :: directory
  character (len=500), intent(in)   :: gen_sta_weights     ! flag for generating station weight file
  character (len=500), intent(in)   :: sta_weight_name     ! station weight file name
  character (len=500), intent(in)   :: use_stn_weights     ! flag for doing distance weighted regression
  character (len=2000),intent(in)   :: nwp_input_list      ! file containing list of NWP predictor input files
  character (len=100), intent(in)   :: nwp_vars(:)         ! list of nwp predictor variables
  character (len=100), intent(in)   :: nwp_prcp_var        ! name of nwp predictor variable for precipitation

  character (len=2000)              ::  nwp_timestep_file  ! name of a NWP predictor file

  integer(I4B), intent(inout)       :: nTotPredictors      ! number of total predictors
  integer(I4B), intent(in)          :: nDynPRedictors      ! number of NWP predictors
  integer(I4B), intent(in)          :: sta_limit           ! total number of stations in regression sample size
  integer(I4B), intent(in)          :: kfold_trials        ! number of kfold xval trials

  real (sp), allocatable, intent (out) :: pcp (:, :), pop (:, :), pcperr (:, :) ! output variables for precipitation
  real (sp), allocatable, intent (out) :: tmean (:, :), tmean_err (:, :)        ! OLS tmean estimate and error
  real (sp), allocatable, intent (out) :: trange (:, :), trange_err (:, :)      ! OLS trange estimate and error

  real (sp), allocatable, intent (out) :: tmean_2 (:, :), tmean_err_2 (:, :)    ! OLS tmean estimate and error
  real (sp), allocatable, intent (out) :: trange_2 (:, :), trange_err_2 (:, :)  ! OLS trange estimate and error
  real (sp), allocatable, intent (out) :: pcp_2 (:, :), pop_2 (:, :), pcperr_2 (:, :)

  integer, intent (out)   :: error              ! integer error flag
  real (dp), intent (out) :: mean_autocorr (:)  ! mean autocorrelation from all stations over entire time period
  real (dp), intent (out) :: mean_tp_corr (:)   ! mean correlation for mean temp and precip

  ! vary at each grid point and time step
  real (dp), intent (out) :: obs_max_pcp (:, :) ! max of normalized time step precipitation

  ! Local declarations
  real (dp), allocatable :: y (:), b (:)
  real (dp), allocatable :: y_red (:)           ! reduced matrix for the predictand
  real (dp), allocatable :: x_red (:, :)        ! reduced matrix for the predictors
  real (dp), allocatable :: tmp_x_red (:, :)    ! reduced matrix for extended predictors
  real (dp), allocatable :: x_red_t (:, :)      ! reduced matrix for predictors, transposed (?)
  real (dp), allocatable :: Z_reg(:)            ! final grid point predictor vector

  ! condensed these variables into just 4 that get re-used
  real (dp), allocatable :: twx_red (:, :), tx_red (:, :)        ! reduced matricies (orig)
  real (dp), allocatable :: twx_red_2 (:, :), tx_red_2 (:, :)    ! reduced matricies (dims 2)

  real (dp), allocatable :: w_base (:, :)                        ! initial distance weight matrix
  real (dp), allocatable :: w_pcp_red (:, :), w_temp_red (:, :)  ! reduced distance weight matricies

  real (dp), allocatable :: y_tmean (:), y_trange (:)            ! transformed station data arrays
  real (dp), allocatable :: y_tmean_red (:), y_trange_red (:)    ! transformed station data arrays
  real (dp), allocatable :: stn_prcp (:), prcp_data (:,:), tair_data (:,:,:), stn_tair (:,:)  ! orig stn data arrays
  real (dp), allocatable :: auto_corr (:)      ! lag-1 autocorrelation for stations over entire time period used
  real (dp), allocatable :: t_p_corr (:)       ! correlation between temp and precip
  integer (i4b), allocatable :: yp (:)         ! binary for logistic regression
  integer (i4b), allocatable :: yp_red (:)     ! reduced binary for logistic regression

  logical, allocatable :: stn_miss (:), stn_miss_t (:) ! missing value logical arrays

  real (dp) :: auto_corr_sum, tp_corr_sum

  real (dp) :: step_max                               ! timestep statistics  
  real (dp) :: ss_tot, ss_res                         ! r-squared and variance correction   ! not yet used

  integer (i4b) :: xsize                              ! size of second dimension of input X array
  integer (i4b) :: ntimes, nstns
  integer (i4b) :: t, i, j, g, ndata, nodata, cnt
  integer (i4b) :: ndata_t, nodata_t
  integer (i4b) :: lag, window
  integer (i4b) :: auto_cnt, tp_cnt

  integer(I4B)  :: nBase                              ! number of geophysical predictors
  integer(I4B)  :: prcpPredictInd                     ! index of precipitation predictor in predictor vector
  integer(I4B),allocatable  :: noSlopePredicts(:)
  integer(I4B),allocatable  :: noPrcpPredicts(:)
  integer(I4B),allocatable  :: noPrcpNoSlopePredicts(:)
  logical                   :: drop_precip = .false.
  logical                   :: no_precip = .false.

  ! variables for tracking closest N stations for precipitation
  !integer (i4b), parameter   :: sta_limit = 30
  integer (i4b), allocatable :: close_loc (:, :)
  integer (i4b), allocatable :: close_count (:)
  real (dp), allocatable     :: close_weights (:, :)
  real (dp), allocatable     :: close_meta (:, :, :)
  real (dp)                  :: max_distance
  real (dp), parameter       :: search_distance = 1000.0

  ! variables for tracking closest N stations for temperature
  integer (i4b), allocatable :: close_loc_t (:, :)
  integer (i4b), allocatable :: close_count_t (:)
  real (dp), allocatable     :: close_weights_t (:, :)
  real (dp), allocatable     :: close_meta_t (:, :, :)
  real (dp)                  :: max_distance_t
  
  real (dp)                  :: tmp_weight, wgtsum_pcp, wgtsum_temp, errsum

  integer(I4B), allocatable :: nearestGridpoint(:)

  integer (i4b) :: slope_flag_pcp
  integer (i4b) :: slope_flag_temp

  ! variables to check for singular matrix
  real (dp), allocatable :: tmp (:, :)
  real (dp), allocatable :: vv (:),vv_temp(:)

  ! variables for timing code AJN
  integer (i4b) :: t1, t2, count_rate
  integer (i4b) :: tg1, tg2

  ! variables for keeping track of max_distance modifications
  integer (i4b), allocatable :: expand_flag (:), expand_flag_t (:)
  real (dp), allocatable     :: expand_dist (:), expand_dist_t (:)

  ! kfold cross validation variables 
  integer(i4b), allocatable :: xval_combinations(:,:)         ! array of sampling combination indices
  integer(i4b), allocatable :: sampling_order(:)              ! vector of sampling indices
  integer(i4B)              :: kfold_nsamp                    ! number of samples to draw from !AW same as nstation?
  
  integer(I4B)              :: n_total, n_test, n_train       ! renamed loop limit variables  
  integer(i4b)              :: end_test_stn, start_test_stn   ! endpoints of range of indices in test 
  integer(i4b)              :: trial_n_test, trial_n_train    ! test & train sample size in each trial (varies)
  integer(i4b)              :: max_trial_n_test               ! size of largest test sample in trials
  integer(i4b)              :: nte, ntr                       ! counters for start/end of train/test indices in array
  
  real (dp), parameter      :: transform_exp = 4.0            ! power law transform exponent; must match setting in ens_generation code

  !==============================================================!
  !                     code starts below here                   !
  !==============================================================!

  nstns       = size (stnid)
  ntimes      = size (times)
  xsize       = size (x, 2)
  kfold_nsamp = sta_limit

  ! allocate variable memory
  allocate (y(nstns))
  allocate (y_tmean(nstns), y_trange(nstns))
  allocate (prcp_data(nstns, ntimes))
  allocate (tair_data(2, nstns, ntimes))
  allocate (w_pcp_red(sta_limit, sta_limit))
  allocate (w_temp_red(sta_limit, sta_limit))
  allocate (y_red(sta_limit))
  allocate (y_tmean_red(sta_limit), y_trange_red(sta_limit))
  allocate (x_red(sta_limit, xsize))
  allocate (tmp_x_red(sta_limit, xsize))
  allocate (x_red_t(sta_limit, xsize))
  allocate (tmp(nTotPredictors, nTotPredictors))
  allocate (vv(nTotPredictors))
  allocate (vv_temp(nTotPredictors))
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
  allocate (close_meta(5, ngrid, sta_limit))    ! 
  allocate (close_count(ngrid))

  ! station limit arrays (temp)
  allocate (close_weights_t(ngrid, sta_limit))
  allocate (close_loc_t(ngrid, sta_limit))
  allocate (close_meta_t(5, ngrid, sta_limit))
  allocate (close_count_t(ngrid))

  ! base weight array
  allocate (w_base(ngrid, nstns))

  ! nearest grid index to stations
  allocate (nearestGridpoint(nstns))
 
  ! predictor index array
  allocate (noSlopePredicts(nTotPredictors-2))
  allocate (noPrcpPredicts(nTotPredictors-1))
  allocate (noPrcpNoSlopePredicts(nTotPredictors-3))

  ! max_dist tracking variables
  allocate (expand_dist(ngrid), expand_flag(ngrid))
  allocate (expand_dist_t(ngrid), expand_flag_t(ngrid))

  ! kfold xval variables
  allocate(xval_combinations(kfold_trials, kfold_nsamp+1))
  allocate(sampling_order(kfold_nsamp))

  ! initializations
  pcp = 0.0d0; pop = 0.0d0; pcperr = 0.0d0; auto_corr = 0.0d0
  tmean = 0.0d0; trange = 0.0d0; tmean_err = 0.0d0; trange_err = 0.0d0
  auto_corr_sum = 0.0d0; auto_cnt = 0; tp_corr_sum = 0.0d0; tp_cnt = 0
  w_base = 0.0d0
  expand_dist = 0.0d0; expand_flag = 0; expand_dist_t = 0.0d0; expand_flag_t = 0
  xval_combinations = -999
  sampling_order    = -999

  ! ================= LOOP OVER STATIONS ============
  ! this part calls subroutines that calculate various correlations
  ! can do autocorrelations and correlation between temperature and precipitation
  ! uses an n-day moving average (window) to remove "monthly" cycle from temp
  ! and computes autocorrelation on the anomalies

  do i = 1, nstns, 1

    call read_station (stnvar, stnid(i), directory, st_rec, end_rec, stn_prcp, stn_tair, &
   & stn_miss, stn_miss_t, error)

    print*, "first value check (output period): "
    print*, "   prcp(1): ", stn_prcp(1)
    print*, "   tmin(1): ", stn_tair(1,1)
    print*, "   tmax(1): ", stn_tair(2,1)

    prcp_data (i, :)    = stn_prcp
    tair_data (1, i, :) = stn_tair (1, :)
    tair_data (2, i, :) = stn_tair (2, :)

    ! compute mean autocorrelation for all stations and all times
    lag = 1
    window = 31 ! AWW:  hardwired parameter; should bring out
    call generic_corr (prcp_data(i, :), tair_data(:, i, :), lag, window, auto_corr(i), t_p_corr(i))

    ! check for correlation value outside of -1 to 1
    ! stations with incomplete data are set to -999
    print *, 'station id: ', trim(stnid(i))
    print *, 'auto_corr: ', auto_corr(i)
    print *, 't-p_corr: ', t_p_corr(i)

    if (auto_corr(i) .ge. -1.0 .and. auto_corr(i) .le. 1.0) then
      auto_corr_sum = auto_corr_sum + auto_corr (i)
      auto_cnt = auto_cnt + 1
    end if
    if (t_p_corr(i) .ge. -1.0 .and. t_p_corr(i) .le. 1.0) then
      tp_corr_sum = tp_corr_sum + t_p_corr (i)
      tp_cnt = tp_cnt + 1
    end if

    deallocate (stn_miss_t)  ! must be allocated within read_station
    deallocate (stn_miss)
    deallocate (stn_prcp)
    deallocate (stn_tair)

  end do
  ! =========== end station read loop ============

  error = 0 ! initialize error code returned for function / subroutine calls

  ! output some checks
  print *, 'auto_cnt, tp_cnt=', auto_cnt, tp_cnt
  if (auto_cnt == 0 .or. tp_cnt == 0) then
    print *, 'ERROR:  autocorr or crosscorr (TxP) could not be calculated due to lack of matching pairs'
    stop
  end if
  mean_autocorr = auto_corr_sum / real (auto_cnt, kind(dp))
  mean_tp_corr = tp_corr_sum / real (tp_cnt, kind(dp))

  print *, ' '
  print *, '=========================================================='
  print *, 'Mean correlation stats across all stations'
  print *, '  Temp. lag 1 auto-correlation: ', mean_autocorr(1)
  print *, '  Temp-precip correlation: ', mean_tp_corr(1)
  print *, 'Note: Values may be uninformative for short output periods'
  print *, '=========================================================='
  print *, ' '

  call system_clock (t1, count_rate)

  ! ========= PERFORM REGRESSION ==========================

  ! Create station-grid cell weight matrices before time stepping
  print *, 'Finding nearest stations for each gridpoint and (optionally) generating base weight matrix'

  if(gen_sta_weights .eq. "TRUE" .or. gen_sta_weights .eq. "true") then
    call compute_station_weights(sta_weight_name,ngrid,nstns,X,Z,search_distance, &             ! input
                                 sta_limit,prcp_data,tair_data, &                               ! input
                                 close_meta,close_meta_t,close_loc,close_loc_t, &               ! output
                                 close_count,close_count_t,close_weights,close_weights_t,error) ! output
    if(error /= 0) then
       return
    endif
    call write_station_weights(sta_weight_name, &                                               ! input
                               close_meta,close_meta_t,close_loc,close_loc_t,close_weights,&    ! input
                               close_weights_t,close_count,close_count_t,error)                 ! input
    if(error /= 0) then
       return
    endif
  else
    call read_station_weights(sta_weight_name, &                                                ! input
                              close_meta,close_meta_t,close_loc,close_loc_t,close_weights,&     ! output
                              close_weights_t,close_count,close_count_t,error)                  ! output
    if(error /= 0) then
       return
    endif
  endif

  call system_clock (t2, count_rate)
  print *, 'Elapsed time for weight generation: ', real (t2-t1) / real (count_rate)


  ! create or read in index array of closest grid point to every station
  call station_grid_correspondence(X,Z,close_weights,close_weights_t,close_loc,close_loc_t,&
                                   & nstns,nearestGridpoint)

  nBase = nTotPredictors - nDynPRedictors    ! number of basic static predictors

  if(nDynPRedictors > 0) then
    ! open NWP predictor file list
    open(unit=55,file=trim(nwp_input_list),form='formatted',status='old',iostat=error)
    if( error /= 0) then
      print *,'Error opening NWP predictor file list: ',trim(nwp_input_list)
      stop
    end if
  end if

  ! if using kfold cross-validation, pre-generate subset sampling indices for each fold
  if(kfold_trials > 0) then
    ! Create a station subset sampling indexes for kfold cross-validation
    sampling_order = scrambled_indices(kfold_nsamp)        ! vector
    n_test  = floor(kfold_nsamp*1.0/(kfold_trials*1.0))    ! number of nominal test records per fold, rounded down if not int
    n_train = kfold_nsamp - n_test                         ! number of nominal training records per fold
    print*, ' '; print*, 'Using cross-validation for regression uncertainty estimation'
    print*, 'k_folds = ', kfold_trials
    print*, 'n_total = ', kfold_nsamp
    print*, 'n_train = ', n_train, ' (nominal)'
    print*, 'n_test  = ', n_test, ' (nominal)'

    ! now parcel out test/train samples into array
    max_trial_n_test = 0
    do i = 1, kfold_trials
      ! starting index for each sample
      start_test_stn = (i-1)*n_test + 1

      ! calculate the ending index and size of each sample
      if(i .lt. kfold_trials) then 
        ! all but the last sample have the nominal train/test sizes
        end_test_stn = i * n_test
      else 
        ! if the number of folds is not a multiple of the total station number, the last test sample 
        !    may be larger than the nominal n_test.  
        end_test_stn = kfold_nsamp
      end if  
      
      trial_n_test  = end_test_stn - start_test_stn + 1  
      trial_n_train = kfold_nsamp - trial_n_test

      ! store indices for each trial in an array (rows=trials, cols=indices)
      ! (this isn't vectorized because training subset may be non-contiguous)
      xval_combinations(i,1) = trial_n_test   ! store test sample size as first column
      ntr = 2                                 ! start of training indices:  go in 2:(n_train+1) (fixed)
      nte = trial_n_train + 2                 ! test indices start at n_train+2 
      do j = 1, kfold_nsamp
        if(j .lt. start_test_stn .or. j .gt. end_test_stn) then
          xval_combinations(i,ntr) = sampling_order(j)  ! store training indices in next n_train cols
          ntr = ntr + 1
        else
          xval_combinations(i,nte) = sampling_order(j)  ! store test indices in final cols
          nte = nte + 1
        end if
      end do

    end do  ! end loop over trials
    
    max_trial_n_test  = maxval(xval_combinations(:, 1))
    print*, 'n_test  = ', max_trial_n_test, ' (max)'; print*, ' ' 
    
  end if  ! end IF block for calculating k-fold indices

  ! message about whether station dist weights are used
  if(use_stn_weights .eq. "TRUE" .or. use_stn_weights .eq. "true") then
    print*, 'Using station distance weights'; print*, ' '
  else 
    print*, ' '; print*, 'setting station distance weights equal'; print*, ' '
  end if

  ! =========== LOOP over all TIME steps and populate grids ===============
  print*, '----- Looping over daily timesteps -----'; print*, ' '
  do t = 1, ntimes, 1

    call system_clock (tg1, count_rate)
    print *, "TIME STEP = ", times (t), " (", t, "/", ntimes, ")"

    ! --- assign vectors of station values for prcp, temp, for current time step
    do i = 1, nstns, 1
      y (i) = prcp_data (i, t)
      y_tmean (i) = (tair_data(1, i, t)+tair_data(2, i, t)) / 2.0d0
      y_trange (i) = abs (tair_data(2, i, t)-tair_data(1, i, t))
    end do

    ! do power-law transformation on current time's precipitation values
    call normalize_y (transform_exp, y)

    ! if dynamic predictors are used
    if(nDynPRedictors > 0) then

      ! read from nwp predictor file list
     read(unit=55,fmt='(A)',iostat=error) nwp_timestep_file
      print *,'Adding predictors from: ', trim(nwp_timestep_file)
      if(error /= 0) then
        print *, 'Error reading from NWP predictor file list: ',trim(nwp_input_list)
        stop
      end if

      ! read NWP predictor file for current timestep
      call read_nwp(times(t),nwp_timestep_file,nTotPredictors,nDynPRedictors,nwp_vars,nearestGridpoint,x,z,error)
      if(error /= 0) then
        print *, 'Error in read_nwp'
        stop
      end if
    end if

    ! -------- loop over all grid cells, doing regression for each, for each timestep --------
    do g = 1, ngrid, 1
    
      ! reset various predictor arrangements (use of slope, no slope, no nwp prcp)
      ! AW:  note, does this need to be in the time/grid loop?  
      !      which list to use varies but lists do not (no g or t in settings)
      nTotPredictors = nBase + nDynPRedictors
      noSlopePredicts(1:4)       = (/1,2,3,4/)
      noPrcpPredicts(1:6)        = (/1,2,3,4,5,6/)
      noPrcpNoSlopePredicts(1:4) = (/1,2,3,4/)

      ! if dynamic predictors are used 
      if(nDynPRedictors > 0) then

        do i = 1,nDynPRedictors,1
          noSlopePredicts(nBase-2+i) = nBase+i
        end do

        ! create predictor set without precipitation if it is used
        if(trim(nwp_prcp_var) .eq. "") no_precip = .true.
        if(.not. no_precip) then
          do i = 1,nDynPRedictors,1
            if(nwp_vars(i) .eq. nwp_prcp_var) prcpPredictInd = i+nBase
          end do
          cnt = 1
          do i = 1,nDynPRedictors,1
            if(nwp_vars(i) .ne. nwp_prcp_var) then
              noPrcpPredicts(nBase+cnt) = nBase+i
              noPrcpNoSlopePredicts(nBase-2+cnt) = nBase+i
              cnt = cnt + 1
            end if
          end do
        end if
      end if

      ! call system_clock(tg1,count_rate)

      ! let's reset all memory so we don't have any unforseen memory bugs 
      ! AW: this is inefficient, leads to millions of 10e8 memory calls in big domains - every timestep, every grid
      !     it was removed in the prior SHARP version, TODO to remove again.  Arrays have fixed size, can be allocated once. 
      
      if(allocated(x_red)) then
        deallocate(x_red)
        allocate(x_red(sta_limit, xsize))
      end if

      if(allocated(x_red_t)) then
        deallocate(x_red_t)
        allocate(x_red_t(sta_limit, xsize))
      end if

      if(allocated(twx_red)) then
        deallocate(twx_red)
      allocate (twx_red(nTotPredictors, sta_limit))
      end if

      if(allocated(tx_red)) then
        deallocate(tx_red)
        allocate (tx_red(nTotPredictors, sta_limit))
      end if

      if(allocated(twx_red_2)) then
        deallocate(twx_red_2)
        allocate (twx_red_2(nTotPredictors, sta_limit))
      end if

      if(allocated(twx_red_2)) then
        deallocate(tx_red_2)
        allocate (tx_red_2(nTotPredictors, sta_limit))
      end if

      ! --- IF the elevation is valid for this grid cell
      ! (this starts a long section working first on precip, then temp)
      if (z(g, 4) .gt.-200) then
      
        !create final Z array (grid predictors) for regressions
        allocate(Z_reg(xsize))
        Z_reg = Z(g,:)

        ! call system_clock(t1,count_rate)

        ! need to reset weights for closest sta_limit stations...
        ! recalc calc_distance_weight function for selected stations
        ! set max_distance equal to the farthest station distance

        ! ---- first, PRECIP ----

        ! set data count integers and initialize reduced arrays to zero
        ndata = 0; nodata = 0
        y_red = 0.0; x_red = 0.0; yp_red = 0; w_pcp_red = 0.0

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

        ! calc reduced matrices for precip (loop over nearest stations)
        slope_flag_pcp = 0
        wgtsum_pcp     = 0
        do i = 1, (close_count(g)-1)

          ! AWW 2020 allow for not using station weights
          if(use_stn_weights .eq. "TRUE" .or. use_stn_weights .eq. "true") then
            call calc_distance_weight (max_distance, close_meta(1, g, i), close_meta(2, g, i), &
              & close_meta(3, g, i), close_meta(4, g, i), &
              & tmp_weight)       ! output
            w_pcp_red (i, i) = tmp_weight   ! assign to diagonal in square weight matrix (nstn X nstn)
          else 
            ! weights are equal
            w_pcp_red (i, i) = 1.0/(close_count(g)-1)
          end if
          wgtsum_pcp = wgtsum_pcp + w_pcp_red (i, i)    ! use later for normalization
          !if(t .eq. 1 .and. g .eq. 1) then
          !  print*, 'close loc ',i,' weight ', w_pcp_red(i,i)
          !end if

          y_red (i)    = y(close_loc(g, i))    ! sample of nearby station precip
          x_red (i, :) = x(close_loc(g, i), :)

          if (prcp_data(close_loc(g, i), t) .gt. 0.0) then
            ndata = ndata + 1    ! count data points with non-zero precip
            yp_red (i) = 1
          else
            nodata = nodata + 1  ! count data points with zero precip
          end if
        end do

        ! output max of sampled obs transformed precip for limiting back-transformed precip; propagated to ens_generation step
        !call max_x (y_red, step_max)
        !obs_max_pcp(g, t) = step_max     
        obs_max_pcp(g, t) = max(maxval(y_red), 0.0)

        ! ---- second, TEMPERATURES ----

        ! start with initializations
        ndata_t = 0; nodata_t = 0
        w_temp_red = 0.0; y_tmean_red = 0.0; y_trange_red = 0.0; x_red_t = 0.0

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
        wgtsum_temp     = 0
        do i = 1, (close_count_t(g)-1)

          ! AWW 2020 allow for not using station weights
          if(use_stn_weights .eq. "TRUE" .or. use_stn_weights .eq. "true") then
            call calc_distance_weight (max_distance_t, close_meta_t(1, g, i), close_meta_t(2, g, i), &
              & close_meta_t(3, g, i), close_meta_t(4, g, i), tmp_weight)
            w_temp_red (i, i) = tmp_weight
          else
            w_temp_red (i, i) = 1.0/(close_count_t(g)-1)
          end if
          wgtsum_temp = wgtsum_temp + w_temp_red (i, i)

          y_tmean_red (i) = y_tmean (close_loc_t(g, i))
          y_trange_red (i) = y_trange (close_loc_t(g, i))
          x_red_t (i, :) = x (close_loc_t(g, i), :)

          if (y_tmean(close_loc_t(g, i)) .gt. -100.0) then
            ndata_t = ndata_t + 1    ! count data points with valid temperature
          else
            nodata_t = nodata_t + 1  ! count data points with invalid temperature
          end if
        end do

        ! ---- checks on station availability for precip and temp

        if (ndata == 0 .and. nodata == 0) then
          print *, "WARNING:  No stations with data within max distance of grid cell!"
          ! this should not happen if station data are filled
          pop(g, t)   = 0.0;   pcp(g, t) = 0.0;   pcperr(g, t) = 0.0
          pop_2(g, t) = 0.0; pcp_2(g, t) = 0.0; pcperr_2(g, t) = 0.0
        end if

        ! added AJN Sept 2013
        if (ndata_t == 0 .and. nodata_t == 0) then
          if (t .gt. 1) then
            tmean (g, t)      = tmean (g, t-1)
            trange (g, t)     = trange (g, t-1)
            tmean_err (g, t)  = tmean_err (g, t-1)
            trange_err (g, t) = trange_err (g, t-1)

            tmean_2 (g, t)      = tmean_2 (g, t-1)
            trange_2 (g, t)     = trange_2 (g, t-1)
            tmean_err_2 (g, t)  = tmean_err_2 (g, t-1)
            trange_err_2 (g, t) = trange_err_2 (g, t-1)
          else
            tmean (g, t)      = -999
            trange (g, t)     = -999
            tmean_err (g, t)  = 0.0
            trange_err (g, t) = 0.0

            tmean_2 (g, t)      = -999
            trange_2 (g, t)     = -999
            tmean_err_2 (g, t)  = 0.0
            trange_err_2 (g, t) = 0.0
          end if
        end if

        ! check to see if the station predictor matrix will be well behaved
        !   if not, we need to change the predictor set
        !   concern here is the dynamic predictor:  precipitation=[0]
        !   test precipitation dynamic predictor
        if(.not. no_precip) then
          twx_red = matmul (transpose(x_red), w_pcp_red)
          tmp = matmul (twx_red, x_red)
          vv = maxval (abs(tmp), dim=2)

          twx_red = matmul (transpose(x_red_t), w_temp_red)
          tmp = matmul (twx_red, x_red_t)
          vv_temp = maxval (abs(tmp), dim=2)

          ! full predictor set is badly behaved
          if (any(vv == 0.0) .or. any(vv_temp == 0.0)) then
            ! drop precipitation
            drop_precip = .true.
            ! redefine reduced station predictor arrays
            tmp_x_red = x_red
            deallocate(x_red)
            ! reallocate x_red
            allocate(x_red(sta_limit, xsize-1))
            x_red = tmp_x_red(:,noPrcpPredicts)
            ! x_red_t
            tmp_x_red = x_red_t
            deallocate(x_red_t)
            allocate(x_red_t(sta_limit,xsize-1))
            x_red_t = tmp_x_red(:,noPrcpPredicts)

            ! create final Z array (grid predictors) for regressions
            Z_reg = Z(g,noPrcpPredicts)
            ! update nTotPredictors
            nTotPredictors = nTotPredictors - 1
            ! update noSlopePredicts
            where(noPrcpNoSlopePredicts > prcpPredictInd) noPrcpNoSlopePredicts = noPrcpNoSlopePredicts - 1
            noSlopePredicts = noPrcpNoSlopePredicts
          end if
        end if
        
        ! ========= Precip & temp regressions (done sequentially) ===========
        ! this is the start of the PRECIP processing block ---

        if (ndata >= 1) then  ! at least one station close by has pcp > 0

          ! original call
          ! TWX = matmul(TX, w_pcp)
          ! tmp needs to be matmul(TX, X) where TX = TWX_red and X = X_red

          ! check to see if the station predictor matrix will be well behaved
          !   if not, we need to change the predictor set
          !   concerns are the two slope terms and also zero precip if nwp fields are used
          
          ! test predictor set for slope terms
          twx_red = matmul (transpose(x_red), w_pcp_red)
          tmp     = matmul (twx_red, x_red)
          vv      = maxval (abs(tmp), dim=2)
          ! if badly behaved, drop slope
          if (any(vv == 0.0)) then
            slope_flag_pcp = 0          ! drop slope terms
          else
            slope_flag_pcp = 1
          end if

          ! -------------- 1. CALCULATING POP -----------------
          ! Note:  POP is not currently cross-validated
          
          if (nodata == 0) then
            ! print *, "All stations have precip>0, POP = 1.0"
            pop(g, t)   = 1.0
            pop_2(g, t) = 1.0

          else
            ! some stations have pcp = 0
            if (slope_flag_pcp .eq. 0) then
              pop (g, t) = -999.   ! will not use using slope predictors
            else
              ! --- regression with slope ---
              tx_red  = transpose (x_red)
              twx_red = matmul (tx_red, w_pcp_red)
              call logistic_regression (x_red, y_red, twx_red, yp_red, b)

              if(-dot_product(Z_reg,B) < 25.) then
                pop(g, t) = real (1.0/(1.0+exp(-dot_product(Z_reg, b))), kind(sp))
              else
                POP(g,t)  = 0.0
              end if

              deallocate (b)
            end if              

            ! --- also calculate the regression without slope (decide use in scrf) ---
            ! AWW note that these now use the 2nd set of T* variables (different dimension)

            if(.not. allocated(tx_red_2)) allocate(tx_red_2(nTotPredictors-2, sta_limit))
            if(.not. allocated(tx_red_2)) allocate(twx_red_2(nTotPredictors-2, sta_limit))

            tx_red_2  = x_red(:,noSlopePredicts)
            tx_red_2  = transpose(tx_red_2)
            twx_red_2 = matmul(tx_red_2, w_pcp_red)
            call logistic_regression (x_red(:, noSlopePredicts), y_red, twx_red_2, yp_red, b)
            if(-dot_product(Z_reg(noSlopePredicts),B) < 25.) then
              pop_2(g, t) = real (1.0/(1.0 + exp(-dot_product(Z_reg(noSlopePredicts), b))), kind(sp))
            else
              POP_2(g,t)  = 0.0
            endif

            deallocate (b)! B must get allocated in logistic reg.; could this also be allocated just once?
            
           end if          

          ! -------------- 2. NOW CALCULATING PCP -----------------

          if(pop(g,t) .gt. 0.0) then
            if(allocated(twx_red)) deallocate(twx_red); allocate(twx_red(nTotPredictors, sta_limit))
            if(allocated(tx_red))  deallocate(tx_red);  allocate(tx_red( nTotPredictors, sta_limit))

            if(slope_flag_pcp .eq. 0) then
              pcp(g, t) = -999.
            else
              ! regression with slope
              tx_red  = transpose (x_red)
              twx_red = matmul (tx_red, w_pcp_red)
              call least_squares (x_red, y_red, twx_red, b)        ! solve regression for stations

              pcp(g, t) = real (dot_product(Z_reg, b), kind(sp))   ! apply regression for grids

              ! calculate error/uncertainty either with or without cross-validation
              if(kfold_trials == 0) then

                errsum = 0.0
                do i = 1, (close_count(g)-1), 1
                  ! calculate error for station locations (weighted sample std dev)
                  errsum = errsum + (w_pcp_red(i, i)*(real (dot_product(x_red(i, :), b), kind(sp))-y_red(i))**2)
                end do
                ! normalize and convert error variance into error standard deviation 
                pcperr(g, t) = real((errsum/wgtsum_pcp)**0.5, kind(sp))
                !pcperr(g, t) = real(errsum**0.5, kind(sp))/wgtsum_pcp

              else
                ! cross-validate to calculate error for station locations
                call kfold_crossval(x_red, y_red, w_pcp_red, kfold_trials, &
                                    kfold_nsamp, n_train, max_trial_n_test, xval_combinations, pcperr(g,t)) 
              end if
              deallocate (b) ! b is allocated in least_squares()
            end if

            ! ---- now do regression without slope predictors

            if(allocated(tx_red_2))  deallocate(tx_red_2);  allocate(tx_red_2( nTotPredictors-2, sta_limit))    
            if(allocated(twx_red_2)) deallocate(twx_red_2); allocate(twx_red_2(nTotPredictors-2, sta_limit))
            
            ! AWW note that these use the 2nd set of T* variables (different dimension)
            tx_red_2  = transpose (x_red(:, noSlopePredicts))
            twx_red_2 = matmul (tx_red_2, w_pcp_red)
            call least_squares (x_red(:, noSlopePredicts), y_red, twx_red_2, b)

            pcp_2(g, t) = real (dot_product(Z_reg(noSlopePredicts), b), kind(sp)) 
            
            ! calculate error/uncertainty either with or without cross-validation
            if(kfold_trials == 0) then
              ! don't cross-validate, calc error from weighted sample std dev
              errsum = 0.0
              do i = 1, (close_count(g)-1), 1
                errsum = errsum + (w_pcp_red(i, i)*(real (dot_product(x_red(i, noSlopePredicts), b), &
                                  & kind(sp)) - y_red(i))**2)
              end do
              ! normalize and convert error variance into error standard deviation 
              pcperr_2(g, t) = real((errsum/wgtsum_pcp)**0.5, kind(sp))
              !pcperr_2(g, t) = real(errsum**0.5, kind(sp))/wgtsum_pcp

            else
              ! cross-validate
              call kfold_crossval(x_red(:,noSlopePredicts), y_red, w_pcp_red, kfold_trials, &
                                kfold_nsamp, n_train, max_trial_n_test, xval_combinations, pcperr_2(g,t)) 
            end if
            deallocate (b)   ! b is allocated in least_squares() process

          else
            ! this means pop = 0 for this grid cell and timestep
            pcp (g, t)   = 0.0; pcperr (g, t)   = 0.0
            pcp_2 (g, t) = 0.0; pcperr_2 (g, t) = 0.0
          endif

        else
          ! this means ndata = 0 for this grid cell and timestep
          ! print *, "INFO:  All stations nearby have pcp = 0, so precip for this cell being set to zero"
          pop(g, t)   = 0.0; pcp(g, t)   = 0.0; pcperr(g, t)   = 0.0
          pop_2(g, t) = 0.0; pcp_2(g, t) = 0.0; pcperr_2(g, t) = 0.0

        end if ! done with precip if (ndata>=1) block

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Temperature OLS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (ndata_t .ge. 1) then ! AJN

          ! regression with slope
          ! AWW note that these use the 1st set of T* variables (6 dim)
          tx_red  = transpose (x_red_t)
          twx_red = matmul (tx_red, w_temp_red)
          call least_squares (x_red_t, y_tmean_red, twx_red, b)
          
          tmean(g, t) = real (dot_product(Z_reg, b), kind(sp))

          ! calculate error/uncertainty either with or without cross-validation
          if(kfold_trials == 0) then
            ! don't cross-validate, calc error from weighted sample std dev
            errsum = 0.0
            do i = 1, (close_count_t(g)-1), 1
              errsum = errsum + (w_temp_red(i, i) * (real (dot_product(x_red_t(i, :), b), &
                                & kind(sp)) - y_tmean_red(i))**2)                                 
            end do
            ! normalize and convert error variance into error standard deviation 
            tmean_err(g, t) = real((errsum/wgtsum_temp)**0.5, kind(sp))
            !tmean_err(g, t) = real(errsum**0.5, kind(sp))/wgtsum_temp

          else
            ! cross-validate
            call kfold_crossval(x_red_t, y_tmean_red, w_temp_red, kfold_trials, &
                              kfold_nsamp, n_train, max_trial_n_test, xval_combinations, tmean_err(g,t)) 
          end if
          deallocate (b)

          ! ---- regression without slope predictors
          
          if(.not. allocated(tx_red_2)) allocate(tx_red_2( nTotPredictors-2, sta_limit))
          if(.not. allocated(tx_red_2)) allocate(twx_red_2(nTotPredictors-2, sta_limit))

          ! AWW note that these use the 2nd set of T* variables
          tx_red_2  = transpose (x_red_t(:, noSlopePredicts))
          twx_red_2 = matmul (tx_red_2, w_temp_red)
          call least_squares (x_red_t(:, noSlopePredicts), y_tmean_red, twx_red_2, b)
          tmean_2(g, t) = real (dot_product(Z_reg(noSlopePredicts), b), kind(sp))
          
          ! calculate error/uncertainty either with or without cross-validation
          if(kfold_trials == 0) then
            ! don't cross-validate, calc error from weighted sample std dev
            errsum = 0.0
            do i = 1, (close_count_t(g)-1), 1
              errsum = errsum + (w_temp_red(i, i) * (real (dot_product(x_red_t(i, noSlopePredicts), b), &
                                & kind(sp)) - y_tmean_red(i))**2)                                 
            end do
            ! normalize and convert error variance into error standard deviation 
            tmean_err_2(g, t) = real((errsum/wgtsum_temp)**0.5, kind(sp))
            !tmean_err_2(g, t) = real(errsum**0.5, kind(sp))/wgtsum_temp

          else
            ! cross-validate
            call kfold_crossval(x_red_t(:,noSlopePredicts), y_tmean_red, w_temp_red, kfold_trials, &
                                kfold_nsamp, n_train, max_trial_n_test, xval_combinations, tmean_err_2(g,t)) 
          end if
          deallocate (b)

          ! ===== NOW do TRANGE ============

          ! ---- regression with slope
          if(allocated(tx_red))  deallocate(tx_red);  allocate( tx_red(nTotPredictors, sta_limit))
          if(allocated(twx_red)) deallocate(twx_red); allocate(twx_red(nTotPredictors, sta_limit))
           
          ! note that these use the 1st set of T* predictor arrays
          tx_red  = transpose (x_red_t)
          twx_red = matmul (tx_red, w_temp_red)
          call least_squares (x_red_t, y_trange_red, twx_red, b)

          trange(g, t) = real (dot_product(Z_reg, b), kind(sp))

          ! calculate error/uncertainty either with or without cross-validation
          if(kfold_trials == 0) then
            ! don't cross-validate, calc error from weighted sample std dev
            errsum = 0.0
            do i = 1, (close_count_t(g)-1), 1
              errsum = errsum + (w_temp_red(i, i) * (real (dot_product(x_red_t(i, :), b), &
                                & kind(sp)) - y_trange_red(i))**2)                                 
            end do
            ! normalize and convert error variance into error standard deviation 
            trange_err(g, t) = real((errsum/wgtsum_temp)**0.5, kind(sp))
            !trange_err(g, t) = real(errsum**0.5, kind(sp))/wgtsum_temp

          else
            ! cross-validate
            call kfold_crossval(x_red_t, y_trange_red, w_temp_red, kfold_trials, &
                                kfold_nsamp, n_train, max_trial_n_test, xval_combinations, trange_err(g,t)) 
          end if
          deallocate (b)

          ! ---- regression without slope ---
          if(allocated(tx_red_2))  deallocate(tx_red_2);  allocate(tx_red_2( nTotPredictors-2, sta_limit))
          if(allocated(twx_red_2)) deallocate(twx_red_2); allocate(twx_red_2(nTotPredictors-2, sta_limit))

          ! note that these use the 2nd set of T* variables
          tx_red_2  = transpose (x_red_t(:, noSlopePredicts))
          twx_red_2 = matmul (tx_red_2, w_temp_red)
          call least_squares (x_red_t(:, noSlopePredicts), y_trange_red, twx_red_2, b)

          trange_2 (g, t) = real (dot_product(Z_reg(noSlopePredicts), b), kind(sp))

          ! calculate error/uncertainty either with or without cross-validation
          if(kfold_trials == 0) then
            ! don't cross-validate, calc error from weighted sample std dev
            errsum = 0.0
            do i = 1, (close_count_t(g)-1), 1
              errsum = errsum + (w_temp_red(i, i) * (real (dot_product(x_red_t(i, noSlopePredicts), b), &
                                & kind(sp)) - y_trange_red(i))**2)                                 
            end do
            ! normalize and convert error variance into error standard deviation 
            trange_err_2(g, t) = real((errsum/wgtsum_temp)**0.5, kind(sp))
            !trange_err_2(g, t) = real(errsum**0.5, kind(sp))/wgtsum_temp

          else
            ! cross-validate
            call kfold_crossval(x_red_t(:,noSlopePredicts), y_trange_red, w_temp_red, kfold_trials, &
                              kfold_nsamp, n_train, max_trial_n_test, xval_combinations, trange_err_2(g,t)) 
          end if
          deallocate (b)

        else ! alternative to having (ndata_t <= 1)

          ! if not enough stations with data
          ! just use value from previous grid point for now AJN
          print *, 'WARNING:  not enough data stations for current point for temperature; using previous grid point'

          if (g .gt. 1) then
            trange (g, t)     = trange (g-1, t)
            trange_err (g, t) = trange_err (g-1, t)
            tmean (g, t)      = tmean (g-1, t)
            tmean_err (g, t)  = tmean_err (g-1, t)

            trange_2 (g, t)     = trange_2 (g-1, t)
            trange_err_2 (g, t) = trange_err_2 (g-1, t)
            tmean_2 (g, t)      = tmean_2 (g-1, t)
            tmean_err_2 (g, t)  = tmean_err_2 (g-1, t)
          else
            trange (g, t)     = trange (g, t-1)
            trange_err (g, t) = trange_err (g-1, t-1)
            tmean (g, t)      = tmean (g, t-1)
            tmean_err (g, t)  = tmean_err (g, t-1)

            trange_2 (g, t)     = trange_2 (g, t-1)
            trange_err_2 (g, t) = trange_err_2 (g-1, t-1)
            tmean_2 (g, t)      = tmean_2 (g, t-1)
            tmean_err_2 (g, t)  = tmean_err_2 (g, t-1)
          end if
        end if ! end data check if statement for temperature
        deallocate(Z_reg)

        !reset nTotPredictors if needed
        if(drop_precip) nTotPredictors = nTotPredictors + 1
        !reset drop_precip
        drop_precip = .false.

      end if ! end check for valid elevation

      if(allocated(noSlopePredicts)) then
        deallocate(noSlopePredicts)
        allocate(noSlopePredicts(nTotPredictors-2))
      end if
    end do ! end grid loop

    call system_clock (tg2, count_rate)
    print *, 'Elapsed time for one time step: ', real (tg2-tg1) / real (count_rate)

  end do ! end time record loop

  ! AWW -- ideally just deallocate once at end of subroutine, allocate once at start
  deallocate (twx_red)
  deallocate (tx_red)
  deallocate (twx_red_2)
  deallocate (tx_red_2)

end subroutine estimate_forcing_regression
