! AWW-2016Jan, modifications to handle time subsetting and reduce mem alloc, and clean up
!   renamed from estimate_precip; add also 'directory' var, changed some var names

subroutine estimate_forcing_regression (x, z, nsta, ngrid, maxdistance, times, st_rec, end_rec, &
& stnid, stnvar, site_var, site_var_t, site_list, directory, pcp, pop, pcperr, tmean, tmean_err, &
& trange, trange_err, mean_autocorr, mean_tp_corr, y_mean, y_std, y_std_all, y_min, y_max, error, &
& pcp_2, pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2)

  ! ===============================================================================================
  ! This routine is called during MODE 2 usage:  creates gridded ensembles from station/point data
  ! ===============================================================================================
 
  use strings
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
 
    ! subroutine read_station(stnvar, stnid, site_list, vals, tair_vals, vals_miss, vals_miss_t, error)
    subroutine read_station (stnvar, stnid, site_list, directory, times, st_rec, end_rec, vals, tair_vals, &
   & vals_miss, vals_miss_t, error)
      use type
      character (len=100), intent (in) :: stnvar
      character (len=100), intent (in) :: stnid
      character (len=500), intent (in) :: site_list
      character (len=500), intent (in) :: directory !AWW
      real (dp), intent (in) :: times (:)!AWW
      integer (i4b), intent (in) :: st_rec, end_rec ! AWW
      real (dp), allocatable, intent (out) :: vals (:), tair_vals (:, :)
      logical, allocatable, intent (out) :: vals_miss (:), vals_miss_t (:)
      integer, intent (out) :: error
    end subroutine read_station
 
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
 
    subroutine calc_weights (times, tt, x, w)
      use type
      real (dp), intent (in) :: times (:)
      integer (i4b), intent (in) :: tt
      real (dp), intent (in) :: x (:, :)
      real (dp), allocatable, intent (out) :: w (:, :)
    end subroutine calc_weights
 
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
 
  real (dp), intent (in) :: x (:, :), z (:, :)!station and grid point description arrays
  real (dp), intent (in) :: maxdistance !max distance for weight function
  integer (i4b), intent (in) :: nsta, ngrid !nuber of input stations and grid points
  real (dp), intent (in) :: times (:)!time step array
 
  integer (i4b), intent (in) :: st_rec, end_rec ! AWW
 
  character (len=100), intent (in) :: stnid (:)!station id array
  character (len=100), intent (in) :: stnvar, site_var, site_var_t !control file variables
  character (len=500), intent (in) :: site_list !file name of station list
  character (len=500), intent (in) :: directory !AWw feb-2016

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
 
  real (dp), allocatable :: y (:), twx (:, :), b (:), tx (:, :)
 
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
  real (dp), allocatable :: auto_corr (:)!lag-1 autocorrelation for stations over entire time period used
  real (dp), allocatable :: t_p_corr (:)!correlation between temp and precip
  integer (i4b), allocatable :: yp (:)!binary for logistic regression
  integer (i4b), allocatable :: yp_red (:)!reduced binary for logistic regression
 
  logical, allocatable :: stn_miss (:), stn_miss_t (:)!missing value logical arrays
 
  real (dp) :: m_tmean, m_trange !used to fill missing values in stations
  real (dp) :: errsum, wgtsum, p, sta_pcp, sta_temp
  real (dp) :: auto_corr_sum, tp_corr_sum
  real (dp) :: step_mean, step_std, step_std_all, step_min, step_max !timestep statistics
  real (dp) :: rsqr, ss_tot, ss_res, vif !r-squared and variance correction
 
  integer (i4b) :: xsize !size of second dimension of input X array
  integer (i4b) :: ntimes, nstns
  integer (i4b) :: t, i, g, ndata, nodata
  integer (i4b) :: ndata_t, nodata_t
  integer (i4b) :: lag, window, tc, trc
  integer (i4b) :: auto_cnt, tp_cnt
  integer (i4b) :: stn_count
 
  ! variables for tracking closest N stations for precipitation
  integer (i4b) :: out_loc
  integer (i4b), parameter :: sta_limit = 30
  integer (i4b), allocatable :: close_loc (:, :)
  integer (i4b), allocatable :: close_count (:)
  real (dp), allocatable :: close_weights (:, :)
  real (dp), allocatable :: close_meta (:, :, :)
  real (dp) :: min_weight
  real (dp) :: max_distance
  real (dp), parameter :: search_distance = 1000.0
 
  ! variables for tracking closest N stations for temperature
  integer (i4b) :: out_loc_t
  integer (i4b), allocatable :: close_loc_t (:, :)
  integer (i4b), allocatable :: close_count_t (:)
  real (dp), allocatable :: close_weights_t (:, :)
  real (dp), allocatable :: close_meta_t (:, :, :)
  real (dp) :: min_weight_t
  real (dp) :: max_distance_t
  real (dp) :: tmp_pcp
  real (dp) :: tmp_weight
 
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
  allocate (tmp(6, 6))
  allocate (vv(6))
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
 
  ! max_dist tracking variables
  allocate (expand_dist(ngrid), expand_flag(ngrid))
  allocate (expand_dist_t(ngrid), expand_flag_t(ngrid))
 
  ! initializations
  pcp = 0.0d0
  pop = 0.0d0
  pcperr = 0.0d0
  auto_corr = 0.0d0
  stn_count = 1
 
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
    !print *,'station read'
    !print *,stnvar,stnid(i),site_var,site_var_t,site_list,Times  !AWW modified
    call read_station (stnvar, stnid(i), site_list, directory, times, st_rec, end_rec, stn_prcp, stn_tair, &
   & stn_miss, stn_miss_t, error)
 
    stn_count = stn_count + 1 !AWW doesn't appear to be used
 
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
 
  print *, 'Temp lag-1 autocorrelation: ', mean_autocorr (1)
  print *, 'Temp-precip correlation: ', mean_tp_corr (1)
  print *, '===================================================='
  print *, ' '

  call system_clock (t1, count_rate)
 
  ! ========= LOOP OVER GRID CELLS ==================
  print *, 'Generating base weight matrix and finding nearest stations for each gridpoint'
  ! AJN 01/15/2014 -- pulled this weight generation step out of time loop, which is now down below
 
  do g = 1, ngrid, 1
 
    ! AJN
    close_count (g) = 1   ! these are for precip
    min_weight = 0.0d0
    close_weights (g, :) = 0.0d0  ! close station weights
 
    close_count_t (g) = 1  ! these are for temp
    min_weight_t = 0.0d0
    close_weights_t (g, :) = 0.0d0
 
    ! for current grid cell, loop through stations, find distance
    do i = 1, nstns, 1
      ! setup distinct weight matrices for precip and temperature
      ! x() are station lonlat; z() are grid lonlat; returns weight (w_base) for grd-to-stn
      call calc_distance_weight (search_distance, x(i, 2), x(i, 3), z(g, 2), z(g, 3), w_base(g, i))
 
      ! original call
      ! call calc_distance_weight(maxDistance, X(i,2), X(i,3), Z(g,2), Z(g,3), w_temp(i,i))
 
      ! also set some logic to limit the number of stations to the N closest
      min_weight = 0.0d0

      ! THIS uses the first value in the timeseries to decide if a variable is present for the stn
      if (w_base(g, i) .gt. min_weight .and. prcp_data(i, 1) .gt.-100.0d0) then
 
        if (close_count(g) .le. sta_limit) then
          close_weights (g, close_count(g)) = w_base (g, i)
          close_loc (g, close_count(g)) = i
          close_meta (1, g, close_count(g)) = x (i, 2)
          close_meta (2, g, close_count(g)) = x (i, 3)
          close_meta (3, g, close_count(g)) = z (g, 2)
          close_meta (4, g, close_count(g)) = z (g, 3)
          call calc_distance (x(i, 2), x(i, 3), z(g, 2), z(g, 3), close_meta(5, g, close_count(g)))
          close_count (g) = close_count (g) + 1
        else
          min_weight = minval (close_weights(g, :), 1)
          if (w_base(g, i) .gt. min_weight) then
            out_loc = minloc (close_weights(g, :), 1)
            close_weights (g, out_loc) = w_base (g, i)
            close_loc (g, out_loc) = i
            close_meta (1, g, out_loc) = x (i, 2)
            close_meta (2, g, out_loc) = x (i, 3)
            close_meta (3, g, out_loc) = z (g, 2)
            close_meta (4, g, out_loc) = z (g, 3)
            call calc_distance (x(i, 2), x(i, 3), z(g, 2), z(g, 3), close_meta(5, g, out_loc))
          end if
        end if
      end if
 
      ! need to repeat above for temperature since that data is independent of precipitation
      min_weight_t = 0.0d0
 
      if (w_base(g, i) .gt. min_weight_t .and. tair_data(1, i, 1) .gt.-200.0d0) then
        if (close_count_t(g) .le. sta_limit) then
          close_weights_t (g, close_count_t(g)) = w_base (g, i)
          close_loc_t (g, close_count_t(g)) = i
          close_meta_t (1, g, close_count_t(g)) = x (i, 2)
          close_meta_t (2, g, close_count_t(g)) = x (i, 3)
          close_meta_t (3, g, close_count_t(g)) = z (g, 2)
          close_meta_t (4, g, close_count_t(g)) = z (g, 3)
          call calc_distance (x(i, 2), x(i, 3), z(g, 2), z(g, 3), close_meta_t(5, g, &
         & close_count_t(g)))
          close_count_t (g) = close_count_t (g) + 1
        else
          min_weight_t = minval (close_weights(g, :), 1)
          if (w_base(g, i) .gt. min_weight) then
            out_loc_t = minloc (close_weights_t(g, :), 1)
            close_weights_t (g, out_loc_t) = w_base (g, i)
            close_loc_t (g, out_loc_t) = i
            close_meta_t (1, g, out_loc_t) = x (i, 2)
            close_meta_t (2, g, out_loc_t) = x (i, 3)
            close_meta_t (3, g, out_loc_t) = z (g, 2)
            close_meta_t (4, g, out_loc_t) = z (g, 3)
            call calc_distance (x(i, 2), x(i, 3), z(g, 2), z(g, 3), close_meta(5, g, out_loc_t))
          end if
        end if
      end if
 
    end do ! end station loop
  end do 
  ! ============== end grid point loop to find nearest stations ==============
 
  call system_clock (t2, count_rate)
  print *, 'Elapsed time for weight generation: ', real (t2-t1) / real (count_rate)
 
  ! AWW-Feb2016:  just allocate grids once time, and re-use in code below
  allocate (twx_red(6, sta_limit))! these have dim1 = 6
  allocate (tx_red(6, sta_limit))
  allocate (twx_red_2(4, sta_limit))! these are for no slope calcs, have dim1 = 4
  allocate (tx_red_2(4, sta_limit))
 
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
    call normalize_y (4.0d0, y)
 
    ! -------- loop through all grid cells for a given time step --------
    do g = 1, ngrid, 1
      ! call system_clock(tg1,count_rate)

      deallocate (twx_red)
      deallocate (tx_red)
      allocate (twx_red(6, sta_limit))! these have dim1 = 6
      allocate (tx_red(6, sta_limit))

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
              call logistic_regression (x_red, y_red, twx_red, yp_red, b)!AJN
              pop (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, :), b))), kind(sp))
 
              deallocate (b)
            end if
 
            ! --- regression without slope ---
            ! AWW note that these now use the 2nd set of T* variables (different dimension)

            deallocate(tx_red_2)   ! just testing
            deallocate(twx_red_2)
            allocate(tx_red_2(4, sta_limit))
            allocate(twx_red_2(4, sta_limit))

            tx_red_2 = transpose(x_red(:, 1:4))
            twx_red_2 = matmul(tx_red_2, w_pcp_red)
            call logistic_regression (x_red(:, 1:4), y_red, twx_red_2, yp_red, b)!AJN
            pop_2 (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, 1:4), b))), kind(sp))
 
            deallocate (b)! B must get allocated in logistic reg.; could this also be allocated just once?
          end if
          ! print *, "POP: ", POP(g,t)
 
          ! -------------- 2. NOW CALCULATING PCP -----------------

          deallocate(twx_red)
          deallocate(tx_red)   ! just testing
          allocate(twx_red(6, sta_limit))
          allocate(tx_red(6, sta_limit))

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
          allocate(tx_red_2(4, sta_limit))
          allocate(twx_red_2(4, sta_limit))

          ! AWW note that these use the 2nd set of T* variables (different dimension)
          tx_red_2 = transpose (x_red(:, 1:4))
          twx_red_2 = matmul (tx_red_2, w_pcp_red)
          call least_squares (x_red(:, 1:4), y_red, twx_red_2, b)
 
          pcp_2 (g, t) = real (dot_product(z(g, 1:4), b), kind(sp))
 
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
          allocate(tx_red_2(4, sta_limit))
          allocate(twx_red_2(4, sta_limit))

          ! AWW note that these use the 2nd set of T* variables
          tx_red_2 = transpose (x_red_t(:, 1:4))
          twx_red_2 = matmul (tx_red_2, w_temp_red)
          call least_squares (x_red_t(:, 1:4), y_tmean_red, twx_red_2, b)
 
          tmean_2 (g, t) = real (dot_product(z(g, 1:4), b), kind(sp))
 
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
          allocate(tx_red(6, sta_limit))
          allocate(twx_red(6, sta_limit))

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
          allocate(tx_red_2(4, sta_limit))
          allocate(twx_red_2(4, sta_limit))

          ! AWW note that these use the 2nd set of T* variables
          tx_red_2 = transpose (x_red_t(:, 1:4))
          twx_red_2 = matmul (tx_red_2, w_temp_red)
          call least_squares (x_red_t(:, 1:4), y_trange_red, twx_red_2, b)
 
          trange_2 (g, t) = real (dot_product(z(g, 1:4), b), kind(sp))
 
          errsum = 0.0
          wgtsum = 0.0
          do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
            sta_temp = real (dot_product(x_red_t(i, 1:4), b), kind(sp))
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