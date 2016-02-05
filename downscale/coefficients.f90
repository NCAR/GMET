Subroutine estimate_coefficients (D, nvars, Lats, Lons, Times, stnid, stnlat, stnlon, stnalt, stnvar, site_var, &
& site_list, C, POC, error)
  Use type
  Implicit None
 
  Interface
 
    Subroutine read_station (stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss, &
   & vals_miss_t, error)
      Use type
      Character (Len=100), Intent (In) :: stnvar
      Character (Len=100), Intent (In) :: stnid
      Character (Len=100), Intent (In) :: site_var, site_var_t
      Character (Len=500), Intent (In) :: site_list
      Real (DP), Intent (In) :: Times (:)
      Real (DP), Allocatable, Intent (Out) :: vals (:), tair_vals (:, :)
      Logical, Allocatable, Intent (Out) :: vals_miss (:), vals_miss_t (:)
      Integer, Intent (Out) :: error
    End Subroutine read_station
 
    Subroutine normalize_X (X)
      Use type
      Real (DP), Intent (Inout) :: X (:, :)
    End Subroutine normalize_X
 
    Subroutine normalize_Y (texp, Y)
      Use type
      Real (DP), Intent (In) :: texp !transform exponent
      Real (DP), Intent (Inout) :: Y (:)
    End Subroutine normalize_Y
 
    Subroutine calc_weights (Times, tt, X, W)
      Use type
      Real (DP), Intent (In) :: Times (:)
      Integer (I4B), Intent (In) :: tt
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Allocatable, Intent (Out) :: W (:, :)
    End Subroutine calc_weights
 
    Subroutine least_squares (X, Y, TX, B)
      Use type
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Intent (In) :: Y (:)
      Real (DP), Intent (In) :: TX (:, :)
      Real (DP), Allocatable, Intent (Out) :: B (:)
    End Subroutine least_squares
 
    Subroutine logistic_regressionrf (X, Y, TX, B)
      Use type
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Intent (In) :: Y (:)
      Real (DP), Intent (In) :: TX (:, :)
      Real (DP), Allocatable, Intent (Out) :: B (:)
    End Subroutine logistic_regressionrf
 
  End Interface
 
  Real (DP), Intent (In) :: D (:, :, :), Lats (:), Lons (:)
  Real (DP), Intent (In) :: Times (:)
  Integer (I4B), Intent (In) :: nvars
  Character (Len=100), Intent (In) :: stnid (:)
  Real (DP), Intent (In) :: stnlat (:), stnlon (:), stnalt (:)
  Character (Len=100), Intent (In) :: stnvar
  Character (Len=100), Intent (In) :: site_var
  Character (Len=500), Intent (In) :: site_list
  Real (DP), Allocatable, Intent (Out) :: C (:, :, :), POC (:, :, :)
  Integer, Intent (Out) :: error
 
  Real (DP), Allocatable :: X (:, :), Y (:), XS (:, :), XP (:, :), TWXS (:, :), YS (:), W (:, :), B (:), tair_vals (:, &
 & :)
 
  Real (DP), Allocatable :: TS (:)
  Logical, Allocatable :: Y_miss (:), Y_miss_t (:)
  Integer (I4B) :: ntimes, nlats, nlons, nstns
  Integer (I4B) :: i, j, k, t, v, tt
  Real (DP) :: minlondis, minlatdis
  Integer (I4B) :: minlatk, minlonj, gridindex
  Integer (I4B) :: vshape (2)
  Character (Len=100) :: site_var_t
  Real (DP) :: transform_exp !precipitation transform variable, the exponent of the transform  norm_pcp = pcp^(1/transform_exp)
 
 
  error = 0
  transform_exp = 4.0d0
 
  ntimes = size (Times)
  nlats = size (Lats)
  nlons = size (Lons)
  nstns = size (stnlat)
 
  Allocate (X(ntimes, nvars+1))
  If (trim(stnvar) .Eq. "PRCP") Then
    Allocate (POC(nstns, ntimes, nvars+1))
    POC = 0.0d0
  End If
  Allocate (C(nstns, ntimes, nvars+1))
  C = 0.0d0
 
  Do i = 1, nstns, 1
 
    minlondis = 360.0
    minlonj = - 1
    Do j = 1, nlons, 1
      If (Abs(stnlon(i)-Lons(j)) < minlondis) Then
        minlondis = Abs (stnlon(i)-Lons(j))
        minlonj = j
      End If
    End Do
    minlatdis = 180.0
    minlatk = - 1
    Do k = 1, nlats, 1
      If (Abs(stnlat(i)-Lats(k)) < minlatdis) Then
        minlatdis = Abs (stnlat(i)-Lats(k))
        minlatk = k
      End If
    End Do
 
    If (minlonj ==-1 .Or. minlatk ==-1) Then
      Print *, "Failed to find closest grid point for station: ", trim (stnid(i))
      error = 1
      Return
    End If
 
    gridindex = ((minlatk-1)*nlons) + minlonj
 
    Print *, "Station: ", trim (stnid(i)), stnlat (i), stnlon (i)
    Print *, "Closest Grid point: ", gridindex, Lats (minlatk), Lons (minlonj)
 
    X (:, 1) = 1.0
    Do v = 1, nvars, 1
      X (:, v+1) = D (v, gridindex, :)
      Do t = 1, ntimes, 1
        If (X(t, v+1) /= 0.0) Exit
        If (t == ntimes) Then
          Print *, "ERROR: var ", v, " is all zero for station ", i
          X (:, v+1) = 1.0
        End If
      End Do
    End Do
 
    Call read_station (stnvar, stnid(i), site_var, site_var_t, site_list, Times, Y, tair_vals, Y_miss, Y_miss_t, error)
 
    If (count(Y_miss) < ntimes) Then
      Allocate (TS(count(Y_miss)))
      Allocate (XS(count(Y_miss), nvars+1))
      Allocate (XP(count(Y_miss), nvars+1))
      Allocate (TWXS(nvars+1, count(Y_miss)))
      Allocate (YS(count(Y_miss)))
      TS (:) = pack (Times, Y_miss)
      Do v = 1, nvars + 1, 1
        XS (:, v) = pack (X(:, v), Y_miss)
      End Do
      XP (:, :) = XS (:, :)
      YS (:) = pack (Y, Y_miss)
    Else
      Allocate (TS(ntimes))
      Allocate (XS(ntimes, nvars+1))
      Allocate (XP(ntimes, nvars+1))
      Allocate (TWXS(nvars+1, ntimes))
      Allocate (YS(ntimes))
      TS (:) = Times (:)
      XS (:, :) = X (:, :)
      XP (:, :) = XS (:, :)
      YS (:) = Y (:)
    End If
 
    Call normalize_X (XS)
    Call normalize_Y (transform_exp, YS)
 
    tt = 1
    Do t = 1, ntimes, 1
 
      If (Y_miss(t) .Eqv. .True.) Then
 
        Call calc_weights (TS, tt, XS, W)
        TWXS = matmul (transpose(XS), W)
        If (trim(stnvar) .Eq. "PRCP") Then
 
          Call logistic_regressionrf (XP, YS, TWXS, B)
 
          POC (i, t, :) = B (:)
          Deallocate (B)
        End If
 
        Call least_squares (XP, YS, TWXS, B)
        C (i, t, :) = B (:)
 
        Deallocate (B)
 
        tt = tt + 1
      Else
        If (trim(stnvar) .Eq. "PRCP") Then
          POC (i, t, :) = - 999.99
        End If
        C (i, t, :) = - 999.99
      End If
    End Do
 
    Deallocate (YS)
    Deallocate (TWXS)
    Deallocate (XS)
    Deallocate (TS)
    Deallocate (XP)
 
  End Do
 
End Subroutine estimate_coefficients
 
 
Subroutine estimate_precip (X, Z, nsta, ngrid, maxDistance, Times, stnid, stnvar, site_var, site_var_t, site_list, PCP, &
& POP, PCPERR, tmean, tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, Y_mean, Y_std, y_std_all, y_min, &
& y_max, error, PCP_2, POP_2, PCPERR_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
 
  Use strings
  Use utim
  Use type
  Implicit None
 
  Interface
 
    Subroutine read_station (stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss, &
   & vals_miss_t, error)
      Use type
      Character (Len=100), Intent (In) :: stnvar
      Character (Len=100), Intent (In) :: stnid
      Character (Len=100), Intent (In) :: site_var
      Character (Len=100), Intent (In) :: site_var_t
      Character (Len=500), Intent (In) :: site_list
      Real (DP), Intent (In) :: Times (:)
      Real (DP), Allocatable, Intent (Out) :: vals (:), tair_vals (:, :)
      Logical, Allocatable, Intent (Out) :: vals_miss (:), vals_miss_t (:)
      Integer, Intent (Out) :: error
    End Subroutine read_station
 
    Subroutine normalize_X (X)
      Use type
      Real (DP), Intent (Inout) :: X (:, :)
    End Subroutine normalize_X
 
    Subroutine normalize_Xv (X, weight, mean, stdev, stdev_all, smin, smax, Yp)
      Use type
      Real (DP), Intent (Inout) :: X (:)
      Real (DP), Intent (In) :: weight (:)
      Real (DP), Intent (Out) :: mean
      Real (DP), Intent (Out) :: stdev
      Real (DP), Intent (Out) :: stdev_all
      Real (DP), Intent (Out) :: smin
      Real (DP), Intent (Out) :: smax
      Integer (I4B), Intent (Out) :: Yp (:)
    End Subroutine normalize_Xv
 
    Subroutine normalize_Y (texp, Y)
      Use type
      Real (DP), Intent (In) :: texp !transform exponent
      Real (DP), Intent (Inout) :: Y (:)
    End Subroutine normalize_Y
 
    Subroutine calc_weights (Times, tt, X, W)
      Use type
      Real (DP), Intent (In) :: Times (:)
      Integer (I4B), Intent (In) :: tt
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Allocatable, Intent (Out) :: W (:, :)
    End Subroutine calc_weights
 
    Subroutine least_squares (X, Y, TX, B)
      Use type
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Intent (In) :: Y (:)
      Real (DP), Intent (In) :: TX (:, :)
      Real (DP), Allocatable, Intent (Out) :: B (:)
    End Subroutine least_squares
 
    Subroutine logistic_regression (X, Y, TX, Yp, B)
      Use type
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Intent (In) :: Y (:)
      Real (DP), Intent (In) :: TX (:, :)
      Integer (I4B), Intent (In) :: Yp (:)
      Real (DP), Allocatable, Intent (Out) :: B (:)
    End Subroutine logistic_regression
 
    Subroutine generic_corr (stn_data, tair_data, lag, window, auto_corr, t_p_corr)
      Use type
 
      Real (DP), Intent (In) :: stn_data (:)
      Real (DP), Intent (In) :: tair_data (:, :)
      Integer (I4B), Intent (In) :: lag
      Integer (I4B), Intent (In) :: window
      Real (DP), Intent (Out) :: auto_corr
      Real (DP), Intent (Out) :: t_p_corr
 
    End Subroutine generic_corr
 
    Subroutine calc_distance (lat1, lon1, lat2, lon2, dist)
      Use type
      Implicit None
 
      Real (DP), Intent (In) :: lat1, lon1, lat2, lon2
      Real (DP), Intent (Out) :: dist
 
    End Subroutine calc_distance
 
 
    Subroutine heapsort (n, ra, rn)
      Use type
      Implicit None
 
      Integer (I4B), Intent (In) :: n
      Integer (I4B), Dimension (:), Intent (Inout) :: rn
      Real (DP), Dimension (:), Intent (Inout) :: ra
 
    End Subroutine heapsort
 
  End Interface
 
  Real (DP), Intent (In) :: X (:, :), Z (:, :)!station and grid point description arrays
  Real (DP), Intent (In) :: maxDistance !max distance for weight function
  Integer (I4B), Intent (In) :: nsta, ngrid !nuber of input stations and grid points
  Real (DP), Intent (In) :: Times (:)!time step array
  Character (Len=100), Intent (In) :: stnid (:)!station id array
  Character (Len=100), Intent (In) :: stnvar, site_var, site_var_t !control file variables
  Character (Len=500), Intent (In) :: site_list !file name of station list
  Real (SP), Allocatable, Intent (Out) :: PCP (:, :), POP (:, :), PCPERR (:, :)!output variables for precipitation
  Real (SP), Allocatable, Intent (Out) :: tmean (:, :), tmean_err (:, :)!OLS tmean estimate and error
  Real (SP), Allocatable, Intent (Out) :: trange (:, :), trange_err (:, :)!OLS trange estimate and error
 
  Real (SP), Allocatable, Intent (Out) :: tmean_2 (:, :), tmean_err_2 (:, :)!OLS tmean estimate and error
  Real (SP), Allocatable, Intent (Out) :: trange_2 (:, :), trange_err_2 (:, :)!OLS trange estimate and error
  Real (SP), Allocatable, Intent (Out) :: PCP_2 (:, :), POP_2 (:, :), PCPERR_2 (:, :)
 
 
  Integer, Intent (Out) :: error !integer error flag
  Real (DP), Intent (Out) :: mean_autocorr (:)!mean autocorrelation from all stations over entire time period
  Real (DP), Intent (Out) :: mean_tp_corr (:)!mean correlation for mean temp and precip
 
   !vary at each grid point and time step
  Real (DP), Intent (Out) :: Y_mean (:, :), Y_std (:, :)!std and mean of time step precipitation
  Real (DP), Intent (Out) :: y_std_all (:, :)!std of time step precip including stations with zero precip
  Real (DP), Intent (Out) :: y_min (:, :), y_max (:, :)!min & max  of normalized time step precipitation
 
  Real (DP), Allocatable :: Y (:), TWX (:, :), B (:), TX (:, :)
 
  Real (DP), Allocatable :: Y_red (:), TWX_red (:, :), TX_red (:, :)!reduced matricies
  Real (DP), Allocatable :: X_red (:, :)!reduced matricies
 
  Real (DP), Allocatable :: TWX_red_t (:, :), TX_red_t (:, :)!reduced matricies
  Real (DP), Allocatable :: X_red_t (:, :)!reduced matricies
 
  Real (DP), Allocatable :: TWX_red_tr (:, :), TX_red_tr (:, :)!reduced matricies
  Real (DP), Allocatable :: X_red_tr (:, :)!reduced matricies
 
  Real (DP), Allocatable :: w_base (:, :)!initial distance weight matrix
  Real (DP), Allocatable :: w_pcp_1d (:), w_temp_1d (:)
  Integer (I4B), Allocatable :: w_pcp_1d_loc (:), w_temp_1d_loc (:)
!  real(DP), allocatable :: w_pcp(:,:), w_temp(:,:) !distance weight matrices for a specific grid point
  Real (DP), Allocatable :: w_pcp_red (:, :), w_temp_red (:, :)!reduced distance weigth matricies
 
 
  Real (DP), Allocatable :: Y_tmean (:), Y_trange (:)!transformed station data arrays
  Real (DP), Allocatable :: Y_tmean_red (:), Y_trange_red (:)!transformed station data arrays
  Real (DP), Allocatable :: stn_vals (:), stn_data (:, :), tair_data (:, :, :), stn_tair (:, :)!original station data arrays
  Real (DP), Allocatable :: auto_corr (:)!lag-1 autocorrelation for stations over entire time period used
  Real (DP), Allocatable :: t_p_corr (:)!correlation between temp and precip
  Integer (I4B), Allocatable :: Yp (:)!binary for logistic regression
  Integer (I4B), Allocatable :: Yp_red (:)!reduced binary for logistic regression
 
 
  Logical, Allocatable :: stn_miss (:), stn_miss_t (:)!missing value logical arrays
 
  Real (DP) :: m_tmean, m_trange !used to fill missing values in stations
 
  Real (DP) :: errsum, wgtsum, p, sta_pcp, sta_temp
  Real (DP) :: auto_corr_sum, tp_corr_sum
  Real (DP) :: step_mean, step_std, step_std_all, step_min, step_max !timestep statistics
 
  Real (DP) :: rsqr, ss_tot, ss_res, vif !r-squared and variance correction
 
  Integer (I4B) :: xsize !size of second dimension of input X array
 
  Integer (I4B) :: ntimes, nstns
  Integer (I4B) :: t, i, g, ndata, nodata
  Integer (I4B) :: ndata_t, nodata_t
  Integer (I4B) :: lag, window, tc, trc
  Integer (I4B) :: auto_cnt, tp_cnt
 
  Integer (I4B) :: stn_count
 
 
  !variables for tracking closest N stations for precipitation
  Integer (I4B) :: out_loc
  Integer (I4B), Parameter :: sta_limit = 30
  Integer (I4B), Allocatable :: close_loc (:, :)
  Integer (I4B), Allocatable :: close_count (:)
 
  Real (DP), Allocatable :: close_weights (:, :)
  Real (DP), Allocatable :: close_meta (:, :, :)
  Real (DP) :: min_weight
  Real (DP) :: max_distance
  Real (DP), Parameter :: search_distance = 1000.0
 
  !variables for tracking closest N stations for temperature
  Integer (I4B) :: out_loc_t
  Integer (I4B), Allocatable :: close_loc_t (:, :)
  Integer (I4B), Allocatable :: close_count_t (:)
 
  Real (DP), Allocatable :: close_weights_t (:, :)
  Real (DP), Allocatable :: close_meta_t (:, :, :)
  Real (DP) :: min_weight_t
  Real (DP) :: max_distance_t
 
 
  Real (DP) :: tmp_pcp
  Real (DP) :: tmp_weight
 
  Integer (I4B) :: slope_flag_pcp
  Integer (I4B) :: slope_flag_temp
 
!variables to check for singular matrix
  Real (DP), Allocatable :: tmp (:, :)
  Real (DP), Allocatable :: vv (:)
 
!variables for timing code
  Integer (I4B) :: t1, t2, count_rate
  Integer (I4B) :: tg1, tg2
 
 
!variables for keeping track of max_distance modifications
  Integer (I4B), Allocatable :: expand_flag (:), expand_flag_t (:)
  Real (DP), Allocatable :: expand_dist (:), expand_dist_t (:)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! code starts below here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  nstns = size (stnid)
!  nstns = 100
  ntimes = size (Times)
 
  xsize = size (X, 2)
 
!allocate variables
  Allocate (Y(nstns))
  Allocate (Y_tmean(nstns), Y_trange(nstns))
  Allocate (stn_data(nstns, ntimes))
  Allocate (tair_data(2, nstns, ntimes))
 
!original allocations
!  allocate(w_pcp(nstns,nstns))
!  allocate(w_temp(nstns,nstns))
!  allocate(TWX(4,nstns))
!  allocate(TX(4,nstns))
 
  Allocate (w_pcp_red(sta_limit, sta_limit))
  Allocate (w_temp_red(sta_limit, sta_limit))
  Allocate (Y_red(sta_limit))
  Allocate (Y_tmean_red(sta_limit), Y_trange_red(sta_limit))
  Allocate (X_red(sta_limit, xsize))
  Allocate (X_red_t(sta_limit, xsize))
 
  Allocate (w_pcp_1d(sta_limit))
  Allocate (w_temp_1d(sta_limit))
  Allocate (w_pcp_1d_loc(sta_limit))
  Allocate (w_temp_1d_loc(sta_limit))
 
  Allocate (tmp(6, 6))
  Allocate (vv(6))
 
  Allocate (PCP(ngrid, ntimes))
  Allocate (POP(ngrid, ntimes))
  Allocate (PCPERR(ngrid, ntimes))
 
  Allocate (tmean(ngrid, ntimes))
  Allocate (tmean_err(ngrid, ntimes))
  Allocate (trange(ngrid, ntimes))
  Allocate (trange_err(ngrid, ntimes))
 
 
  Allocate (PCP_2(ngrid, ntimes))
  Allocate (POP_2(ngrid, ntimes))
  Allocate (PCPERR_2(ngrid, ntimes))
 
  Allocate (tmean_2(ngrid, ntimes))
  Allocate (tmean_err_2(ngrid, ntimes))
  Allocate (trange_2(ngrid, ntimes))
  Allocate (trange_err_2(ngrid, ntimes))
 
 
  Allocate (auto_corr(nstns))
  Allocate (t_p_corr(nstns))
 
  Allocate (Yp(nstns))
  Allocate (Yp_red(sta_limit))
 
  !station limit arrays
  Allocate (close_weights(ngrid, sta_limit))
  Allocate (close_loc(ngrid, sta_limit))
  Allocate (close_meta(5, ngrid, sta_limit))
  Allocate (close_count(ngrid))
 
  Allocate (close_weights_t(ngrid, sta_limit))
  Allocate (close_loc_t(ngrid, sta_limit))
  Allocate (close_meta_t(5, ngrid, sta_limit))
  Allocate (close_count_t(ngrid))
 
 
  !base weight array
  Allocate (w_base(ngrid, nstns))
 
  !max_dist tracking variables
  Allocate (expand_dist(ngrid), expand_flag(ngrid))
  Allocate (expand_dist_t(ngrid), expand_flag_t(ngrid))
 
 
  PCP = 0.0d0
  POP = 0.0d0
  PCPERR = 0.0d0
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
 
  Do i = 1, nstns, 1
 
    Call read_station (stnvar, stnid(i), site_var, site_var_t, site_list, Times, stn_vals, stn_tair, stn_miss, &
   & stn_miss_t, error)
 
 
    stn_count = stn_count + 1
 
 
    stn_data (i, :) = stn_vals
    tair_data (1, i, :) = stn_tair (1, :)
    tair_data (2, i, :) = stn_tair (2, :)
 
    !call subroutine that does various  correlation calculations
    !can do autocorrelations and correlation between temperature and precipitation
    !uses an n-day moving average (window) to remove "monthly" cycle from temp
    !and computes autocorrelation on the anomalies
 
    lag = 1
    window = 31
    Call generic_corr (stn_data(i, :), tair_data(:, i, :), lag, window, auto_corr(i), t_p_corr(i))
 
    !compute mean autocorrelation for all stations and all times
    !check for values outside of -1 to 1
    !stations with incomplete data are set to -999
    If (auto_corr(i) .Ge.-1.0 .And. auto_corr(i) .Le. 1.0) Then
      auto_corr_sum = auto_corr_sum + auto_corr (i)
      auto_cnt = auto_cnt + 1
    End If
    If (t_p_corr(i) .Ge.-1.0 .And. t_p_corr(i) .Le. 1.0) Then
      tp_corr_sum = tp_corr_sum + t_p_corr (i)
      tp_cnt = tp_cnt + 1
    End If
 
    Deallocate (stn_miss_t)
    Deallocate (stn_miss)
    Deallocate (stn_vals)
    Deallocate (stn_tair)
  End Do !end station read loop
 
 
  error = 0
 
  mean_autocorr = auto_corr_sum / real (auto_cnt, kind(DP))
  mean_tp_corr = tp_corr_sum / real (tp_cnt, kind(DP))
 
  Print *, 'Temp lag-1 autocorrelation: ', mean_autocorr (1)
  Print *, 'Temp-precip correlation: ', mean_tp_corr (1)
 
 
  !pull weight generation outside of time loop
 
 
  Print *, 'Generating base weight matrix & '
  Print *, 'finding nearest stations for each gridpoint'
  Call system_clock (t1, count_rate)
 
 
  Do g = 1, ngrid, 1
 
 
    close_count (g) = 1
    min_weight = 0.0d0
    close_weights (g, :) = 0.0d0
 
    close_count_t (g) = 1
    min_weight_t = 0.0d0
    close_weights_t (g, :) = 0.0d0
 
    Do i = 1, nstns, 1
      !setup distinct weight matrices for precip and temperature
      Call calc_distance_weight (search_distance, X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), w_base(g, i))
 
!original call
!      call calc_distance_weight(maxDistance, X(i,2), X(i,3), Z(g,2), Z(g,3), w_temp(i,i))
 
  !also set some logic to limit the number of stations to the N closest
      min_weight = 0.0d0
 
      If (w_base(g, i) .Gt. min_weight .And. stn_data(i, 1) .Gt.-100.0d0) Then
        If (close_count(g) .Le. sta_limit) Then
 
          close_weights (g, close_count(g)) = w_base (g, i)
          close_loc (g, close_count(g)) = i
 
          close_meta (1, g, close_count(g)) = X (i, 2)
          close_meta (2, g, close_count(g)) = X (i, 3)
          close_meta (3, g, close_count(g)) = Z (g, 2)
          close_meta (4, g, close_count(g)) = Z (g, 3)
          Call calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta(5, g, close_count(g)))
 
          close_count (g) = close_count (g) + 1
        Else
          min_weight = minval (close_weights(g, :), 1)
          If (w_base(g, i) .Gt. min_weight) Then
            out_loc = minloc (close_weights(g, :), 1)
            close_weights (g, out_loc) = w_base (g, i)
            close_loc (g, out_loc) = i
 
            close_meta (1, g, out_loc) = X (i, 2)
            close_meta (2, g, out_loc) = X (i, 3)
            close_meta (3, g, out_loc) = Z (g, 2)
            close_meta (4, g, out_loc) = Z (g, 3)
            Call calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta(5, g, out_loc))
          End If
 
        End If
      End If
 
  !need to repeat above for temperature since that data is independent of precipitation
      min_weight_t = 0.0d0
 
      If (w_base(g, i) .Gt. min_weight_t .And. tair_data(1, i, 1) .Gt.-200.0d0) Then
        If (close_count_t(g) .Le. sta_limit) Then
 
          close_weights_t (g, close_count_t(g)) = w_base (g, i)
          close_loc_t (g, close_count_t(g)) = i
 
          close_meta_t (1, g, close_count_t(g)) = X (i, 2)
          close_meta_t (2, g, close_count_t(g)) = X (i, 3)
          close_meta_t (3, g, close_count_t(g)) = Z (g, 2)
          close_meta_t (4, g, close_count_t(g)) = Z (g, 3)
          Call calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta_t(5, g, close_count_t(g)))
 
          close_count_t (g) = close_count_t (g) + 1
        Else
          min_weight_t = minval (close_weights(g, :), 1)
          If (w_base(g, i) .Gt. min_weight) Then
            out_loc_t = minloc (close_weights_t(g, :), 1)
            close_weights_t (g, out_loc_t) = w_base (g, i)
 
            close_loc_t (g, out_loc_t) = i
 
            close_meta_t (1, g, out_loc_t) = X (i, 2)
            close_meta_t (2, g, out_loc_t) = X (i, 3)
            close_meta_t (3, g, out_loc_t) = Z (g, 2)
            close_meta_t (4, g, out_loc_t) = Z (g, 3)
            Call calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta(5, g, out_loc_t))
          End If
        End If
      End If
 
    End Do !end station loop
 
 
  End Do !end grid point loop
  Call system_clock (t2, count_rate)
  Print *, 'Elapsed time for weight generation: ', real (t2-t1) / real (count_rate)
 
  Do t = 1, ntimes, 1
 
    Call system_clock (tg1, count_rate)
    Print *, "TIME STEP = ", Times (t), " (", t, "/", ntimes, ")"
 
 
    Do i = 1, nstns, 1
 
      Y (i) = stn_data (i, t)
 
      Y_tmean (i) = (tair_data(1, i, t)+tair_data(2, i, t)) / 2.0d0
      Y_trange (i) = Abs (tair_data(2, i, t)-tair_data(1, i, t))
 
    End Do
 
!do power transformation on input vector
    Call normalize_Y (4.0d0, Y)
 
 
    Do g = 1, ngrid, 1
 
! call system_clock(tg1,count_rate)
 
      If (Z(g, 4) .Gt.-200) Then
 
        Allocate (TWX_red(6, sta_limit))
        Allocate (TX_red(6, sta_limit))
        Allocate (TWX_red_t(6, sta_limit))
        Allocate (TX_red_t(6, sta_limit))
 
!want to reset weights for closest sta_limit stations...
!recalc calc_distance_weight function for selected stations
!set max_distance equal to the farthest station distance
 
   !set data count integers and initialize reduced arrays to zero
! call system_clock(t1,count_rate)
        ndata = 0
        nodata = 0
        w_pcp_red = 0.0
        Y_red = 0.0
        X_red = 0.0
        Yp_red = 0
 
        max_distance = 0.0d0
        Do i = 1, (close_count(g)-1)
          If (close_meta(5, g, i) .Gt. max_distance) Then
            max_distance = close_meta (5, g, i)
          End If
        End Do
 
 
 
        If (max_distance .Le. maxDistance) Then
          expand_dist (g) = max_distance
          expand_flag (g) = 0
          max_distance = maxDistance
        Else
          max_distance = max_distance + 1.0d0
          expand_flag (g) = 1
          expand_dist (g) = max_distance
        End If
 
    !reduced matrices for precip
        slope_flag_pcp = 0
        Do i = 1, (close_count(g)-1)
          Call calc_distance_weight (max_distance, close_meta(1, g, i), close_meta(2, g, i), close_meta(3, g, i), &
         & close_meta(4, g, i), tmp_weight)
 
          w_pcp_red (i, i) = tmp_weight
          w_pcp_1d (i) = tmp_weight
          w_pcp_1d_loc (i) = close_loc (g, i)
          Y_red (i) = Y (close_loc(g, i))
          X_red (i, :) = X (close_loc(g, i), :)
 
 
 
          If (stn_data(close_loc(g, i), t) .Gt. 0.0) Then
            ndata = ndata + 1
            Yp_red (i) = 1
          Else
            nodata = nodata + 1
          End If
 
 
        End Do
 
 
        Call normalize_Xv (Y_red, w_pcp_1d, step_mean, step_std, step_std_all, step_min, step_max, Yp_red)
 
        Y_mean (g, t) = step_mean
        Y_std (g, t) = step_std
        y_std_all (g, t) = step_std_all
        y_min (g, t) = step_min
        y_max (g, t) = step_max
 
 
        ndata_t = 0
        nodata_t = 0
        w_temp_red = 0.0
        Y_tmean_red = 0.0
        Y_trange_red = 0.0
        X_red_t = 0.0
 
 
  ! max_distance_t = maxval(close_meta_t(5,g,:))
        max_distance_t = 0.0d0
        Do i = 1, (close_count_t(g)-1)
          If (close_meta_t(5, g, i) .Gt. max_distance_t) Then
            max_distance_t = close_meta_t (5, g, i)
          End If
        End Do
 
        If (max_distance_t .Le. maxDistance) Then
          max_distance_t = maxDistance
        Else
          max_distance_t = max_distance_t + 1.0d0
          expand_flag_t (g) = 1
          expand_dist_t (g) = max_distance_t
        End If
 
    !reduced matrices for temperature
        slope_flag_temp = 0
        Do i = 1, (close_count_t(g)-1)
          Call calc_distance_weight (max_distance_t, close_meta_t(1, g, i), close_meta_t(2, g, i), close_meta_t(3, g, &
         & i), close_meta_t(4, g, i), tmp_weight)
 
          w_temp_red (i, i) = tmp_weight
          w_temp_1d (i) = tmp_weight
          Y_tmean_red (i) = Y_tmean (close_loc_t(g, i))
          Y_trange_red (i) = Y_trange (close_loc_t(g, i))
          X_red_t (i, :) = X (close_loc_t(g, i), :)
 
          If (Y_tmean(close_loc_t(g, i)) .Gt.-100.0) Then
            ndata_t = ndata_t + 1
          Else
            nodata_t = nodata_t + 1
          End If
 
 
        End Do
 
        If (ndata == 0 .And. nodata == 0) Then
     !print *, "No stations within max distance of grid cell!"
          POP (g, t) = 0.0
          PCP (g, t) = 0.0
          PCPERR (g, t) = 0.0
 
          POP_2 (g, t) = 0.0
          PCP_2 (g, t) = 0.0
          PCPERR_2 (g, t) = 0.0
 
        End If
 
        If (ndata_t == 0 .And. nodata_t == 0) Then
          If (t .Gt. 1) Then
            tmean (g, t) = tmean (g, t-1)
            trange (g, t) = trange (g, t-1)
            tmean_err (g, t) = tmean_err (g, t-1)
            trange_err (g, t) = trange_err (g, t-1)
 
            tmean_2 (g, t) = tmean_2 (g, t-1)
            trange_2 (g, t) = trange_2 (g, t-1)
            tmean_err_2 (g, t) = tmean_err_2 (g, t-1)
            trange_err_2 (g, t) = trange_err_2 (g, t-1)
          Else
            tmean (g, t) = - 999
            trange (g, t) = - 999
            tmean_err (g, t) = 0.0
            trange_err (g, t) = 0.0
 
            tmean_2 (g, t) = - 999
            trange_2 (g, t) = - 999
            tmean_err_2 (g, t) = 0.0
            trange_err_2 (g, t) = 0.0
          End If
        End If
 
        If (ndata >= 1) Then
 
            !tmp needs to be matmul(TX,X) where TX = TWX_red and X = X_red
          TWX_red = matmul (transpose(X_red), w_pcp_red)
          tmp = matmul (TWX_red, X_red)
          vv = maxval (Abs(tmp), dim=2)
 
          If (any(vv == 0.0)) Then
            slope_flag_pcp = 0
          Else
            slope_flag_pcp = 1
          End If
 
          If (nodata == 0) Then
       !print *, "All stations have precip, POP = 1.0"
            POP (g, t) = 1.0
 
            POP_2 (g, t) = 1.0
 
          Else
 
 
 
            If (slope_flag_pcp .Eq. 0) Then
              POP (g, t) = - 999.
            Else
       !regression with slope
              TX_red = transpose (X_red)
              TWX_red = matmul (TX_red, w_pcp_red)
 
              Call logistic_regression (X_red, Y_red, TWX_red, Yp_red, B)
              POP (g, t) = real (1.0/(1.0+Exp(-dot_product(Z(g, :), B))), kind(SP))
 
              Deallocate (B)
            End If
 
 
       !regression without slope
            Deallocate (TWX_red)
            Deallocate (TX_red)
            Allocate (TWX_red(4, sta_limit))
            Allocate (TX_red(4, sta_limit))
            TX_red = transpose (X_red(:, 1:4))
 
            TWX_red = matmul (TX_red, w_pcp_red)
 
 
 
            Call logistic_regression (X_red(:, 1:4), Y_red, TWX_red, Yp_red, B)
            POP_2 (g, t) = real (1.0/(1.0+Exp(-dot_product(Z(g, 1:4), B))), kind(SP))
 
 
            Deallocate (B)
 
          End If
 
          Deallocate (TWX_red)
          Deallocate (TX_red)
          Allocate (TWX_red(6, sta_limit))
          Allocate (TX_red(6, sta_limit))
 
 
          If (slope_flag_pcp .Eq. 0) Then
            PCP (g, t) = - 999.
          Else
       !regression with slope
            TX_red = transpose (X_red)
            TWX_red = matmul (TX_red, w_pcp_red)
 
            Call least_squares (X_red, Y_red, TWX_red, B)
            PCP (g, t) = real (dot_product(Z(g, :), B), kind(SP))
 
 
            wgtsum = 0.0
            errsum = 0.0
            ss_tot = 0.0
            ss_res = 0.0
            Do i = 1, (close_count(g)-1), 1
              wgtsum = wgtsum + w_pcp_red (i, i)
 
 
              errsum = errsum + (w_pcp_red(i, i)*(PCP(g, t)-Y_red(i))**2)
 
            End Do
            PCPERR (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
          End If
 
 
      !regression without slope
 
          Deallocate (TWX_red)
          Deallocate (TX_red)
          Allocate (TWX_red(4, sta_limit))
          Allocate (TX_red(4, sta_limit))
 
          TX_red = transpose (X_red(:, 1:4))
          TWX_red = matmul (TX_red, w_pcp_red)
 
          Call least_squares (X_red(:, 1:4), Y_red, TWX_red, B)
          PCP_2 (g, t) = real (dot_product(Z(g, 1:4), B), kind(SP))
 
          wgtsum = 0.0
          errsum = 0.0
          ss_tot = 0.0
          ss_res = 0.0
          Do i = 1, (close_count(g)-1), 1
            wgtsum = wgtsum + w_pcp_red (i, i)
 
            errsum = errsum + (w_pcp_red(i, i)*(PCP_2(g, t)-Y_red(i))**2)
 
          End Do
          PCPERR_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
 
          Deallocate (B)
 
        Else
     !print *, "Not enough stations with data within max distance"
          POP (g, t) = 0.0
          PCP (g, t) = 0.0
          PCPERR (g, t) = 0.0
 
          POP_2 (g, t) = 0.0
          PCP_2 (g, t) = 0.0
          PCPERR_2 (g, t) = 0.0
        End If !precip if statement
 
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 !  do temperature ols now
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
        If (ndata_t .Ge. 1) Then
 
 
   !regression with slope
          TX_red_t = transpose (X_red_t)
          TWX_red_t = matmul (TX_red_t, w_temp_red)
          Call least_squares (X_red_t, Y_tmean_red, TWX_red_t, B)
          tmean (g, t) = real (dot_product(Z(g, :), B), kind(SP))
 
          errsum = 0.0
          wgtsum = 0.0
          Do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
 
            errsum = errsum + (w_temp_red(i, i)*(tmean(g, t)-Y_tmean_red(i))**2)
          End Do
          tmean_err (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
 
   !regression without slope
          Deallocate (B)
          Deallocate (TWX_red_t)
          Deallocate (TX_red_t)
          Allocate (TWX_red_t(4, sta_limit))
          Allocate (TX_red_t(4, sta_limit))
          TX_red_t = transpose (X_red_t(:, 1:4))
          TWX_red_t = matmul (TX_red_t, w_temp_red)
 
          Call least_squares (X_red_t(:, 1:4), Y_tmean_red, TWX_red_t, B)
          tmean_2 (g, t) = real (dot_product(Z(g, 1:4), B), kind(SP))
 
          errsum = 0.0
          wgtsum = 0.0
   !do i = 1, nstns, 1
          Do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
 
 
            errsum = errsum + (w_temp_red(i, i)*(tmean_2(g, t)-Y_tmean_red(i))**2)
 
          End Do
          tmean_err_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
          Deallocate (B)
 
   !then do trange
 
          Deallocate (TWX_red_t)
          Deallocate (TX_red_t)
          Allocate (TWX_red_t(6, sta_limit))
          Allocate (TX_red_t(6, sta_limit))
          TX_red_t = transpose (X_red_t)
          TWX_red_t = matmul (TX_red_t, w_temp_red)
 
          Call least_squares (X_red_t, Y_trange_red, TWX_red_t, B)
          trange (g, t) = real (dot_product(Z(g, :), B), kind(SP))
 
          errsum = 0.0
          wgtsum = 0.0
          Do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
 
            errsum = errsum + (w_temp_red(i, i)*(trange(g, t)-Y_trange_red(i))**2)
 
          End Do
          trange_err (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
          !regression without slope
          Deallocate (B)
          Deallocate (TWX_red_t)
          Deallocate (TX_red_t)
          Allocate (TWX_red_t(4, sta_limit))
          Allocate (TX_red_t(4, sta_limit))
          TX_red_t = transpose (X_red_t(:, 1:4))
          TWX_red_t = matmul (TX_red_t, w_temp_red)
 
          Call least_squares (X_red_t(:, 1:4), Y_trange_red, TWX_red_t, B)
          trange_2 (g, t) = real (dot_product(Z(g, 1:4), B), kind(SP))
 
          errsum = 0.0
          wgtsum = 0.0
          Do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
 
            sta_temp = real (dot_product(X_red_t(i, 1:4), B), kind(SP))
 
            errsum = errsum + (w_temp_red(i, i)*(trange_2(g, t)-Y_trange_red(i))**2)
 
          End Do
          trange_err_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
        Else
 
     !if not enough stations with data
     !just use value from previous grid point for now
          If (g .Gt. 1) Then
            trange (g, t) = trange (g-1, t)
            trange_err (g, t) = trange_err (g-1, t)
            tmean (g, t) = tmean (g-1, t)
            tmean_err (g, t) = tmean_err (g-1, t)
 
            trange_2 (g, t) = trange_2 (g-1, t)
            trange_err_2 (g, t) = trange_err_2 (g-1, t)
            tmean_2 (g, t) = tmean_2 (g-1, t)
            tmean_err_2 (g, t) = tmean_err_2 (g-1, t)
          Else
            trange (g, t) = trange (g, t-1)
            trange_err (g, t) = trange_err (g-1, t-1)
            tmean (g, t) = tmean (g, t-1)
            tmean_err (g, t) = tmean_err (g, t-1)
 
            trange_2 (g, t) = trange_2 (g, t-1)
            trange_err_2 (g, t) = trange_err_2 (g-1, t-1)
            tmean_2 (g, t) = tmean_2 (g, t-1)
            tmean_err_2 (g, t) = tmean_err_2 (g, t-1)
          End If
        End If !end data check if statement for temperature
 
 
      End If !end check for valid elevation
 
      If (allocated(TWX_red)) Then
        Deallocate (TWX_red)
      End If
      If (allocated(TX_red)) Then
        Deallocate (TX_red)
      End If
      If (allocated(TWX_red_t)) Then
        Deallocate (TWX_red_t)
      End If
      If (allocated(TX_red_t)) Then
        Deallocate (TX_red_t)
      End If
 
    End Do !end grid loop
 
    Call system_clock (tg2, count_rate)
    Print *, 'Elapsed time for one time step: ', real (tg2-tg1) / real (count_rate)
  End Do !end time loop
 
End Subroutine estimate_precip
 
 
Subroutine normalize_X (X)
  Use type
  Implicit None
 
  Real (DP), Intent (Inout) :: X (:, :)
 
  Real (DP) :: mean, stdev, sum_x, sum_x2
  Integer (I4B) :: v, t
  Integer (I4B) :: nvars, ntimes
 
  nvars = size (X, 2) - 1
  ntimes = size (X, 1)
 
  Do v = 2, nvars + 1, 1
    sum_x = 0.0d0
    sum_x2 = 0.0d0
    Do t = 1, ntimes, 1
      sum_x = sum_x + X (t, v)
      sum_x2 = sum_x2 + X (t, v) ** 2
    End Do
    mean = sum_x / real (ntimes)
    stdev = Sqrt ((real(ntimes)*sum_x2-sum_x**2)/(real(ntimes)*real(ntimes-1)))
    Do t = 1, ntimes, 1
      If (stdev .Eq. 0.0) Then
        X (t, v) = X (t, v)
      Else
        X (t, v) = (X(t, v)-mean) / stdev
      End If
    End Do
  End Do
 
End Subroutine normalize_X
 
 
Subroutine normalize_Xv (X, weight, mean, stdev, stdev_all, smin, smax, Yp)
  Use type
  Implicit None
 
  Real (DP), Intent (Inout) :: X (:)
  Real (DP), Intent (In) :: weight (:)
  Real (DP), Intent (Out) :: mean
  Real (DP), Intent (Out) :: stdev
  Real (DP), Intent (Out) :: stdev_all
  Real (DP), Intent (Out) :: smin
  Real (DP), Intent (Out) :: smax
  Integer (I4B), Intent (Out) :: Yp (:)
 
!  real(DP) :: mean, stdev, sum_x2, sum_x
  Real (DP) :: sum_x2, sum_x
  Real (DP) :: sum_stdev, sum_std
 
  Real (DP) :: mean_all, sum_weight, sum_xw
 
  Integer (I4B) :: t
  Integer (I4B) :: ntimes
 
  ntimes = size (X)
 
  Yp = 0
  smin = 9999.0
  smax = 0.0
  stdev = 0.01
  mean = 0.0
  sum_std = 0.0
  sum_stdev = 0.0
 
  sum_xw = 0.0d0
  sum_weight = 0.0d0
 
  sum_x = 0.0d0
  sum_x2 = 0.0d0
  Do t = 1, ntimes, 1
    If (X(t) > 0.0) Then
      sum_x = sum_x + X (t)
      sum_x2 = sum_x2 + X (t) ** 2
      sum_xw = sum_xw + weight (t) * X (t)
      sum_weight = sum_weight + weight (t)
      Yp (t) = 1
 
      If (X(t) .Le. smin) Then
        smin = X (t)
      End If
      If (X(t) .Ge. smax) Then
        smax = X (t)
      End If
 
    End If
  End Do
 
  mean_all = sum_x / real (ntimes)
  mean = sum_x / real (sum(Yp))
  Do t = 1, ntimes, 1
    sum_stdev = sum_stdev + (X(t)-mean_all) ** 2
    If (Yp(t) .Eq. 1) Then
      sum_std = sum_std + (X(t)-mean) ** 2
    End If
  End Do
!  stdev_all = sqrt(sum_stdev/(ntimes-1))
 
 
  If (sum(Yp) .Ge. 2) Then
    stdev = Sqrt (sum_std/real(sum(Yp)-1.0))
    stdev_all = Sqrt (sum_std/real(ntimes-1.0))
 
  Else
    mean = sum_x
    stdev = 0.01
!    stdev_all = sum_xw
    stdev_all = 0.01
  End If
 
 
  If (stdev .Gt. 0.0) Then
    Do t = 1, ntimes, 1
      If (Yp(t) .Gt. 0.0) Then
        X (t) = (X(t)-mean) / stdev
      End If
    End Do
  End If
 
 
  Return
 
End Subroutine normalize_Xv
 
 
!
! Normalize a vector.
! Input :
!  Y   = An n-element vector.
!  exp = Integer exponent type (0 to 2).
!          Exponent types available:
!              0 = Exponent of 1, do not normalize.
!              1 = Exponent of 1/2.
!              2 = Exponent of 1/3rd.
! Output:
!  Y   = Input vector is also output vector.
Subroutine normalize_Y (texp, Y)
  Use type
  Implicit None
 
  Real (DP), Intent (In) :: texp !transform exponent
  Real (DP), Intent (Inout) :: Y (:)
  Integer (I4B) :: t
  Integer (I4B) :: ntimes
 
  ntimes = size (Y)
 
  Do t = 1, ntimes, 1
!     Y(t) = Y(t) ** (1.0d0/2.5d0)
!    Y(t) = Y(t) ** (1.0d0/4d0)
    Y (t) = Y (t) ** (1.0d0/texp)
  End Do
 
End Subroutine normalize_Y
 
!
! Calculate the weights for time step tt.
! Dates within exclusion, same day/week/month as time tt, will have zero weight.
!
! Input :
!   Times = A t-element vector containing the unix time of each time step.
!   tt    = An integer containing the index of the time to calculate weights at.
!   X     = A t by n array containing the input values.
!   excl  = Integer Exclusion method (0 to 3).
!             Exclusion methods available:
!                 0 = no exclusion, set all weights to 1
!                 1 = exclude day, 24Hr period
!                 2 = exclude week, 7 Day period
!                 3 = exclude month, 30 Day period
! Output:
!   W     = A t by t array containing the diagonal weights for each time step.
Subroutine calc_weights (Times, tt, X, W)
  Use type
  Implicit None
 
  Real (DP), Intent (In) :: Times (:)
  Integer (I4B), Intent (In) :: tt
  Real (DP), Intent (In) :: X (:, :)
  Real (DP), Allocatable, Intent (Out) :: W (:, :)
  Real (DP) :: sum
  Integer (I4B) :: v, t
  Integer (I4B) :: nvars, ntimes
 
  nvars = size (X, 2) - 1
  ntimes = size (X, 1)
 
  Allocate (W(ntimes, ntimes))
 
  Do t = 1, ntimes, 1
    W (t, :) = 0.0d0
    If (t /= tt) Then
      sum = 0.0d0
      Do v = 2, nvars + 1, 1
        sum = sum + (X(t, v)-X(tt, v)) ** 2
      End Do
      W (t, :) = 0.0d0
 
      W (t, t) = 1 / Sqrt (sum/nvars)
    End If
  End Do
 
End Subroutine calc_weights
 
! Great circle distance calculation
! Output in nm
Subroutine calc_distance_weight (maxd, lat1, lon1, lat2, lon2, weight)
  Use type
  Implicit None
 
  Real (DP), Intent (In) :: maxd, lat1, lon1, lat2, lon2
  Real (DP), Intent (Out) :: weight
 
  Real (DP) :: dist, lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927
 
  lat1r = lat1 * Pi / 180.0d0
  lon1r = lon1 * Pi / 180.0d0
  lat2r = lat2 * Pi / 180.0d0
  lon2r = lon2 * Pi / 180.0d0
  dist = ((180*60)/Pi) * (2*Asin(Sqrt((Sin((lat1r-lat2r)/2))**2+Cos(lat1r)*Cos(lat2r)*(Sin((lon1r-lon2r)/2))**2)))
  If (dist .Gt. maxd) Then
    weight = 0.0d0
  Else
    weight = (1.0d0-(dist/maxd)**3) ** 3
!    weight = 1.0d0 - (dist/maxd)**0.5
  End If
 
End Subroutine calc_distance_weight
 
! Great circle distance calculation
! Output in nm
Subroutine calc_distance (lat1, lon1, lat2, lon2, dist)
  Use type
  Implicit None
 
  Real (DP), Intent (In) :: lat1, lon1, lat2, lon2
  Real (DP), Intent (Out) :: dist
 
  Real (DP) :: lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927
 
  lat1r = lat1 * Pi / 180.0d0
  lon1r = lon1 * Pi / 180.0d0
  lat2r = lat2 * Pi / 180.0d0
  lon2r = lon2 * Pi / 180.0d0
  dist = ((180*60)/Pi) * (2*Asin(Sqrt((Sin((lat1r-lat2r)/2))**2+Cos(lat1r)*Cos(lat2r)*(Sin((lon1r-lon2r)/2))**2)))
 
End Subroutine calc_distance
 
!
! Solve linear equation for x (Ax = b => x = bA^-1) using LU decomposition and back substitution.
! Input:
!   X  = An m by n array.
!   TX = Precalculated transpose array of X, size n by m
!   Y  = An m-element vector containing the right-hand side of the linear system Ax = b.
! Output:
!   B  = An n-element vector.
Subroutine least_squares (X, Y, TX, B)
  Use type
  Implicit None
 
  Interface
    Subroutine ludcmp (a, indx, D)
      Use type
      Real (DP), Dimension (:, :), Intent (Inout) :: a
      Integer (I4B), Dimension (:), Intent (Out) :: indx
      Real (DP), Intent (Out) :: D
    End Subroutine ludcmp
 
    Subroutine lubksb (a, indx, B)
      Use type
      Real (DP), Dimension (:, :), Intent (In) :: a
      Integer (I4B), Dimension (:), Intent (In) :: indx
      Real (DP), Dimension (:), Intent (Inout) :: B
    End Subroutine lubksb
  End Interface
 
  Real (DP), Intent (In) :: X (:, :)
  Real (DP), Intent (In) :: Y (:)
  Real (DP), Intent (In) :: TX (:, :)
  Real (DP), Allocatable, Intent (Out) :: B (:)
 
  Real (DP), Allocatable :: a (:, :)
  Integer (I4B), Allocatable :: indx (:)
  Integer (I4B) :: nvars, ntimes
  Real (DP) :: D
 
  nvars = size (X, 2) - 1
  ntimes = size (Y)
 
  Allocate (B(nvars+1))
  Allocate (a(nvars+1, nvars+1))
  Allocate (indx(nvars+1))
 
  B = matmul (TX, Y)
  a = matmul (TX, X)
 
  Call ludcmp (a, indx, D)
  If (any(Abs(a) < 9.99999968E-15)) Then
    B (:) = 0.0d0
    Print *, "Warning, LUdcmp produced a zero."
    Return
  End If
 
  Call lubksb (a, indx, B)
 
  Deallocate (a)
  Deallocate (indx)
 
End Subroutine least_squares
 
 
Subroutine logistic_regression (X, Y, TX, Yp, B)
  Use type
  Implicit None
 
  Interface
    Subroutine least_squares (X, Y, TX, B)
      Use type
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Intent (In) :: Y (:)
      Real (DP), Intent (In) :: TX (:, :)
      Real (DP), Allocatable, Intent (Out) :: B (:)
    End Subroutine least_squares
  End Interface
 
  Real (DP), Intent (In) :: X (:, :)
  Real (DP), Intent (In) :: Y (:)
  Real (DP), Intent (In) :: TX (:, :)
  Integer (I4B), Intent (In) :: Yp (:)
  Real (DP), Allocatable, Intent (Out) :: B (:)
 
  Real (DP), Allocatable :: Ypd (:), p (:), YN (:), BN (:), v (:, :), XV (:, :)
!  real(DP), allocatable :: P(:), YN(:), BN(:), V(:,:), XV(:,:)
  Integer (I4B) :: nvars, ntimes, i, t, f, it
  Real (DP) :: D
 
  nvars = size (X, 2) - 1
  ntimes = size (Y)
 
  Allocate (B(nvars+1))
  Allocate (Ypd(ntimes))
  Allocate (YN(ntimes))
  Allocate (p(ntimes))
  Allocate (v(ntimes, ntimes))
  Allocate (XV(ntimes, nvars+1))
 
  Do t = 1, ntimes, 1
    If (Yp(t) .Gt. 0.0) Then
      Ypd (t) = 1.0d0
    Else
      Ypd (t) = 0.0d0
    End If
  End Do
 
  B = 0.0d0
  i = 0
  it = 0
 
  Do while (f /=  1)
!     print *, "Iteration ", it
    p = 1.0d0 / (1.0d0+Exp(-matmul(X, B)))
    If (any(p > 0.97)) Then
!    PRINT *, "WARNING: logistic regression diverging"
      f = 1
    Else
 
      YN = Ypd - p
      v = 0.0d0
      Do t = 1, ntimes, 1
        v (t, t) = p (t) * (1.0d0-p(t))
      End Do
      XV = matmul (v, X)
      Call least_squares (XV, YN, TX, BN)
 
      f = 1
      Do i = 1, nvars + 1, 1
        If (BN(i) .Gt. 1.0E-04 .Or. BN(i) .Lt.-1.0E-04) Then
          f = 0
        End If
      End Do
      If (it > 8) Then
!     PRINT *, "WARNING: logistic regression failed to converge"
        f = 1
      End If
 
      B = B + BN
!        print *, "Bnew = ", B
      Deallocate (BN)
 
    End If
    it = it + 1
  End Do
 
End Subroutine logistic_regression
 
 
Subroutine logistic_regressionrf (X, Y, TX, B)
  Use type
  Implicit None
 
  Interface
    Subroutine least_squares (X, Y, TX, B)
      Use type
      Real (DP), Intent (In) :: X (:, :)
      Real (DP), Intent (In) :: Y (:)
      Real (DP), Intent (In) :: TX (:, :)
      Real (DP), Allocatable, Intent (Out) :: B (:)
    End Subroutine least_squares
  End Interface
 
  Real (DP), Intent (In) :: X (:, :)
  Real (DP), Intent (In) :: Y (:)
  Real (DP), Intent (In) :: TX (:, :)
  Real (DP), Allocatable, Intent (Out) :: B (:)
 
  Real (DP), Allocatable :: Ypd (:), p (:), YN (:), BN (:), v (:, :), XV (:, :)
  Integer (I4B) :: nvars, ntimes, i, t, f, it
  Real (DP) :: D
 
  nvars = size (X, 2) - 1
  ntimes = size (Y)
 
  Allocate (B(nvars+1))
  Allocate (Ypd(ntimes))
  Allocate (YN(ntimes))
  Allocate (p(ntimes))
  Allocate (v(ntimes, ntimes))
  Allocate (XV(ntimes, nvars+1))
 
  Do t = 1, ntimes, 1
    If (Y(t) .Ne. 0.0) Then
      Ypd (t) = 1.0d0
    Else
      Ypd (t) = 0.0d0
    End If
  End Do
 
  B = 0.0d0
  i = 0
  it = 0
  !print *, "B = ", B
  Do while (f /=  1)
!     print *, "Iteration ", it
    p = 1.0d0 / (1.0d0+Exp(-matmul(X, B)))
!     print *, "Pie = ", P
    If (any(p > 0.97)) Then
!    PRINT *, "WARNING: logistic regression diverging"
      f = 1
    Else
 
      YN = Ypd - p
      v = 0.0d0
      Do t = 1, ntimes, 1
        v (t, t) = p (t) * (1.0d0-p(t))
      End Do
      XV = matmul (v, X)
 
      Call least_squares (XV, YN, TX, BN)
 
      f = 1
      Do i = 1, nvars + 1, 1
        If (BN(i) .Gt. 1.0E-04 .Or. BN(i) .Lt.-1.0E-04) Then
          f = 0
        End If
      End Do
      If (it > 8) Then
!     PRINT *, "WARNING: logistic regression failed to converge"
        f = 1
      End If
 
      B = B + BN
      Deallocate (BN)
 
    End If
    it = it + 1
  End Do
 
 
End Subroutine logistic_regressionrf
 
Subroutine generic_corr (stn_data, tair_data, lag, window, auto_corr, t_p_corr)
  Use type
 
  Implicit None
 
!input
  Real (DP), Intent (In) :: stn_data (:)
  Real (DP), Intent (In) :: tair_data (:, :)
  Integer (I4B), Intent (In) :: lag
  Integer (I4B), Intent (In) :: window
 
 
!output
  Real (DP), Intent (Out) :: auto_corr
  Real (DP), Intent (Out) :: t_p_corr
 
!local variables
  Real (DP), Allocatable :: tmean (:), trange (:)
  Real (DP), Allocatable :: moving_avg (:, :)
  Real (DP) :: lag_0_mean
  Real (DP) :: lag_n_mean
  Real (DP) :: lag_0_var
  Real (DP) :: lag_n_var
  Real (DP) :: lag_0_sum
  Real (DP) :: lag_n_sum
  Real (DP) :: cov
  Real (DP) :: lag_0_pmean, lag_0_pvar, lag_0_psum
  Real (DP) :: trange_mean, trange_sum, trange_var
 
  Real (DP) :: tmp_tmean, tmp_trange
 
  Integer (I4B) :: i, j, tmp_cnt
  Integer (I4B) :: ntimes
  Integer (I4B) :: half_window
  Integer (I4B) :: cnt_sums
  Integer (I4B) :: data_cnt
 
!code
  ntimes = size (stn_data)
 
  Allocate (tmean(ntimes))
  Allocate (trange(ntimes))
  Allocate (moving_avg(2, ntimes))
 
  data_cnt = 0
  Do i = 1, ntimes, 1
    If (tair_data(1, i) .Gt.-100.0 .And. tair_data(2, i) .Gt.-100.0) Then
      tmean (i) = ((tair_data(1, i)+tair_data(2, i))/2.d0) + 273.15d0
      trange (i) = (tair_data(2, i)-tair_data(1, i)/2.d0)
      data_cnt = data_cnt + 1
    Else
      tmean (i) = - 999.0d0
      trange (i) = - 999.0d0
    End If
  End Do
 
 
  half_window = floor (window/2.0d0)
 
!do the lag correlation for temperature
!first compute the moving average for climo removal
!need to check for missing values....
 
  Do i = 1, ntimes, 1
    If (i .Lt. half_window) Then
      tmp_tmean = 0.0
      tmp_trange = 0.0
      tmp_cnt = 0
      Do j = 1, window, 1
        If (tmean(j) .Gt.-100.0) Then
          tmp_tmean = tmp_tmean + tmean (j)
          tmp_trange = tmp_trange + trange (j)
          tmp_cnt = tmp_cnt + 1
        End If
        If (tmp_cnt .Gt. 0) Then
          moving_avg (1, i) = tmp_tmean / real (tmp_cnt, kind(DP))
          moving_avg (2, i) = tmp_trange / real (tmp_cnt, kind(DP))
        Else
          moving_avg (1, i) = - 999.0
          moving_avg (2, i) = - 999.0
        End If
      End Do
 
    Else If (i .Gt. ntimes-half_window) Then
      tmp_tmean = 0.0
      tmp_trange = 0.0
      tmp_cnt = 0
      Do j = 1, window, 1
        If (tmean(ntimes-j) .Gt.-100.0) Then
          tmp_tmean = tmp_tmean + tmean (ntimes-j)
          tmp_trange = tmp_trange + trange (ntimes-j)
          tmp_cnt = tmp_cnt + 1
        End If
      End Do
      If (tmp_cnt .Gt. 0) Then
        moving_avg (1, i) = tmp_tmean / real (tmp_cnt, kind(DP))
        moving_avg (2, i) = tmp_trange / real (tmp_cnt, kind(DP))
      Else
        moving_avg (1, i) = - 999.0
        moving_avg (2, i) = - 999.0
      End If
    Else
      tmp_tmean = 0.0
      tmp_trange = 0.0
      tmp_cnt = 0
      Do j = - half_window, half_window, 1
        If (tmean(i+j) .Gt.-100.0) Then
          tmp_tmean = tmp_tmean + tmean (i+j)
          tmp_trange = tmp_trange + trange (i+j)
          tmp_cnt = tmp_cnt + 1
        End If
      End Do
      If (tmp_cnt .Gt. 0) Then
        moving_avg (1, i) = tmp_tmean / real (tmp_cnt, kind(DP))
        moving_avg (2, i) = tmp_trange / real (tmp_cnt, kind(DP))
      Else
        moving_avg (1, i) = - 999.0
        moving_avg (2, i) = - 999.0
      End If
    End If
  End Do
 
 
!only use portions of timeseries to compute auto_corr
!need to go through and check to see if values exist for lag-0 and lag-n and moving_avg
!if values do not exist for any of the three, don't add to running sums
 
!compute means
  lag_0_sum = 0.0d0
  lag_n_sum = 0.0d0
  lag_0_var = 0.0d0
  lag_n_var = 0.0d0
  cov = 0.0d0
  cnt_sums = 0
 
  Do i = lag + 1, ntimes, 1
    If (tmean(i) .Gt.-100.0 .And. tmean(i-lag) .Gt.-100.0 .And. moving_avg(1, i) .Gt.-100.0 .And. moving_avg(1, i-lag) &
   & .Gt.-100.0) Then
      lag_n_sum = lag_n_sum + tmean (i-lag) - moving_avg (1, i-lag)
      lag_0_sum = lag_0_sum + tmean (i) - moving_avg (1, i)
      cnt_sums = cnt_sums + 1
    End If
  End Do
 
  lag_0_mean = lag_0_sum / real (cnt_sums, kind(DP))
  lag_n_mean = lag_n_sum / real (cnt_sums, kind(DP))
 
!compute variance,covariance
  Do i = lag + 1, ntimes, 1
    If (tmean(i) .Gt.-100.0 .And. tmean(i-lag) .Gt.-100.0 .And. moving_avg(1, i) .Gt.-100.0 .And. moving_avg(1, i-lag) &
   & .Gt.-100.0) Then
      lag_n_var = lag_n_var + ((tmean(i-lag)-moving_avg(1, i-lag))-lag_n_mean) ** 2
      lag_0_var = lag_0_var + ((tmean(i)-moving_avg(1, i))-lag_0_mean) ** 2
      cov = cov + ((tmean(i-lag)-moving_avg(1, i-lag))-lag_n_mean) * ((tmean(i)-moving_avg(1, i))-lag_0_mean)
    End If
  End Do
 
!compute autocorrelation
  auto_corr = cov / (Sqrt(lag_0_var)*Sqrt(lag_n_var))
 
 
!!!!!!!!!!!!!!!!
!
! now do the t - p correlation
! do correlation on trange, not tmean
!
!!!!!!!!!!!!!!!!!!!
 
  lag_0_pmean = 0.0d0
  lag_0_pvar = 0.0d0
  lag_0_psum = 0.0d0
  trange_sum = 0.0d0
  trange_mean = 0.0d0
  trange_var = 0.0d0
  cov = 0.0d0
  tmp_cnt = 0
 
 
!again need to check for missing values....
  Do i = 1, ntimes, 1
    If (trange(i) .Gt.-100 .And. stn_data(i) .Gt.-100.0) Then
      !compute for precip mean
      lag_0_psum = lag_0_psum + stn_data (i)
      !anomaly means of trange
      trange_sum = trange_sum + (trange(i)-moving_avg(2, i))
      tmp_cnt = tmp_cnt + 1
    End If
  End Do
  lag_0_pmean = lag_0_psum / real (tmp_cnt, kind(DP))
  trange_mean = trange_sum / real (tmp_cnt, kind(DP))
 
!compute variance and covariance
  Do i = 1, ntimes, 1
    If (trange(i) .Gt.-100 .And. stn_data(i) .Gt.-100.0) Then
      lag_0_pvar = lag_0_pvar + (stn_data(i)-lag_0_pmean) ** 2
      trange_var = trange_var + ((trange(i)-moving_avg(2, i))-trange_mean) ** 2
      cov = cov + ((trange(i)-moving_avg(2, i))-trange_mean) * (stn_data(i)-lag_0_pmean)
    End If
  End Do
 
  t_p_corr = cov / (Sqrt(lag_0_pvar)*Sqrt(trange_var))
  If (Sqrt(lag_0_pvar)*Sqrt(trange_var) .Le. 0.00001) Then
    t_p_corr = 0.0
  End If
!in some situations, there are very limited data used for calculations
!set those cases to missing value
  If (data_cnt .Lt. real(ntimes)*0.25) Then
    auto_corr = - 999.0
    t_p_corr = - 999.0
  End If
 
End Subroutine generic_corr
