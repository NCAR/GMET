SUBROUTINE estimate_coefficients (D, nvars, Lats, Lons, Times, stnid, stnlat, stnlon, stnalt, stnvar, site_var, &
& site_list, C, POC, error)
  USE type
  IMPLICIT NONE
 
  INTERFACE
 
   SUBROUTINE read_station (stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss, &
  & vals_miss_t, error)
    USE type
    CHARACTER (LEN=100), INTENT (IN) :: stnvar
    CHARACTER (LEN=100), INTENT (IN) :: stnid
    CHARACTER (LEN=100), INTENT (IN) :: site_var, site_var_t
    CHARACTER (LEN=500), INTENT (IN) :: site_list
    REAL (DP), INTENT (IN) :: Times (:)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: vals (:), tair_vals (:, :)
    LOGICAL, ALLOCATABLE, INTENT (OUT) :: vals_miss (:), vals_miss_t (:)
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_station
 
   SUBROUTINE normalize_X (X)
    USE type
    REAL (DP), INTENT (INOUT) :: X (:, :)
   END SUBROUTINE normalize_X
 
   SUBROUTINE normalize_Y (texp, Y)
    USE type
    REAL (DP), INTENT (IN) :: texp !transform exponent
    REAL (DP), INTENT (INOUT) :: Y (:)
   END SUBROUTINE normalize_Y
 
   SUBROUTINE calc_weights (Times, tt, X, W)
    USE type
    REAL (DP), INTENT (IN) :: Times (:)
    INTEGER (I4B), INTENT (IN) :: tt
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: W (:, :)
   END SUBROUTINE calc_weights
 
   SUBROUTINE least_squares (X, Y, TX, B)
    USE type
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), INTENT (IN) :: Y (:)
    REAL (DP), INTENT (IN) :: TX (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
   END SUBROUTINE least_squares
 
   SUBROUTINE logistic_regressionrf (X, Y, TX, B)
    USE type
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), INTENT (IN) :: Y (:)
    REAL (DP), INTENT (IN) :: TX (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
   END SUBROUTINE logistic_regressionrf
 
  END INTERFACE
 
  REAL (DP), INTENT (IN) :: D (:, :, :), Lats (:), Lons (:)
  REAL (DP), INTENT (IN) :: Times (:)
  INTEGER (I4B), INTENT (IN) :: nvars
  CHARACTER (LEN=100), INTENT (IN) :: stnid (:)
  REAL (DP), INTENT (IN) :: stnlat (:), stnlon (:), stnalt (:)
  CHARACTER (LEN=100), INTENT (IN) :: stnvar
  CHARACTER (LEN=100), INTENT (IN) :: site_var
  CHARACTER (LEN=500), INTENT (IN) :: site_list
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: C (:, :, :), POC (:, :, :)
  INTEGER, INTENT (OUT) :: error
 
  REAL (DP), ALLOCATABLE :: X (:, :), Y (:), XS (:, :), XP (:, :), TWXS (:, :), YS (:), W (:, :), B (:), tair_vals (:, &
 & :)
 
  REAL (DP), ALLOCATABLE :: TS (:)
  LOGICAL, ALLOCATABLE :: Y_miss (:), Y_miss_t (:)
  INTEGER (I4B) :: ntimes, nlats, nlons, nstns
  INTEGER (I4B) :: i, j, k, t, v, tt
  REAL (DP) :: minlondis, minlatdis
  INTEGER (I4B) :: minlatk, minlonj, gridindex
  INTEGER (I4B) :: vshape (2)
  CHARACTER (LEN=100) :: site_var_t
  REAL (DP) :: transform_exp !precipitation transform variable, the exponent of the transform  norm_pcp = pcp^(1/transform_exp)
 
 
  error = 0
  transform_exp = 4.0d0
 
  ntimes = size (Times)
  nlats = size (Lats)
  nlons = size (Lons)
  nstns = size (stnlat)
 
  ALLOCATE (X(ntimes, nvars+1))
  IF (trim(stnvar) .EQ. "PRCP") THEN
   ALLOCATE (POC(nstns, ntimes, nvars+1))
   POC = 0.0d0
  END IF
  ALLOCATE (C(nstns, ntimes, nvars+1))
  C = 0.0d0
 
  DO i = 1, nstns, 1
 
   minlondis = 360.0
   minlonj = - 1
   DO j = 1, nlons, 1
    IF (Abs(stnlon(i)-Lons(j)) < minlondis) THEN
     minlondis = Abs (stnlon(i)-Lons(j))
     minlonj = j
    END IF
   END DO
   minlatdis = 180.0
   minlatk = - 1
   DO k = 1, nlats, 1
    IF (Abs(stnlat(i)-Lats(k)) < minlatdis) THEN
     minlatdis = Abs (stnlat(i)-Lats(k))
     minlatk = k
    END IF
   END DO
 
   IF (minlonj ==-1 .OR. minlatk ==-1) THEN
    PRINT *, "Failed to find closest grid point for station: ", trim (stnid(i))
    error = 1
    RETURN
   END IF
 
   gridindex = ((minlatk-1)*nlons) + minlonj
 
   PRINT *, "Station: ", trim (stnid(i)), stnlat (i), stnlon (i)
   PRINT *, "Closest Grid point: ", gridindex, Lats (minlatk), Lons (minlonj)
 
   X (:, 1) = 1.0
   DO v = 1, nvars, 1
    X (:, v+1) = D (v, gridindex, :)
    DO t = 1, ntimes, 1
     IF (X(t, v+1) /= 0.0) EXIT
     IF (t == ntimes) THEN
      PRINT *, "ERROR: var ", v, " is all zero for station ", i
      X (:, v+1) = 1.0
     END IF
    END DO
   END DO
 
   CALL read_station (stnvar, stnid(i), site_var, site_var_t, site_list, Times, Y, tair_vals, Y_miss, Y_miss_t, error)
 
   IF (count(Y_miss) < ntimes) THEN
    ALLOCATE (TS(count(Y_miss)))
    ALLOCATE (XS(count(Y_miss), nvars+1))
    ALLOCATE (XP(count(Y_miss), nvars+1))
    ALLOCATE (TWXS(nvars+1, count(Y_miss)))
    ALLOCATE (YS(count(Y_miss)))
    TS (:) = pack (Times, Y_miss)
    DO v = 1, nvars + 1, 1
     XS (:, v) = pack (X(:, v), Y_miss)
    END DO
    XP (:, :) = XS (:, :)
    YS (:) = pack (Y, Y_miss)
   ELSE
    ALLOCATE (TS(ntimes))
    ALLOCATE (XS(ntimes, nvars+1))
    ALLOCATE (XP(ntimes, nvars+1))
    ALLOCATE (TWXS(nvars+1, ntimes))
    ALLOCATE (YS(ntimes))
    TS (:) = Times (:)
    XS (:, :) = X (:, :)
    XP (:, :) = XS (:, :)
    YS (:) = Y (:)
   END IF
 
   CALL normalize_X (XS)
   CALL normalize_Y (transform_exp, YS)
 
   tt = 1
   DO t = 1, ntimes, 1
 
    IF (Y_miss(t) .EQV. .TRUE.) THEN
 
     CALL calc_weights (TS, tt, XS, W)
     TWXS = matmul (transpose(XS), W)
     IF (trim(stnvar) .EQ. "PRCP") THEN
 
      CALL logistic_regressionrf (XP, YS, TWXS, B)
 
      POC (i, t, :) = B (:)
      DEALLOCATE (B)
     END IF
 
     CALL least_squares (XP, YS, TWXS, B)
     C (i, t, :) = B (:)
 
     DEALLOCATE (B)
 
     tt = tt + 1
    ELSE
     IF (trim(stnvar) .EQ. "PRCP") THEN
      POC (i, t, :) = - 999.99
     END IF
     C (i, t, :) = - 999.99
    END IF
   END DO
 
   DEALLOCATE (YS)
   DEALLOCATE (TWXS)
   DEALLOCATE (XS)
   DEALLOCATE (TS)
   DEALLOCATE (XP)
 
  END DO
 
END SUBROUTINE estimate_coefficients
 
 
SUBROUTINE estimate_precip (X, Z, nsta, ngrid, maxDistance, Times, stnid, stnvar, site_var, site_var_t, site_list, PCP, &
& POP, PCPERR, tmean, tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, Y_mean, Y_std, y_std_all, y_min, &
& y_max, error, PCP_2, POP_2, PCPERR_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
 
  USE strings
  USE utim
  USE type
  IMPLICIT NONE
 
  INTERFACE
 
   SUBROUTINE read_station (stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss, &
  & vals_miss_t, error)
    USE type
    CHARACTER (LEN=100), INTENT (IN) :: stnvar
    CHARACTER (LEN=100), INTENT (IN) :: stnid
    CHARACTER (LEN=100), INTENT (IN) :: site_var
    CHARACTER (LEN=100), INTENT (IN) :: site_var_t
    CHARACTER (LEN=500), INTENT (IN) :: site_list
    REAL (DP), INTENT (IN) :: Times (:)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: vals (:), tair_vals (:, :)
    LOGICAL, ALLOCATABLE, INTENT (OUT) :: vals_miss (:), vals_miss_t (:)
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_station
 
   SUBROUTINE normalize_X (X)
    USE type
    REAL (DP), INTENT (INOUT) :: X (:, :)
   END SUBROUTINE normalize_X
 
   SUBROUTINE normalize_Xv (X, weight, mean, stdev, stdev_all, smin, smax, Yp)
    USE type
    REAL (DP), INTENT (INOUT) :: X (:)
    REAL (DP), INTENT (IN) :: weight (:)
    REAL (DP), INTENT (OUT) :: mean
    REAL (DP), INTENT (OUT) :: stdev
    REAL (DP), INTENT (OUT) :: stdev_all
    REAL (DP), INTENT (OUT) :: smin
    REAL (DP), INTENT (OUT) :: smax
    INTEGER (I4B), INTENT (OUT) :: Yp (:)
   END SUBROUTINE normalize_Xv
 
   SUBROUTINE normalize_Y (texp, Y)
    USE type
    REAL (DP), INTENT (IN) :: texp !transform exponent
    REAL (DP), INTENT (INOUT) :: Y (:)
   END SUBROUTINE normalize_Y
 
   SUBROUTINE calc_weights (Times, tt, X, W)
    USE type
    REAL (DP), INTENT (IN) :: Times (:)
    INTEGER (I4B), INTENT (IN) :: tt
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: W (:, :)
   END SUBROUTINE calc_weights
 
   SUBROUTINE least_squares (X, Y, TX, B)
    USE type
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), INTENT (IN) :: Y (:)
    REAL (DP), INTENT (IN) :: TX (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
   END SUBROUTINE least_squares
 
   SUBROUTINE logistic_regression (X, Y, TX, Yp, B)
    USE type
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), INTENT (IN) :: Y (:)
    REAL (DP), INTENT (IN) :: TX (:, :)
    INTEGER (I4B), INTENT (IN) :: Yp (:)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
   END SUBROUTINE logistic_regression
 
   SUBROUTINE generic_corr (stn_data, tair_data, lag, window, auto_corr, t_p_corr)
    USE type
 
    REAL (DP), INTENT (IN) :: stn_data (:)
    REAL (DP), INTENT (IN) :: tair_data (:, :)
    INTEGER (I4B), INTENT (IN) :: lag
    INTEGER (I4B), INTENT (IN) :: window
    REAL (DP), INTENT (OUT) :: auto_corr
    REAL (DP), INTENT (OUT) :: t_p_corr
 
   END SUBROUTINE generic_corr
 
   SUBROUTINE calc_distance (lat1, lon1, lat2, lon2, dist)
    USE type
    IMPLICIT NONE
 
    REAL (DP), INTENT (IN) :: lat1, lon1, lat2, lon2
    REAL (DP), INTENT (OUT) :: dist
 
   END SUBROUTINE calc_distance
 
 
   SUBROUTINE heapsort (n, ra, rn)
    USE type
    IMPLICIT NONE
 
    INTEGER (I4B), INTENT (IN) :: n
    INTEGER (I4B), DIMENSION (:), INTENT (INOUT) :: rn
    REAL (DP), DIMENSION (:), INTENT (INOUT) :: ra
 
   END SUBROUTINE heapsort
 
  END INTERFACE
 
  REAL (DP), INTENT (IN) :: X (:, :), Z (:, :)!station and grid point description arrays
  REAL (DP), INTENT (IN) :: maxDistance !max distance for weight function
  INTEGER (I4B), INTENT (IN) :: nsta, ngrid !nuber of input stations and grid points
  REAL (DP), INTENT (IN) :: Times (:)!time step array
  CHARACTER (LEN=100), INTENT (IN) :: stnid (:)!station id array
  CHARACTER (LEN=100), INTENT (IN) :: stnvar, site_var, site_var_t !control file variables
  CHARACTER (LEN=500), INTENT (IN) :: site_list !file name of station list
  REAL (SP), ALLOCATABLE, INTENT (OUT) :: PCP (:, :), POP (:, :), PCPERR (:, :)!output variables for precipitation
  REAL (SP), ALLOCATABLE, INTENT (OUT) :: tmean (:, :), tmean_err (:, :)!OLS tmean estimate and error
  REAL (SP), ALLOCATABLE, INTENT (OUT) :: trange (:, :), trange_err (:, :)!OLS trange estimate and error
 
  REAL (SP), ALLOCATABLE, INTENT (OUT) :: tmean_2 (:, :), tmean_err_2 (:, :)!OLS tmean estimate and error
  REAL (SP), ALLOCATABLE, INTENT (OUT) :: trange_2 (:, :), trange_err_2 (:, :)!OLS trange estimate and error
  REAL (SP), ALLOCATABLE, INTENT (OUT) :: PCP_2 (:, :), POP_2 (:, :), PCPERR_2 (:, :)
 
 
  INTEGER, INTENT (OUT) :: error !integer error flag
  REAL (DP), INTENT (OUT) :: mean_autocorr (:)!mean autocorrelation from all stations over entire time period
  REAL (DP), INTENT (OUT) :: mean_tp_corr (:)!mean correlation for mean temp and precip
 
   !vary at each grid point and time step
  REAL (DP), INTENT (OUT) :: Y_mean (:, :), Y_std (:, :)!std and mean of time step precipitation
  REAL (DP), INTENT (OUT) :: y_std_all (:, :)!std of time step precip including stations with zero precip
  REAL (DP), INTENT (OUT) :: y_min (:, :), y_max (:, :)!min & max  of normalized time step precipitation
 
  REAL (DP), ALLOCATABLE :: Y (:), TWX (:, :), B (:), TX (:, :)
 
  REAL (DP), ALLOCATABLE :: Y_red (:), TWX_red (:, :), TX_red (:, :)!reduced matricies
  REAL (DP), ALLOCATABLE :: X_red (:, :)!reduced matricies
 
  REAL (DP), ALLOCATABLE :: TWX_red_t (:, :), TX_red_t (:, :)!reduced matricies
  REAL (DP), ALLOCATABLE :: X_red_t (:, :)!reduced matricies
 
  REAL (DP), ALLOCATABLE :: TWX_red_tr (:, :), TX_red_tr (:, :)!reduced matricies
  REAL (DP), ALLOCATABLE :: X_red_tr (:, :)!reduced matricies
 
  REAL (DP), ALLOCATABLE :: w_base (:, :)!initial distance weight matrix
  REAL (DP), ALLOCATABLE :: w_pcp_1d (:), w_temp_1d (:)
  INTEGER (I4B), ALLOCATABLE :: w_pcp_1d_loc (:), w_temp_1d_loc (:)
!  real(DP), allocatable :: w_pcp(:,:), w_temp(:,:) !distance weight matrices for a specific grid point
  REAL (DP), ALLOCATABLE :: w_pcp_red (:, :), w_temp_red (:, :)!reduced distance weigth matricies
 
 
  REAL (DP), ALLOCATABLE :: Y_tmean (:), Y_trange (:)!transformed station data arrays
  REAL (DP), ALLOCATABLE :: Y_tmean_red (:), Y_trange_red (:)!transformed station data arrays
  REAL (DP), ALLOCATABLE :: stn_vals (:), stn_data (:, :), tair_data (:, :, :), stn_tair (:, :)!original station data arrays
  REAL (DP), ALLOCATABLE :: auto_corr (:)!lag-1 autocorrelation for stations over entire time period used
  REAL (DP), ALLOCATABLE :: t_p_corr (:)!correlation between temp and precip
  INTEGER (I4B), ALLOCATABLE :: Yp (:)!binary for logistic regression
  INTEGER (I4B), ALLOCATABLE :: Yp_red (:)!reduced binary for logistic regression
 
 
  LOGICAL, ALLOCATABLE :: stn_miss (:), stn_miss_t (:)!missing value logical arrays
 
  REAL (DP) :: m_tmean, m_trange !used to fill missing values in stations
 
  REAL (DP) :: errsum, wgtsum, p, sta_pcp, sta_temp
  REAL (DP) :: auto_corr_sum, tp_corr_sum
  REAL (DP) :: step_mean, step_std, step_std_all, step_min, step_max !timestep statistics
 
  REAL (DP) :: rsqr, ss_tot, ss_res, vif !r-squared and variance correction
 
  INTEGER (I4B) :: xsize !size of second dimension of input X array
 
  INTEGER (I4B) :: ntimes, nstns
  INTEGER (I4B) :: t, i, g, ndata, nodata
  INTEGER (I4B) :: ndata_t, nodata_t
  INTEGER (I4B) :: lag, window, tc, trc
  INTEGER (I4B) :: auto_cnt, tp_cnt
 
  INTEGER (I4B) :: stn_count
 
 
  !variables for tracking closest N stations for precipitation
  INTEGER (I4B) :: out_loc
  INTEGER (I4B), PARAMETER :: sta_limit = 30
  INTEGER (I4B), ALLOCATABLE :: close_loc (:, :)
  INTEGER (I4B), ALLOCATABLE :: close_count (:)
 
  REAL (DP), ALLOCATABLE :: close_weights (:, :)
  REAL (DP), ALLOCATABLE :: close_meta (:, :, :)
  REAL (DP) :: min_weight
  REAL (DP) :: max_distance
  REAL (DP), PARAMETER :: search_distance = 1000.0
 
  !variables for tracking closest N stations for temperature
  INTEGER (I4B) :: out_loc_t
  INTEGER (I4B), ALLOCATABLE :: close_loc_t (:, :)
  INTEGER (I4B), ALLOCATABLE :: close_count_t (:)
 
  REAL (DP), ALLOCATABLE :: close_weights_t (:, :)
  REAL (DP), ALLOCATABLE :: close_meta_t (:, :, :)
  REAL (DP) :: min_weight_t
  REAL (DP) :: max_distance_t
 
 
  REAL (DP) :: tmp_pcp
  REAL (DP) :: tmp_weight
 
  INTEGER (I4B) :: slope_flag_pcp
  INTEGER (I4B) :: slope_flag_temp
 
!variables to check for singular matrix
  REAL (DP), ALLOCATABLE :: tmp (:, :)
  REAL (DP), ALLOCATABLE :: vv (:)
 
!variables for timing code
  INTEGER (I4B) :: t1, t2, count_rate
  INTEGER (I4B) :: tg1, tg2
 
 
!variables for keeping track of max_distance modifications
  INTEGER (I4B), ALLOCATABLE :: expand_flag (:), expand_flag_t (:)
  REAL (DP), ALLOCATABLE :: expand_dist (:), expand_dist_t (:)
 
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
  ALLOCATE (Y(nstns))
  ALLOCATE (Y_tmean(nstns), Y_trange(nstns))
  ALLOCATE (stn_data(nstns, ntimes))
  ALLOCATE (tair_data(2, nstns, ntimes))
 
!original allocations
!  allocate(w_pcp(nstns,nstns))
!  allocate(w_temp(nstns,nstns))
!  allocate(TWX(4,nstns))
!  allocate(TX(4,nstns))
 
  ALLOCATE (w_pcp_red(sta_limit, sta_limit))
  ALLOCATE (w_temp_red(sta_limit, sta_limit))
  ALLOCATE (Y_red(sta_limit))
  ALLOCATE (Y_tmean_red(sta_limit), Y_trange_red(sta_limit))
  ALLOCATE (X_red(sta_limit, xsize))
  ALLOCATE (X_red_t(sta_limit, xsize))
 
  ALLOCATE (w_pcp_1d(sta_limit))
  ALLOCATE (w_temp_1d(sta_limit))
  ALLOCATE (w_pcp_1d_loc(sta_limit))
  ALLOCATE (w_temp_1d_loc(sta_limit))
 
  ALLOCATE (tmp(6, 6))
  ALLOCATE (vv(6))
 
  ALLOCATE (PCP(ngrid, ntimes))
  ALLOCATE (POP(ngrid, ntimes))
  ALLOCATE (PCPERR(ngrid, ntimes))
 
  ALLOCATE (tmean(ngrid, ntimes))
  ALLOCATE (tmean_err(ngrid, ntimes))
  ALLOCATE (trange(ngrid, ntimes))
  ALLOCATE (trange_err(ngrid, ntimes))
 
 
  ALLOCATE (PCP_2(ngrid, ntimes))
  ALLOCATE (POP_2(ngrid, ntimes))
  ALLOCATE (PCPERR_2(ngrid, ntimes))
 
  ALLOCATE (tmean_2(ngrid, ntimes))
  ALLOCATE (tmean_err_2(ngrid, ntimes))
  ALLOCATE (trange_2(ngrid, ntimes))
  ALLOCATE (trange_err_2(ngrid, ntimes))
 
 
  ALLOCATE (auto_corr(nstns))
  ALLOCATE (t_p_corr(nstns))
 
  ALLOCATE (Yp(nstns))
  ALLOCATE (Yp_red(sta_limit))
 
  !station limit arrays
  ALLOCATE (close_weights(ngrid, sta_limit))
  ALLOCATE (close_loc(ngrid, sta_limit))
  ALLOCATE (close_meta(5, ngrid, sta_limit))
  ALLOCATE (close_count(ngrid))
 
  ALLOCATE (close_weights_t(ngrid, sta_limit))
  ALLOCATE (close_loc_t(ngrid, sta_limit))
  ALLOCATE (close_meta_t(5, ngrid, sta_limit))
  ALLOCATE (close_count_t(ngrid))
 
 
  !base weight array
  ALLOCATE (w_base(ngrid, nstns))
 
  !max_dist tracking variables
  ALLOCATE (expand_dist(ngrid), expand_flag(ngrid))
  ALLOCATE (expand_dist_t(ngrid), expand_flag_t(ngrid))
 
 
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
 
  DO i = 1, nstns, 1
 
   CALL read_station (stnvar, stnid(i), site_var, site_var_t, site_list, Times, stn_vals, stn_tair, stn_miss, &
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
   CALL generic_corr (stn_data(i, :), tair_data(:, i, :), lag, window, auto_corr(i), t_p_corr(i))
 
    !compute mean autocorrelation for all stations and all times
    !check for values outside of -1 to 1
    !stations with incomplete data are set to -999
   IF (auto_corr(i) .GE.-1.0 .AND. auto_corr(i) .LE. 1.0) THEN
    auto_corr_sum = auto_corr_sum + auto_corr (i)
    auto_cnt = auto_cnt + 1
   END IF
   IF (t_p_corr(i) .GE.-1.0 .AND. t_p_corr(i) .LE. 1.0) THEN
    tp_corr_sum = tp_corr_sum + t_p_corr (i)
    tp_cnt = tp_cnt + 1
   END IF
 
   DEALLOCATE (stn_miss_t)
   DEALLOCATE (stn_miss)
   DEALLOCATE (stn_vals)
   DEALLOCATE (stn_tair)
  END DO !end station read loop
 
 
  error = 0
 
  mean_autocorr = auto_corr_sum / real (auto_cnt, kind(DP))
  mean_tp_corr = tp_corr_sum / real (tp_cnt, kind(DP))
 
  PRINT *, 'Temp lag-1 autocorrelation: ', mean_autocorr (1)
  PRINT *, 'Temp-precip correlation: ', mean_tp_corr (1)
 
 
  !pull weight generation outside of time loop
 
 
  PRINT *, 'Generating base weight matrix & '
  PRINT *, 'finding nearest stations for each gridpoint'
  CALL system_clock (t1, count_rate)
 
 
  DO g = 1, ngrid, 1
 
 
   close_count (g) = 1
   min_weight = 0.0d0
   close_weights (g, :) = 0.0d0
 
   close_count_t (g) = 1
   min_weight_t = 0.0d0
   close_weights_t (g, :) = 0.0d0
 
   DO i = 1, nstns, 1
      !setup distinct weight matrices for precip and temperature
    CALL calc_distance_weight (search_distance, X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), w_base(g, i))
 
!original call
!      call calc_distance_weight(maxDistance, X(i,2), X(i,3), Z(g,2), Z(g,3), w_temp(i,i))
 
  !also set some logic to limit the number of stations to the N closest
    min_weight = 0.0d0
 
    IF (w_base(g, i) .GT. min_weight .AND. stn_data(i, 1) .GT.-100.0d0) THEN
     IF (close_count(g) .LE. sta_limit) THEN
 
      close_weights (g, close_count(g)) = w_base (g, i)
      close_loc (g, close_count(g)) = i
 
      close_meta (1, g, close_count(g)) = X (i, 2)
      close_meta (2, g, close_count(g)) = X (i, 3)
      close_meta (3, g, close_count(g)) = Z (g, 2)
      close_meta (4, g, close_count(g)) = Z (g, 3)
      CALL calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta(5, g, close_count(g)))
 
      close_count (g) = close_count (g) + 1
     ELSE
      min_weight = minval (close_weights(g, :), 1)
      IF (w_base(g, i) .GT. min_weight) THEN
       out_loc = minloc (close_weights(g, :), 1)
       close_weights (g, out_loc) = w_base (g, i)
       close_loc (g, out_loc) = i
 
       close_meta (1, g, out_loc) = X (i, 2)
       close_meta (2, g, out_loc) = X (i, 3)
       close_meta (3, g, out_loc) = Z (g, 2)
       close_meta (4, g, out_loc) = Z (g, 3)
       CALL calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta(5, g, out_loc))
      END IF
 
     END IF
    END IF
 
  !need to repeat above for temperature since that data is independent of precipitation
    min_weight_t = 0.0d0
 
    IF (w_base(g, i) .GT. min_weight_t .AND. tair_data(1, i, 1) .GT.-200.0d0) THEN
     IF (close_count_t(g) .LE. sta_limit) THEN
 
      close_weights_t (g, close_count_t(g)) = w_base (g, i)
      close_loc_t (g, close_count_t(g)) = i
 
      close_meta_t (1, g, close_count_t(g)) = X (i, 2)
      close_meta_t (2, g, close_count_t(g)) = X (i, 3)
      close_meta_t (3, g, close_count_t(g)) = Z (g, 2)
      close_meta_t (4, g, close_count_t(g)) = Z (g, 3)
      CALL calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta_t(5, g, close_count_t(g)))
 
      close_count_t (g) = close_count_t (g) + 1
     ELSE
      min_weight_t = minval (close_weights(g, :), 1)
      IF (w_base(g, i) .GT. min_weight) THEN
       out_loc_t = minloc (close_weights_t(g, :), 1)
       close_weights_t (g, out_loc_t) = w_base (g, i)
 
       close_loc_t (g, out_loc_t) = i
 
       close_meta_t (1, g, out_loc_t) = X (i, 2)
       close_meta_t (2, g, out_loc_t) = X (i, 3)
       close_meta_t (3, g, out_loc_t) = Z (g, 2)
       close_meta_t (4, g, out_loc_t) = Z (g, 3)
       CALL calc_distance (X(i, 2), X(i, 3), Z(g, 2), Z(g, 3), close_meta(5, g, out_loc_t))
      END IF
     END IF
    END IF
 
   END DO !end station loop
 
 
  END DO !end grid point loop
  CALL system_clock (t2, count_rate)
  PRINT *, 'Elapsed time for weight generation: ', real (t2-t1) / real (count_rate)
 
  DO t = 1, ntimes, 1
 
   CALL system_clock (tg1, count_rate)
   PRINT *, "TIME STEP = ", Times (t), " (", t, "/", ntimes, ")"
 
 
   DO i = 1, nstns, 1
 
    Y (i) = stn_data (i, t)
 
    Y_tmean (i) = (tair_data(1, i, t)+tair_data(2, i, t)) / 2.0d0
    Y_trange (i) = Abs (tair_data(2, i, t)-tair_data(1, i, t))
 
   END DO
 
!do power transformation on input vector
   CALL normalize_Y (4.0d0, Y)
 
 
   DO g = 1, ngrid, 1
 
! call system_clock(tg1,count_rate)
 
    IF (Z(g, 4) .GT.-200) THEN
 
     ALLOCATE (TWX_red(6, sta_limit))
     ALLOCATE (TX_red(6, sta_limit))
     ALLOCATE (TWX_red_t(6, sta_limit))
     ALLOCATE (TX_red_t(6, sta_limit))
 
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
     DO i = 1, (close_count(g)-1)
      IF (close_meta(5, g, i) .GT. max_distance) THEN
       max_distance = close_meta (5, g, i)
      END IF
     END DO
 
 
 
     IF (max_distance .LE. maxDistance) THEN
      expand_dist (g) = max_distance
      expand_flag (g) = 0
      max_distance = maxDistance
     ELSE
      max_distance = max_distance + 1.0d0
      expand_flag (g) = 1
      expand_dist (g) = max_distance
     END IF
 
    !reduced matrices for precip
     slope_flag_pcp = 0
     DO i = 1, (close_count(g)-1)
      CALL calc_distance_weight (max_distance, close_meta(1, g, i), close_meta(2, g, i), close_meta(3, g, i), &
     & close_meta(4, g, i), tmp_weight)
 
      w_pcp_red (i, i) = tmp_weight
      w_pcp_1d (i) = tmp_weight
      w_pcp_1d_loc (i) = close_loc (g, i)
      Y_red (i) = Y (close_loc(g, i))
      X_red (i, :) = X (close_loc(g, i), :)
 
 
 
      IF (stn_data(close_loc(g, i), t) .GT. 0.0) THEN
       ndata = ndata + 1
       Yp_red (i) = 1
      ELSE
       nodata = nodata + 1
      END IF
 
 
     END DO
 
 
     CALL normalize_Xv (Y_red, w_pcp_1d, step_mean, step_std, step_std_all, step_min, step_max, Yp_red)
 
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
     DO i = 1, (close_count_t(g)-1)
      IF (close_meta_t(5, g, i) .GT. max_distance_t) THEN
       max_distance_t = close_meta_t (5, g, i)
      END IF
     END DO
 
     IF (max_distance_t .LE. maxDistance) THEN
      max_distance_t = maxDistance
     ELSE
      max_distance_t = max_distance_t + 1.0d0
      expand_flag_t (g) = 1
      expand_dist_t (g) = max_distance_t
     END IF
 
    !reduced matrices for temperature
     slope_flag_temp = 0
     DO i = 1, (close_count_t(g)-1)
      CALL calc_distance_weight (max_distance_t, close_meta_t(1, g, i), close_meta_t(2, g, i), close_meta_t(3, g, i), &
     & close_meta_t(4, g, i), tmp_weight)
 
      w_temp_red (i, i) = tmp_weight
      w_temp_1d (i) = tmp_weight
      Y_tmean_red (i) = Y_tmean (close_loc_t(g, i))
      Y_trange_red (i) = Y_trange (close_loc_t(g, i))
      X_red_t (i, :) = X (close_loc_t(g, i), :)
 
      IF (Y_tmean(close_loc_t(g, i)) .GT.-100.0) THEN
       ndata_t = ndata_t + 1
      ELSE
       nodata_t = nodata_t + 1
      END IF
 
 
     END DO
 
     IF (ndata == 0 .AND. nodata == 0) THEN
     !print *, "No stations within max distance of grid cell!"
      POP (g, t) = 0.0
      PCP (g, t) = 0.0
      PCPERR (g, t) = 0.0
 
      POP_2 (g, t) = 0.0
      PCP_2 (g, t) = 0.0
      PCPERR_2 (g, t) = 0.0
 
     END IF
 
     IF (ndata_t == 0 .AND. nodata_t == 0) THEN
      IF (t .GT. 1) THEN
       tmean (g, t) = tmean (g, t-1)
       trange (g, t) = trange (g, t-1)
       tmean_err (g, t) = tmean_err (g, t-1)
       trange_err (g, t) = trange_err (g, t-1)
 
       tmean_2 (g, t) = tmean_2 (g, t-1)
       trange_2 (g, t) = trange_2 (g, t-1)
       tmean_err_2 (g, t) = tmean_err_2 (g, t-1)
       trange_err_2 (g, t) = trange_err_2 (g, t-1)
      ELSE
       tmean (g, t) = - 999
       trange (g, t) = - 999
       tmean_err (g, t) = 0.0
       trange_err (g, t) = 0.0
 
       tmean_2 (g, t) = - 999
       trange_2 (g, t) = - 999
       tmean_err_2 (g, t) = 0.0
       trange_err_2 (g, t) = 0.0
      END IF
     END IF
 
     IF (ndata >= 1) THEN
 
            !tmp needs to be matmul(TX,X) where TX = TWX_red and X = X_red
      TWX_red = matmul (transpose(X_red), w_pcp_red)
      tmp = matmul (TWX_red, X_red)
      vv = maxval (Abs(tmp), dim=2)
 
      IF (any(vv == 0.0)) THEN
       slope_flag_pcp = 0
      ELSE
       slope_flag_pcp = 1
      END IF
 
      IF (nodata == 0) THEN
       !print *, "All stations have precip, POP = 1.0"
       POP (g, t) = 1.0
 
       POP_2 (g, t) = 1.0
 
      ELSE
 
 
 
       IF (slope_flag_pcp .EQ. 0) THEN
        POP (g, t) = - 999.
       ELSE
       !regression with slope
        TX_red = transpose (X_red)
        TWX_red = matmul (TX_red, w_pcp_red)
 
        CALL logistic_regression (X_red, Y_red, TWX_red, Yp_red, B)
        POP (g, t) = real (1.0/(1.0+Exp(-dot_product(Z(g, :), B))), kind(SP))
 
        DEALLOCATE (B)
       END IF
 
 
       !regression without slope
       DEALLOCATE (TWX_red)
       DEALLOCATE (TX_red)
       ALLOCATE (TWX_red(4, sta_limit))
       ALLOCATE (TX_red(4, sta_limit))
       TX_red = transpose (X_red(:, 1:4))
 
       TWX_red = matmul (TX_red, w_pcp_red)
 
 
 
       CALL logistic_regression (X_red(:, 1:4), Y_red, TWX_red, Yp_red, B)
       POP_2 (g, t) = real (1.0/(1.0+Exp(-dot_product(Z(g, 1:4), B))), kind(SP))
 
 
       DEALLOCATE (B)
 
      END IF
 
      DEALLOCATE (TWX_red)
      DEALLOCATE (TX_red)
      ALLOCATE (TWX_red(6, sta_limit))
      ALLOCATE (TX_red(6, sta_limit))
 
 
      IF (slope_flag_pcp .EQ. 0) THEN
       PCP (g, t) = - 999.
      ELSE
       !regression with slope
       TX_red = transpose (X_red)
       TWX_red = matmul (TX_red, w_pcp_red)
 
       CALL least_squares (X_red, Y_red, TWX_red, B)
       PCP (g, t) = real (dot_product(Z(g, :), B), kind(SP))
 
 
       wgtsum = 0.0
       errsum = 0.0
       ss_tot = 0.0
       ss_res = 0.0
       DO i = 1, (close_count(g)-1), 1
        wgtsum = wgtsum + w_pcp_red (i, i)
 
 
        errsum = errsum + (w_pcp_red(i, i)*(PCP(g, t)-Y_red(i))**2)
 
       END DO
       PCPERR (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
      END IF
 
 
      !regression without slope
 
      DEALLOCATE (TWX_red)
      DEALLOCATE (TX_red)
      ALLOCATE (TWX_red(4, sta_limit))
      ALLOCATE (TX_red(4, sta_limit))
 
      TX_red = transpose (X_red(:, 1:4))
      TWX_red = matmul (TX_red, w_pcp_red)
 
      CALL least_squares (X_red(:, 1:4), Y_red, TWX_red, B)
      PCP_2 (g, t) = real (dot_product(Z(g, 1:4), B), kind(SP))
 
      wgtsum = 0.0
      errsum = 0.0
      ss_tot = 0.0
      ss_res = 0.0
      DO i = 1, (close_count(g)-1), 1
       wgtsum = wgtsum + w_pcp_red (i, i)
 
       errsum = errsum + (w_pcp_red(i, i)*(PCP_2(g, t)-Y_red(i))**2)
 
      END DO
      PCPERR_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
 
      DEALLOCATE (B)
 
     ELSE
     !print *, "Not enough stations with data within max distance"
      POP (g, t) = 0.0
      PCP (g, t) = 0.0
      PCPERR (g, t) = 0.0
 
      POP_2 (g, t) = 0.0
      PCP_2 (g, t) = 0.0
      PCPERR_2 (g, t) = 0.0
     END IF !precip if statement
 
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 !  do temperature ols now
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
     IF (ndata_t .GE. 1) THEN
 
 
   !regression with slope
      TX_red_t = transpose (X_red_t)
      TWX_red_t = matmul (TX_red_t, w_temp_red)
      CALL least_squares (X_red_t, Y_tmean_red, TWX_red_t, B)
      tmean (g, t) = real (dot_product(Z(g, :), B), kind(SP))
 
      errsum = 0.0
      wgtsum = 0.0
      DO i = 1, (close_count_t(g)-1), 1
       wgtsum = wgtsum + w_temp_red (i, i)
 
       errsum = errsum + (w_temp_red(i, i)*(tmean(g, t)-Y_tmean_red(i))**2)
      END DO
      tmean_err (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
 
   !regression without slope
      DEALLOCATE (B)
      DEALLOCATE (TWX_red_t)
      DEALLOCATE (TX_red_t)
      ALLOCATE (TWX_red_t(4, sta_limit))
      ALLOCATE (TX_red_t(4, sta_limit))
      TX_red_t = transpose (X_red_t(:, 1:4))
      TWX_red_t = matmul (TX_red_t, w_temp_red)
 
      CALL least_squares (X_red_t(:, 1:4), Y_tmean_red, TWX_red_t, B)
      tmean_2 (g, t) = real (dot_product(Z(g, 1:4), B), kind(SP))
 
      errsum = 0.0
      wgtsum = 0.0
   !do i = 1, nstns, 1
      DO i = 1, (close_count_t(g)-1), 1
       wgtsum = wgtsum + w_temp_red (i, i)
 
 
       errsum = errsum + (w_temp_red(i, i)*(tmean_2(g, t)-Y_tmean_red(i))**2)
 
      END DO
      tmean_err_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
      DEALLOCATE (B)
 
   !then do trange
 
      DEALLOCATE (TWX_red_t)
      DEALLOCATE (TX_red_t)
      ALLOCATE (TWX_red_t(6, sta_limit))
      ALLOCATE (TX_red_t(6, sta_limit))
      TX_red_t = transpose (X_red_t)
      TWX_red_t = matmul (TX_red_t, w_temp_red)
 
      CALL least_squares (X_red_t, Y_trange_red, TWX_red_t, B)
      trange (g, t) = real (dot_product(Z(g, :), B), kind(SP))
 
      errsum = 0.0
      wgtsum = 0.0
      DO i = 1, (close_count_t(g)-1), 1
       wgtsum = wgtsum + w_temp_red (i, i)
 
       errsum = errsum + (w_temp_red(i, i)*(trange(g, t)-Y_trange_red(i))**2)
 
      END DO
      trange_err (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
          !regression without slope
      DEALLOCATE (B)
      DEALLOCATE (TWX_red_t)
      DEALLOCATE (TX_red_t)
      ALLOCATE (TWX_red_t(4, sta_limit))
      ALLOCATE (TX_red_t(4, sta_limit))
      TX_red_t = transpose (X_red_t(:, 1:4))
      TWX_red_t = matmul (TX_red_t, w_temp_red)
 
      CALL least_squares (X_red_t(:, 1:4), Y_trange_red, TWX_red_t, B)
      trange_2 (g, t) = real (dot_product(Z(g, 1:4), B), kind(SP))
 
      errsum = 0.0
      wgtsum = 0.0
      DO i = 1, (close_count_t(g)-1), 1
       wgtsum = wgtsum + w_temp_red (i, i)
 
       sta_temp = real (dot_product(X_red_t(i, 1:4), B), kind(SP))
 
       errsum = errsum + (w_temp_red(i, i)*(trange_2(g, t)-Y_trange_red(i))**2)
 
      END DO
      trange_err_2 (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(SP))
 
     ELSE
 
     !if not enough stations with data
     !just use value from previous grid point for now
      IF (g .GT. 1) THEN
       trange (g, t) = trange (g-1, t)
       trange_err (g, t) = trange_err (g-1, t)
       tmean (g, t) = tmean (g-1, t)
       tmean_err (g, t) = tmean_err (g-1, t)
 
       trange_2 (g, t) = trange_2 (g-1, t)
       trange_err_2 (g, t) = trange_err_2 (g-1, t)
       tmean_2 (g, t) = tmean_2 (g-1, t)
       tmean_err_2 (g, t) = tmean_err_2 (g-1, t)
      ELSE
       trange (g, t) = trange (g, t-1)
       trange_err (g, t) = trange_err (g-1, t-1)
       tmean (g, t) = tmean (g, t-1)
       tmean_err (g, t) = tmean_err (g, t-1)
 
       trange_2 (g, t) = trange_2 (g, t-1)
       trange_err_2 (g, t) = trange_err_2 (g-1, t-1)
       tmean_2 (g, t) = tmean_2 (g, t-1)
       tmean_err_2 (g, t) = tmean_err_2 (g, t-1)
      END IF
     END IF !end data check if statement for temperature
 
 
    END IF !end check for valid elevation
 
    IF (allocated(TWX_red)) THEN
     DEALLOCATE (TWX_red)
    END IF
    IF (allocated(TX_red)) THEN
     DEALLOCATE (TX_red)
    END IF
    IF (allocated(TWX_red_t)) THEN
     DEALLOCATE (TWX_red_t)
    END IF
    IF (allocated(TX_red_t)) THEN
     DEALLOCATE (TX_red_t)
    END IF
 
   END DO !end grid loop
 
   CALL system_clock (tg2, count_rate)
   PRINT *, 'Elapsed time for one time step: ', real (tg2-tg1) / real (count_rate)
  END DO !end time loop
 
END SUBROUTINE estimate_precip
 
 
SUBROUTINE normalize_X (X)
  USE type
  IMPLICIT NONE
 
  REAL (DP), INTENT (INOUT) :: X (:, :)
 
  REAL (DP) :: mean, stdev, sum_x, sum_x2
  INTEGER (I4B) :: v, t
  INTEGER (I4B) :: nvars, ntimes
 
  nvars = size (X, 2) - 1
  ntimes = size (X, 1)
 
  DO v = 2, nvars + 1, 1
   sum_x = 0.0d0
   sum_x2 = 0.0d0
   DO t = 1, ntimes, 1
    sum_x = sum_x + X (t, v)
    sum_x2 = sum_x2 + X (t, v) ** 2
   END DO
   mean = sum_x / real (ntimes)
   stdev = Sqrt ((real(ntimes)*sum_x2-sum_x**2)/(real(ntimes)*real(ntimes-1)))
   DO t = 1, ntimes, 1
    IF (stdev .EQ. 0.0) THEN
     X (t, v) = X (t, v)
    ELSE
     X (t, v) = (X(t, v)-mean) / stdev
    END IF
   END DO
  END DO
 
END SUBROUTINE normalize_X
 
 
SUBROUTINE normalize_Xv (X, weight, mean, stdev, stdev_all, smin, smax, Yp)
  USE type
  IMPLICIT NONE
 
  REAL (DP), INTENT (INOUT) :: X (:)
  REAL (DP), INTENT (IN) :: weight (:)
  REAL (DP), INTENT (OUT) :: mean
  REAL (DP), INTENT (OUT) :: stdev
  REAL (DP), INTENT (OUT) :: stdev_all
  REAL (DP), INTENT (OUT) :: smin
  REAL (DP), INTENT (OUT) :: smax
  INTEGER (I4B), INTENT (OUT) :: Yp (:)
 
!  real(DP) :: mean, stdev, sum_x2, sum_x
  REAL (DP) :: sum_x2, sum_x
  REAL (DP) :: sum_stdev, sum_std
 
  REAL (DP) :: mean_all, sum_weight, sum_xw
 
  INTEGER (I4B) :: t
  INTEGER (I4B) :: ntimes
 
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
  DO t = 1, ntimes, 1
   IF (X(t) > 0.0) THEN
    sum_x = sum_x + X (t)
    sum_x2 = sum_x2 + X (t) ** 2
    sum_xw = sum_xw + weight (t) * X (t)
    sum_weight = sum_weight + weight (t)
    Yp (t) = 1
 
    IF (X(t) .LE. smin) THEN
     smin = X (t)
    END IF
    IF (X(t) .GE. smax) THEN
     smax = X (t)
    END IF
 
   END IF
  END DO
 
  mean_all = sum_x / real (ntimes)
  mean = sum_x / real (sum(Yp))
  DO t = 1, ntimes, 1
   sum_stdev = sum_stdev + (X(t)-mean_all) ** 2
   IF (Yp(t) .EQ. 1) THEN
    sum_std = sum_std + (X(t)-mean) ** 2
   END IF
  END DO
!  stdev_all = sqrt(sum_stdev/(ntimes-1))
 
 
  IF (sum(Yp) .GE. 2) THEN
   stdev = Sqrt (sum_std/real(sum(Yp)-1.0))
   stdev_all = Sqrt (sum_std/real(ntimes-1.0))
 
  ELSE
   mean = sum_x
   stdev = 0.01
!    stdev_all = sum_xw
   stdev_all = 0.01
  END IF
 
 
  IF (stdev .GT. 0.0) THEN
   DO t = 1, ntimes, 1
    IF (Yp(t) .GT. 0.0) THEN
     X (t) = (X(t)-mean) / stdev
    END IF
   END DO
  END IF
 
 
  RETURN
 
END SUBROUTINE normalize_Xv
 
 
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
SUBROUTINE normalize_Y (texp, Y)
  USE type
  IMPLICIT NONE
 
  REAL (DP), INTENT (IN) :: texp !transform exponent
  REAL (DP), INTENT (INOUT) :: Y (:)
  INTEGER (I4B) :: t
  INTEGER (I4B) :: ntimes
 
  ntimes = size (Y)
 
  DO t = 1, ntimes, 1
!     Y(t) = Y(t) ** (1.0d0/2.5d0)
!    Y(t) = Y(t) ** (1.0d0/4d0)
   Y (t) = Y (t) ** (1.0d0/texp)
  END DO
 
END SUBROUTINE normalize_Y
 
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
SUBROUTINE calc_weights (Times, tt, X, W)
  USE type
  IMPLICIT NONE
 
  REAL (DP), INTENT (IN) :: Times (:)
  INTEGER (I4B), INTENT (IN) :: tt
  REAL (DP), INTENT (IN) :: X (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: W (:, :)
  REAL (DP) :: sum
  INTEGER (I4B) :: v, t
  INTEGER (I4B) :: nvars, ntimes
 
  nvars = size (X, 2) - 1
  ntimes = size (X, 1)
 
  ALLOCATE (W(ntimes, ntimes))
 
  DO t = 1, ntimes, 1
   W (t, :) = 0.0d0
   IF (t /= tt) THEN
    sum = 0.0d0
    DO v = 2, nvars + 1, 1
     sum = sum + (X(t, v)-X(tt, v)) ** 2
    END DO
    W (t, :) = 0.0d0
 
    W (t, t) = 1 / Sqrt (sum/nvars)
   END IF
  END DO
 
END SUBROUTINE calc_weights
 
! Great circle distance calculation
! Output in nm
SUBROUTINE calc_distance_weight (maxd, lat1, lon1, lat2, lon2, weight)
  USE type
  IMPLICIT NONE
 
  REAL (DP), INTENT (IN) :: maxd, lat1, lon1, lat2, lon2
  REAL (DP), INTENT (OUT) :: weight
 
  REAL (DP) :: dist, lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927
 
  lat1r = lat1 * Pi / 180.0d0
  lon1r = lon1 * Pi / 180.0d0
  lat2r = lat2 * Pi / 180.0d0
  lon2r = lon2 * Pi / 180.0d0
  dist = ((180*60)/Pi) * (2*Asin(Sqrt((Sin((lat1r-lat2r)/2))**2+Cos(lat1r)*Cos(lat2r)*(Sin((lon1r-lon2r)/2))**2)))
  IF (dist .GT. maxd) THEN
   weight = 0.0d0
  ELSE
   weight = (1.0d0-(dist/maxd)**3) ** 3
!    weight = 1.0d0 - (dist/maxd)**0.5
  END IF
 
END SUBROUTINE calc_distance_weight
 
! Great circle distance calculation
! Output in nm
SUBROUTINE calc_distance (lat1, lon1, lat2, lon2, dist)
  USE type
  IMPLICIT NONE
 
  REAL (DP), INTENT (IN) :: lat1, lon1, lat2, lon2
  REAL (DP), INTENT (OUT) :: dist
 
  REAL (DP) :: lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927
 
  lat1r = lat1 * Pi / 180.0d0
  lon1r = lon1 * Pi / 180.0d0
  lat2r = lat2 * Pi / 180.0d0
  lon2r = lon2 * Pi / 180.0d0
  dist = ((180*60)/Pi) * (2*Asin(Sqrt((Sin((lat1r-lat2r)/2))**2+Cos(lat1r)*Cos(lat2r)*(Sin((lon1r-lon2r)/2))**2)))
 
END SUBROUTINE calc_distance
 
!
! Solve linear equation for x (Ax = b => x = bA^-1) using LU decomposition and back substitution.
! Input:
!   X  = An m by n array.
!   TX = Precalculated transpose array of X, size n by m
!   Y  = An m-element vector containing the right-hand side of the linear system Ax = b.
! Output:
!   B  = An n-element vector.
SUBROUTINE least_squares (X, Y, TX, B)
  USE type
  IMPLICIT NONE
 
  INTERFACE
   SUBROUTINE ludcmp (a, indx, D)
    USE type
    REAL (DP), DIMENSION (:, :), INTENT (INOUT) :: a
    INTEGER (I4B), DIMENSION (:), INTENT (OUT) :: indx
    REAL (DP), INTENT (OUT) :: D
   END SUBROUTINE ludcmp
 
   SUBROUTINE lubksb (a, indx, B)
    USE type
    REAL (DP), DIMENSION (:, :), INTENT (IN) :: a
    INTEGER (I4B), DIMENSION (:), INTENT (IN) :: indx
    REAL (DP), DIMENSION (:), INTENT (INOUT) :: B
   END SUBROUTINE lubksb
  END INTERFACE
 
  REAL (DP), INTENT (IN) :: X (:, :)
  REAL (DP), INTENT (IN) :: Y (:)
  REAL (DP), INTENT (IN) :: TX (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
 
  REAL (DP), ALLOCATABLE :: a (:, :)
  INTEGER (I4B), ALLOCATABLE :: indx (:)
  INTEGER (I4B) :: nvars, ntimes
  REAL (DP) :: D
 
  nvars = size (X, 2) - 1
  ntimes = size (Y)
 
  ALLOCATE (B(nvars+1))
  ALLOCATE (a(nvars+1, nvars+1))
  ALLOCATE (indx(nvars+1))
 
  B = matmul (TX, Y)
  a = matmul (TX, X)
 
  CALL ludcmp (a, indx, D)
  IF (any(Abs(a) < 9.99999968E-15)) THEN
   B (:) = 0.0d0
   PRINT *, "Warning, LUdcmp produced a zero."
   RETURN
  END IF
 
  CALL lubksb (a, indx, B)
 
  DEALLOCATE (a)
  DEALLOCATE (indx)
 
END SUBROUTINE least_squares
 
 
SUBROUTINE logistic_regression (X, Y, TX, Yp, B)
  USE type
  IMPLICIT NONE
 
  INTERFACE
   SUBROUTINE least_squares (X, Y, TX, B)
    USE type
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), INTENT (IN) :: Y (:)
    REAL (DP), INTENT (IN) :: TX (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
   END SUBROUTINE least_squares
  END INTERFACE
 
  REAL (DP), INTENT (IN) :: X (:, :)
  REAL (DP), INTENT (IN) :: Y (:)
  REAL (DP), INTENT (IN) :: TX (:, :)
  INTEGER (I4B), INTENT (IN) :: Yp (:)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
 
  REAL (DP), ALLOCATABLE :: Ypd (:), p (:), YN (:), BN (:), v (:, :), XV (:, :)
!  real(DP), allocatable :: P(:), YN(:), BN(:), V(:,:), XV(:,:)
  INTEGER (I4B) :: nvars, ntimes, i, t, f, it
  REAL (DP) :: D
 
  nvars = size (X, 2) - 1
  ntimes = size (Y)
 
  ALLOCATE (B(nvars+1))
  ALLOCATE (Ypd(ntimes))
  ALLOCATE (YN(ntimes))
  ALLOCATE (p(ntimes))
  ALLOCATE (v(ntimes, ntimes))
  ALLOCATE (XV(ntimes, nvars+1))
 
  DO t = 1, ntimes, 1
   IF (Yp(t) .GT. 0.0) THEN
    Ypd (t) = 1.0d0
   ELSE
    Ypd (t) = 0.0d0
   END IF
  END DO
 
  B = 0.0d0
  i = 0
  it = 0
 
  DO while (f /=  1)
!     print *, "Iteration ", it
   p = 1.0d0 / (1.0d0+Exp(-matmul(X, B)))
   IF (any(p > 0.97)) THEN
!    PRINT *, "WARNING: logistic regression diverging"
    f = 1
   ELSE
 
    YN = Ypd - p
    v = 0.0d0
    DO t = 1, ntimes, 1
     v (t, t) = p (t) * (1.0d0-p(t))
    END DO
    XV = matmul (v, X)
    CALL least_squares (XV, YN, TX, BN)
 
    f = 1
    DO i = 1, nvars + 1, 1
     IF (BN(i) .GT. 1.0E-04 .OR. BN(i) .LT.-1.0E-04) THEN
      f = 0
     END IF
    END DO
    IF (it > 8) THEN
!     PRINT *, "WARNING: logistic regression failed to converge"
     f = 1
    END IF
 
    B = B + BN
!        print *, "Bnew = ", B
    DEALLOCATE (BN)
 
   END IF
   it = it + 1
  END DO
 
END SUBROUTINE logistic_regression
 
 
SUBROUTINE logistic_regressionrf (X, Y, TX, B)
  USE type
  IMPLICIT NONE
 
  INTERFACE
   SUBROUTINE least_squares (X, Y, TX, B)
    USE type
    REAL (DP), INTENT (IN) :: X (:, :)
    REAL (DP), INTENT (IN) :: Y (:)
    REAL (DP), INTENT (IN) :: TX (:, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
   END SUBROUTINE least_squares
  END INTERFACE
 
  REAL (DP), INTENT (IN) :: X (:, :)
  REAL (DP), INTENT (IN) :: Y (:)
  REAL (DP), INTENT (IN) :: TX (:, :)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: B (:)
 
  REAL (DP), ALLOCATABLE :: Ypd (:), p (:), YN (:), BN (:), v (:, :), XV (:, :)
  INTEGER (I4B) :: nvars, ntimes, i, t, f, it
  REAL (DP) :: D
 
  nvars = size (X, 2) - 1
  ntimes = size (Y)
 
  ALLOCATE (B(nvars+1))
  ALLOCATE (Ypd(ntimes))
  ALLOCATE (YN(ntimes))
  ALLOCATE (p(ntimes))
  ALLOCATE (v(ntimes, ntimes))
  ALLOCATE (XV(ntimes, nvars+1))
 
  DO t = 1, ntimes, 1
   IF (Y(t) .NE. 0.0) THEN
    Ypd (t) = 1.0d0
   ELSE
    Ypd (t) = 0.0d0
   END IF
  END DO
 
  B = 0.0d0
  i = 0
  it = 0
  !print *, "B = ", B
  DO while (f /=  1)
!     print *, "Iteration ", it
   p = 1.0d0 / (1.0d0+Exp(-matmul(X, B)))
!     print *, "Pie = ", P
   IF (any(p > 0.97)) THEN
!    PRINT *, "WARNING: logistic regression diverging"
    f = 1
   ELSE
 
    YN = Ypd - p
    v = 0.0d0
    DO t = 1, ntimes, 1
     v (t, t) = p (t) * (1.0d0-p(t))
    END DO
    XV = matmul (v, X)
 
    CALL least_squares (XV, YN, TX, BN)
 
    f = 1
    DO i = 1, nvars + 1, 1
     IF (BN(i) .GT. 1.0E-04 .OR. BN(i) .LT.-1.0E-04) THEN
      f = 0
     END IF
    END DO
    IF (it > 8) THEN
!     PRINT *, "WARNING: logistic regression failed to converge"
     f = 1
    END IF
 
    B = B + BN
    DEALLOCATE (BN)
 
   END IF
   it = it + 1
  END DO
 
 
END SUBROUTINE logistic_regressionrf
 
SUBROUTINE generic_corr (stn_data, tair_data, lag, window, auto_corr, t_p_corr)
  USE type
 
  IMPLICIT NONE
 
!input
  REAL (DP), INTENT (IN) :: stn_data (:)
  REAL (DP), INTENT (IN) :: tair_data (:, :)
  INTEGER (I4B), INTENT (IN) :: lag
  INTEGER (I4B), INTENT (IN) :: window
 
 
!output
  REAL (DP), INTENT (OUT) :: auto_corr
  REAL (DP), INTENT (OUT) :: t_p_corr
 
!local variables
  REAL (DP), ALLOCATABLE :: tmean (:), trange (:)
  REAL (DP), ALLOCATABLE :: moving_avg (:, :)
  REAL (DP) :: lag_0_mean
  REAL (DP) :: lag_n_mean
  REAL (DP) :: lag_0_var
  REAL (DP) :: lag_n_var
  REAL (DP) :: lag_0_sum
  REAL (DP) :: lag_n_sum
  REAL (DP) :: cov
  REAL (DP) :: lag_0_pmean, lag_0_pvar, lag_0_psum
  REAL (DP) :: trange_mean, trange_sum, trange_var
 
  REAL (DP) :: tmp_tmean, tmp_trange
 
  INTEGER (I4B) :: i, j, tmp_cnt
  INTEGER (I4B) :: ntimes
  INTEGER (I4B) :: half_window
  INTEGER (I4B) :: cnt_sums
  INTEGER (I4B) :: data_cnt
 
!code
  ntimes = size (stn_data)
 
  ALLOCATE (tmean(ntimes))
  ALLOCATE (trange(ntimes))
  ALLOCATE (moving_avg(2, ntimes))
 
  data_cnt = 0
  DO i = 1, ntimes, 1
   IF (tair_data(1, i) .GT.-100.0 .AND. tair_data(2, i) .GT.-100.0) THEN
    tmean (i) = ((tair_data(1, i)+tair_data(2, i))/2.d0) + 273.15d0
    trange (i) = (tair_data(2, i)-tair_data(1, i)/2.d0)
    data_cnt = data_cnt + 1
   ELSE
    tmean (i) = - 999.0d0
    trange (i) = - 999.0d0
   END IF
  END DO
 
 
  half_window = floor (window/2.0d0)
 
!do the lag correlation for temperature
!first compute the moving average for climo removal
!need to check for missing values....
 
  DO i = 1, ntimes, 1
   IF (i .LT. half_window) THEN
    tmp_tmean = 0.0
    tmp_trange = 0.0
    tmp_cnt = 0
    DO j = 1, window, 1
     IF (tmean(j) .GT.-100.0) THEN
      tmp_tmean = tmp_tmean + tmean (j)
      tmp_trange = tmp_trange + trange (j)
      tmp_cnt = tmp_cnt + 1
     END IF
     IF (tmp_cnt .GT. 0) THEN
      moving_avg (1, i) = tmp_tmean / real (tmp_cnt, kind(DP))
      moving_avg (2, i) = tmp_trange / real (tmp_cnt, kind(DP))
     ELSE
      moving_avg (1, i) = - 999.0
      moving_avg (2, i) = - 999.0
     END IF
    END DO
 
   ELSE IF (i .GT. ntimes-half_window) THEN
    tmp_tmean = 0.0
    tmp_trange = 0.0
    tmp_cnt = 0
    DO j = 1, window, 1
     IF (tmean(ntimes-j) .GT.-100.0) THEN
      tmp_tmean = tmp_tmean + tmean (ntimes-j)
      tmp_trange = tmp_trange + trange (ntimes-j)
      tmp_cnt = tmp_cnt + 1
     END IF
    END DO
    IF (tmp_cnt .GT. 0) THEN
     moving_avg (1, i) = tmp_tmean / real (tmp_cnt, kind(DP))
     moving_avg (2, i) = tmp_trange / real (tmp_cnt, kind(DP))
    ELSE
     moving_avg (1, i) = - 999.0
     moving_avg (2, i) = - 999.0
    END IF
   ELSE
    tmp_tmean = 0.0
    tmp_trange = 0.0
    tmp_cnt = 0
    DO j = - half_window, half_window, 1
     IF (tmean(i+j) .GT.-100.0) THEN
      tmp_tmean = tmp_tmean + tmean (i+j)
      tmp_trange = tmp_trange + trange (i+j)
      tmp_cnt = tmp_cnt + 1
     END IF
    END DO
    IF (tmp_cnt .GT. 0) THEN
     moving_avg (1, i) = tmp_tmean / real (tmp_cnt, kind(DP))
     moving_avg (2, i) = tmp_trange / real (tmp_cnt, kind(DP))
    ELSE
     moving_avg (1, i) = - 999.0
     moving_avg (2, i) = - 999.0
    END IF
   END IF
  END DO
 
 
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
 
  DO i = lag + 1, ntimes, 1
   IF (tmean(i) .GT.-100.0 .AND. tmean(i-lag) .GT.-100.0 .AND. moving_avg(1, i) .GT.-100.0 .AND. moving_avg(1, i-lag) &
  & .GT.-100.0) THEN
    lag_n_sum = lag_n_sum + tmean (i-lag) - moving_avg (1, i-lag)
    lag_0_sum = lag_0_sum + tmean (i) - moving_avg (1, i)
    cnt_sums = cnt_sums + 1
   END IF
  END DO
 
  lag_0_mean = lag_0_sum / real (cnt_sums, kind(DP))
  lag_n_mean = lag_n_sum / real (cnt_sums, kind(DP))
 
!compute variance,covariance
  DO i = lag + 1, ntimes, 1
   IF (tmean(i) .GT.-100.0 .AND. tmean(i-lag) .GT.-100.0 .AND. moving_avg(1, i) .GT.-100.0 .AND. moving_avg(1, i-lag) &
  & .GT.-100.0) THEN
    lag_n_var = lag_n_var + ((tmean(i-lag)-moving_avg(1, i-lag))-lag_n_mean) ** 2
    lag_0_var = lag_0_var + ((tmean(i)-moving_avg(1, i))-lag_0_mean) ** 2
    cov = cov + ((tmean(i-lag)-moving_avg(1, i-lag))-lag_n_mean) * ((tmean(i)-moving_avg(1, i))-lag_0_mean)
   END IF
  END DO
 
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
  DO i = 1, ntimes, 1
   IF (trange(i) .GT.-100 .AND. stn_data(i) .GT.-100.0) THEN
      !compute for precip mean
    lag_0_psum = lag_0_psum + stn_data (i)
      !anomaly means of trange
    trange_sum = trange_sum + (trange(i)-moving_avg(2, i))
    tmp_cnt = tmp_cnt + 1
   END IF
  END DO
  lag_0_pmean = lag_0_psum / real (tmp_cnt, kind(DP))
  trange_mean = trange_sum / real (tmp_cnt, kind(DP))
 
!compute variance and covariance
  DO i = 1, ntimes, 1
   IF (trange(i) .GT.-100 .AND. stn_data(i) .GT.-100.0) THEN
    lag_0_pvar = lag_0_pvar + (stn_data(i)-lag_0_pmean) ** 2
    trange_var = trange_var + ((trange(i)-moving_avg(2, i))-trange_mean) ** 2
    cov = cov + ((trange(i)-moving_avg(2, i))-trange_mean) * (stn_data(i)-lag_0_pmean)
   END IF
  END DO
 
  t_p_corr = cov / (Sqrt(lag_0_pvar)*Sqrt(trange_var))
  IF (Sqrt(lag_0_pvar)*Sqrt(trange_var) .LE. 0.00001) THEN
   t_p_corr = 0.0
  END IF
!in some situations, there are very limited data used for calculations
!set those cases to missing value
  IF (data_cnt .LT. real(ntimes)*0.25) THEN
   auto_corr = - 999.0
   t_p_corr = - 999.0
  END IF
 
END SUBROUTINE generic_corr
