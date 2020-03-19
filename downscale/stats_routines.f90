! This source code file has various statistical functions, eg, normalizing, weighting, correlation
! They were originally in coefficients.f90 before it was split up

subroutine normalize_x (x)
  use type
  implicit none
 
  real (dp), intent (inout) :: x (:, :)
 
  real (dp) :: mean, stdev, sum_x, sum_x2
  integer (i4b) :: v, t
  integer (i4b) :: nvars, ntimes
 
  nvars = size (x, 2) - 1
  ntimes = size (x, 1)
 
  do v = 2, nvars + 1, 1
    sum_x = 0.0d0
    sum_x2 = 0.0d0
    do t = 1, ntimes, 1
      sum_x = sum_x + x (t, v)
      sum_x2 = sum_x2 + x (t, v) ** 2
    end do
    mean = sum_x / real (ntimes)
    stdev = sqrt ((real(ntimes)*sum_x2-sum_x**2)/(real(ntimes)*real(ntimes-1)))
    do t = 1, ntimes, 1
      if (stdev .eq. 0.0) then
        x (t, v) = x (t, v)
      else
        x (t, v) = (x(t, v)-mean) / stdev
      end if
    end do
  end do
 
end subroutine normalize_x
 
!added AJN Sept 2013
subroutine normalize_xv (x, weight, yp, smax)
  use type
  implicit none
 
  real (dp), intent (in) :: x (:)
  real (dp), intent (in) :: weight (:)
  integer (i4b), intent (in) :: yp (:)
  real (dp), intent (out) :: smax
 
!  real(DP) :: mean, stdev, sum_x2, sum_x
  real (dp) :: sum_x2, sum_x
  real (dp) :: sum_stdev, sum_std, smin
 
  real (dp) :: mean_all, sum_weight, sum_xw
 
  integer (i4b) :: t
  integer (i4b) :: ntimes
 
  ntimes = size (x)
 
  smin = 9999.0
  smax = 0.0
 
  sum_xw = 0.0d0
  sum_weight = 0.0d0
 
  sum_x = 0.0d0
  sum_x2 = 0.0d0
  do t = 1, ntimes, 1
    if (x(t) >= -100) then
      sum_x = sum_x + x (t)
      sum_x2 = sum_x2 + x (t) ** 2
      sum_xw = sum_xw + weight (t) * x (t)
      sum_weight = sum_weight + weight (t)
 
      if (x(t) .le. smin) then
        smin = x (t)
      end if
      if (x(t) .ge. smax) then
        smax = x (t)
      end if
 
    end if
  end do
  if(smin > 1000.) then
    smin = 0.0
  end if 

  return
 
end subroutine normalize_xv
 
 
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
subroutine normalize_y (texp, y)
  use type
  implicit none
 
  real (dp), intent (in) :: texp !transform exponent
  real (dp), intent (inout) :: y (:)
  integer (i4b) :: t
  integer (i4b) :: ntimes
 
  ntimes = size (y)
 
  do t = 1, ntimes, 1

    if(Y(t) > 0) then
      !power law
      !!y (t) = y (t) ** (1.0d0/texp)

      !box-cox
      y(t) = ((Y(t)**(1.0/texp))-1.0)/(1.0/texp)
    else
!      y (t) = -999.0
      y (t) = -3.0
    end if 
  end do
 
end subroutine normalize_y
 
! SUBROUTINE calc_weights()
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
subroutine calc_weights (tt, x, w)
  use type
  implicit none
 
  integer (i4b), intent (in) :: tt
  real (dp), intent (in) :: x (:, :)
  real (dp), allocatable, intent (out) :: w (:, :)
  real (dp) :: sum
  integer (i4b) :: v, t
  integer (i4b) :: nvars, ntimes
 
  nvars = size (x, 2) - 1
  ntimes = size (x, 1)
 
  allocate (w(ntimes, ntimes))
 
!  do v = 2, nvars+1, 1
!     print *, "X(tt,v) : ",  X(tt,v)
!  enddo
  do t = 1, ntimes, 1
    w (t, :) = 0.0d0
    if (t /= tt) then
      sum = 0.0d0
      do v = 2, nvars + 1, 1
           !print *, "X(t,v) : ",  X(t,v)
        sum = sum + (x(t, v)-x(tt, v)) ** 2
      end do
      w (t, :) = 0.0d0
        !W(t,t) = 1.0
 
!old version AJN
!        W(t,t) = sqrt(sum / nvars)
 
      w (t, t) = 1 / sqrt (sum/nvars)
    end if
!     print *, "W : ",  t, W(t,t)
  end do
 
end subroutine calc_weights
 
! Great circle distance calculation
! Output in nm
subroutine calc_distance_weight (maxd, lat1, lon1, lat2, lon2, weight)
  use type
  implicit none
 
  real (dp), intent (in) :: maxd, lat1, lon1, lat2, lon2
  real (dp), intent (out) :: weight
 
  real (dp) :: dist, lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927
 
  lat1r = lat1 * pi / 180.0d0
  lon1r = lon1 * pi / 180.0d0
  lat2r = lat2 * pi / 180.0d0
  lon2r = lon2 * pi / 180.0d0
  dist = ((180*60)/pi) * &
 & (2*asin(sqrt((sin((lat1r-lat2r)/2))**2+cos(lat1r)*cos(lat2r)*(sin((lon1r-lon2r)/2))**2)))
  if (dist .gt. maxd) then
    weight = 0.0d0
  else
    weight = (1.0d0-(dist/maxd)**3) ** 3    ! inverse cubic for weight fcn
    ! weight = 1.0d0 - (dist/maxd)**0.5
  end if
 
end subroutine calc_distance_weight
 
! Great circle distance calculation
! Output in nm
subroutine calc_distance (lat1, lon1, lat2, lon2, dist)
  use type
  implicit none
 
  real (dp), intent (in) :: lat1, lon1, lat2, lon2
  real (dp), intent (out) :: dist
 
  real (dp) :: lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927
 
  lat1r = lat1 * pi / 180.0d0
  lon1r = lon1 * pi / 180.0d0
  lat2r = lat2 * pi / 180.0d0
  lon2r = lon2 * pi / 180.0d0
  dist = ((180*60)/pi) * &
 & (2*asin(sqrt((sin((lat1r-lat2r)/2))**2+cos(lat1r)*cos(lat2r)*(sin((lon1r-lon2r)/2))**2)))
 
end subroutine calc_distance
 
! added AJN Sept 2013
! modified AWW Feb 2016, change stn_data name to prcp_data
subroutine generic_corr (prcp_data, tair_data, lag, window, auto_corr, t_p_corr)
  use type
 
  implicit none
 
  ! input
  real (dp), intent (in) :: prcp_data (:)
  real (dp), intent (in) :: tair_data (:, :)
  integer (i4b), intent (in) :: lag
  integer (i4b), intent (in) :: window
 
  ! output
  real (dp), intent (out) :: auto_corr
  real (dp), intent (out) :: t_p_corr
 
  ! local variables
  real (dp), allocatable :: tmean (:), trange (:)
  real (dp), allocatable :: moving_avg (:, :)
  real (dp) :: lag_0_mean
  real (dp) :: lag_n_mean
  real (dp) :: lag_0_var
  real (dp) :: lag_n_var
  real (dp) :: lag_0_sum
  real (dp) :: lag_n_sum
  real (dp) :: cov
  real (dp) :: lag_0_pmean, lag_0_pvar, lag_0_psum
  real (dp) :: trange_mean, trange_sum, trange_var
 
  real (dp) :: tmp_tmean, tmp_trange
 
  integer (i4b) :: i, j, tmp_cnt
  integer (i4b) :: ntimes
  integer (i4b) :: half_window
  integer (i4b) :: local_window
  integer (i4b) :: cnt_sums
  integer (i4b) :: data_cnt
 
  ! code
  ntimes = size (prcp_data)
  allocate (tmean(ntimes))
  allocate (trange(ntimes))
  allocate (moving_avg(2, ntimes))

  if(window .gt. ntimes) then
   local_window = ntimes
  else
    local_window = window
  endif
 
  data_cnt = 0
  do i = 1, ntimes, 1
    ! here assume data are in C, switch to Kelvin?  Is that needed?
    ! Calculate Trange
    if (tair_data(1, i) .gt.-100.0 .and. tair_data(2, i) .gt.-100.0) then
      tmean (i) = ((tair_data(1, i)+tair_data(2, i))/2.d0) + 273.15d0
      trange (i) = (tair_data(2, i)-tair_data(1, i)/2.d0)
      data_cnt = data_cnt + 1
    else
      tmean (i) = - 999.0d0
      trange (i) = - 999.0d0
    end if
  end do
 
  half_window = floor (local_window/2.0d0)
 
  ! do the lag correlation for temperature
  ! first compute the moving average for climo removal in window
  ! need to check for missing values....
  !    handle edge cases for window averaging near bounds of 1 or ntimes
  do i = 1, ntimes, 1
      tmp_tmean = 0.0
      tmp_trange = 0.0
      tmp_cnt = 0

    if (i .le. half_window) then
      do j = 1, local_window, 1
        if (tmean(j) .gt.-100.0) then
          tmp_tmean = tmp_tmean + tmean (j)
          tmp_trange = tmp_trange + trange (j)
          tmp_cnt = tmp_cnt + 1
        end if
      end do
    else if (i .gt. ntimes-half_window) then
      do j = 1, local_window, 1
        if (tmean(ntimes-j+1) .gt.-100.0) then
          tmp_tmean = tmp_tmean + tmean (ntimes-j+1)
          tmp_trange = tmp_trange + trange (ntimes-j+1)
          tmp_cnt = tmp_cnt + 1
        end if
      end do
    else
      do j = - half_window, half_window, 1
        if (tmean(i+j) .gt.-100.0) then
          tmp_tmean = tmp_tmean + tmean (i+j)
          tmp_trange = tmp_trange + trange (i+j)
          tmp_cnt = tmp_cnt + 1
        end if
      end do
    end if

    if (tmp_cnt .gt. 0) then
      moving_avg (1, i) = tmp_tmean / real (tmp_cnt, kind(dp))
      moving_avg (2, i) = tmp_trange / real (tmp_cnt, kind(dp))
    else
      moving_avg (1, i) = -999.0
      moving_avg (2, i) = -999.0
    end if

  end do
 
  ! print *,'tmean ',moving_avg(1,1:2)
  ! print *,'trange ',moving_avg(2,1:2)
 
  ! only use portions of timeseries to compute auto_corr
  ! need to go through and check to see if values exist for lag-0 and lag-n and moving_avg
  ! if values do not exist for any of the three, don't add to running sums
 
  ! compute means
  lag_0_sum = 0.0d0
  lag_n_sum = 0.0d0
  lag_0_var = 0.0d0
  lag_n_var = 0.0d0
  cov = 0.0d0
  cnt_sums = 0
 
  do i = lag + 1, ntimes, 1
    if (tmean(i) .gt.-100.0 .and. tmean(i-lag) .gt.-100.0 .and. moving_avg(1, i) .gt.-100.0 .and. &
   & moving_avg(1, i-lag) .gt.-100.0) then
      lag_n_sum = lag_n_sum + tmean (i-lag) - moving_avg (1, i-lag)
      lag_0_sum = lag_0_sum + tmean (i) - moving_avg (1, i)
      cnt_sums = cnt_sums + 1
    end if
  end do
 
  if(cnt_sums .gt. 0)then
    lag_0_mean = lag_0_sum/real(cnt_sums,kind(dp))
    lag_n_mean = lag_n_sum/real(cnt_sums,kind(dp))
  else
    lag_0_mean = 0.0
    lag_n_mean = 0.0
  endif
 
  ! compute variance, covariance
  do i = lag + 1, ntimes, 1
    if (tmean(i) .gt.-100.0 .and. tmean(i-lag) .gt.-100.0 .and. moving_avg(1, i) .gt.-100.0 .and. &
   & moving_avg(1, i-lag) .gt.-100.0) then
      lag_n_var = lag_n_var + ((tmean(i-lag)-moving_avg(1, i-lag))-lag_n_mean) ** 2
      lag_0_var = lag_0_var + ((tmean(i)-moving_avg(1, i))-lag_0_mean) ** 2
      cov = cov + ((tmean(i-lag)-moving_avg(1, i-lag))-lag_n_mean) * ((tmean(i)-moving_avg(1, &
     & i))-lag_0_mean)
    end if
  end do
 
  ! compute autocorrelation
  if(lag_0_var .gt. 0 .and. lag_n_var .gt. 0) then
    auto_corr = cov/(sqrt(lag_0_var)*sqrt(lag_n_var))
  else
    auto_corr = 0.0
    print*, "INFO: lag_0_var or lag_n_var < 0; setting auto_corr to zero"
  endif
 
  !------------------------------------------------------
  ! now do the TxP cross-correlation for a single station
  ! do correlation on trange, not tmean
  !------------------------------------------------------
 
  lag_0_pmean = 0.0d0
  lag_0_pvar = 0.0d0
  lag_0_psum = 0.0d0
  trange_sum = 0.0d0
  trange_mean = 0.0d0
  trange_var = 0.0d0
  cov = 0.0d0
  tmp_cnt = 0
 
  ! again need to check for missing values....
  do i = 1, ntimes, 1
    ! if (trange(i) .gt.-100 .and. prcp_data(i) .gt.-100.0) then
    if (trange(i) .gt.-100 .and. prcp_data(i) .gt.-1.0) then   ! AWW-feb2016
      !compute for precip mean
      lag_0_psum = lag_0_psum + prcp_data (i)
      !anomaly means of trange
      trange_sum = trange_sum + (trange(i)-moving_avg(2, i))
      tmp_cnt = tmp_cnt + 1
    end if
  end do
  if(tmp_cnt .gt. 0) then
    lag_0_pmean = lag_0_psum/real(tmp_cnt,kind(dp))
    trange_mean = trange_sum/real(tmp_cnt,kind(dp))
  else
    lag_0_pmean = 0.0
    trange_mean = 0.0
  end if
 
  ! compute variance and covariance
  do i = 1, ntimes, 1
    if (trange(i) .gt.-100 .and. prcp_data(i) .gt.-100.0) then
      lag_0_pvar = lag_0_pvar + (prcp_data(i)-lag_0_pmean) ** 2
      trange_var = trange_var + ((trange(i)-moving_avg(2, i))-trange_mean) ** 2
      cov = cov + ((trange(i)-moving_avg(2, i))-trange_mean) * (prcp_data(i)-lag_0_pmean)
    end if
  end do
  ! print *,cov,lag_0_pvar,trange_var,trange_mean,lag_0_pmean,lag_0_psum
 
  if(lag_0_pvar .le. 0.00001 .or. trange_var .le. 0.00001) then
    t_p_corr = 0.0
  else
    t_p_corr = cov/(sqrt(lag_0_pvar)*sqrt(trange_var))
  endif

  ! print *,t_p_corr
 
  ! in some situations, there are very limited data used for calculations
  ! set those cases to missing value
  if (data_cnt .lt. real(ntimes)*0.25) then
    auto_corr = - 999.0
    t_p_corr = - 999.0
  end if
 
  ! print *,cov,lag_0_var,lag_n_var,auto_corr
 
end subroutine generic_corr
