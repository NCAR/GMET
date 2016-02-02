subroutine estimate_coefficients(D, nvars, Lats, Lons, Times, stnid, stnlat, stnlon, &
  stnalt, stnvar, site_var, site_list, C, POC, error)
  use type
  implicit none

  interface

     subroutine read_station(stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss,&
                             vals_miss_t, error)
       use type
       character(len=100), intent(in) :: stnvar
       character(len=100), intent(in) :: stnid
       character(len=100), intent(in) :: site_var,site_var_t
       character(len=500), intent(in) :: site_list
       real(DP), intent(in) :: Times(:)
       real(DP), allocatable, intent(out) :: vals(:),tair_vals(:,:)
       logical, allocatable, intent(out) :: vals_miss(:),vals_miss_t(:)
       integer, intent(out) :: error
     end subroutine read_station

     subroutine normalize_X(X)
       use type
       real(DP), intent(inout) :: X(:,:)
     end subroutine normalize_X

     subroutine normalize_Y(texp,Y)
       use type
       real(DP), intent(in)    :: texp  !transform exponent
       real(DP), intent(inout) :: Y(:)
     end subroutine normalize_Y

     subroutine calc_weights(Times, tt, X, W)
       use type
       real(DP), intent(in) :: Times(:)
       integer(I4B), intent(in) :: tt
       real(DP), intent(in) :: X(:,:)
       real(DP), allocatable, intent(out) :: W(:,:)
     end subroutine calc_weights
     
     subroutine least_squares(X, Y, TX, B)
       use type
       real(DP), intent(in) :: X(:,:)
       real(DP), intent(in) :: Y(:)
       real(DP), intent(in) :: TX(:,:)
       real(DP), allocatable, intent(out) :: B(:)
     end subroutine least_squares

     subroutine logistic_regressionrf(X, Y, TX, B)
       use type
       real(DP), intent(in) :: X(:,:)
       real(DP), intent(in) :: Y(:)
       real(DP), intent(in) :: TX(:,:)
       real(DP), allocatable, intent(out) :: B(:)
     end subroutine logistic_regressionrf

  end interface

  real(DP), intent(in) :: D(:,:,:), Lats(:), Lons(:)
  real(DP), intent(in) :: Times(:)
  integer(I4B), intent(in) :: nvars
  character (len = 100), intent(in) :: stnid(:)
  real(DP), intent(in) :: stnlat(:), stnlon(:), stnalt(:)
  character(len=100), intent(in) :: stnvar
  character(len=100), intent(in) :: site_var
  character(len=500), intent(in) :: site_list
  real(DP), allocatable, intent(out) :: C(:,:,:), POC(:,:,:)
  integer, intent(out) :: error

  real(DP), allocatable :: X(:,:), Y(:), XS(:,:), XP(:,:), TWXS(:,:), YS(:), W(:,:), B(:),tair_vals(:,:)

  real(DP), allocatable :: TS(:)
  logical, allocatable :: Y_miss(:),Y_miss_t(:) 
  integer(I4B) :: ntimes, nlats, nlons, nstns
  integer(I4B) :: i, j, k, t, v, tt
  real(DP) :: minlondis, minlatdis
  integer(I4B) :: minlatk, minlonj, gridindex
  integer(I4B) :: vshape(2)
  character(len=100)     :: site_var_t 
  real(DP)	:: transform_exp  !precipitation transform variable, the exponent of the transform  norm_pcp = pcp^(1/transform_exp)


  error = 0
  transform_exp = 4.0d0

  ntimes = size(Times)
  nlats = size(Lats)
  nlons = size(Lons)
  nstns = size(stnlat)
  print *, ntimes, nlats, nlons, nstns

  allocate(X(ntimes, nvars+1))
  if(trim(stnvar) .EQ. "PRCP") then
    allocate(POC(nstns, ntimes, nvars+1))
    POC = 0.0d0
  end if
  allocate(C(nstns, ntimes, nvars+1)) 
  C = 0.0d0

!  do i = 5, 5, 1!nstns, 1
  do i = 1,nstns,1

     minlondis = 360.0
     minlonj = -1
     do j = 1, nlons, 1
        if( ABS(stnlon(i) - Lons(j)) < minlondis ) then
           minlondis = ABS(stnlon(i) - Lons(j))
           minlonj = j
        endif
     enddo
     minlatdis = 180.0
     minlatk = -1
     do k = 1, nlats, 1
        if( ABS(stnlat(i) - Lats(k)) < minlatdis ) then
           minlatdis = ABS(stnlat(i) - Lats(k))
           minlatk = k
        endif
     enddo

     if(minlonj == -1 .OR. minlatk == -1) then
        print *, "Failed to find closest grid point for station: ", trim(stnid(i))
        error = 1
        return
     endif

     gridindex = ((minlatk-1) * nlons) + minlonj

     print *, "Station: ", trim(stnid(i)), stnlat(i), stnlon(i)
     print *, "Closest Grid point: ", gridindex, Lats(minlatk), Lons(minlonj)

     X(:,1) = 1.0
     do v = 1, nvars, 1
        X(:,v+1) = D(v, gridindex, :)
        do t = 1, ntimes, 1
           if(X(t,v+1) /= 0.0) exit
           if(t == ntimes) then
              print *, "ERROR: var ", v, " is all zero for station ", i
              X(:,v+1) = 1.0
           endif
        enddo
     enddo

     call read_station(stnvar, stnid(i), site_var, site_var_t, site_list, Times, Y, tair_vals,&
                       Y_miss, Y_miss_t, error)

     if(count(Y_miss) < ntimes) then
        allocate(TS(count(Y_miss)))
        allocate(XS(count(Y_miss), nvars+1))
	allocate(XP(count(Y_miss), nvars+1))
        allocate(TWXS(nvars+1, count(Y_miss)))
        allocate(YS(count(Y_miss)))
        TS(:) = pack(Times, Y_miss)
        do v = 1, nvars+1, 1
           XS(:,v) = pack(X(:,v), Y_miss)
        enddo
	XP(:,:) = XS(:,:)
        YS(:) = pack(Y, Y_miss)
     else
        allocate(TS(ntimes))
        allocate(XS(ntimes, nvars+1))
	allocate(XP(ntimes, nvars+1))
        allocate(TWXS(nvars+1, ntimes))
        allocate(YS(ntimes))
        TS(:) = Times(:)
        XS(:,:) = X(:,:)
	XP(:,:) = XS(:,:)
        YS(:) = Y(:)
     endif

     call normalize_X(XS)
     call normalize_Y(transform_exp,YS)
!     print *, "XS:", XS
 !    print *, "YS:", YS
     tt = 1
     do t = 1, ntimes, 1

        if(Y_miss(t) .EQV. .TRUE.) then
!           print *, "tt = ", tt
           call calc_weights(TS, tt, XS, W)
           TWXS = matmul(transpose(XS), W)
           if(trim(stnvar) .EQ. "PRCP") then

              call logistic_regressionrf(XP, YS, TWXS, B)

              POC(i,t,:) = B(:)
              deallocate(B)
           endif

           call least_squares(XP, YS, TWXS, B)
           C(i,t,:) = B(:)

           deallocate(B)
!           print *, "Coefficients: ", C(i,t,:)

           tt = tt +1
        else
           if(trim(stnvar) .EQ. "PRCP") then
              POC(i,t,:) = -999.99
           endif
           C(i,t,:) = -999.99
        endif
     end do

     deallocate(YS)
     deallocate(TWXS)
     deallocate(XS)
     deallocate(TS)
     deallocate(XP)

  enddo

end subroutine estimate_coefficients


subroutine estimate_precip(X, Z, nsta, ngrid, maxDistance, Times,  &
     stnid, stnvar, site_var, site_var_t, site_list, PCP, POP, PCPERR, &
     tmean,tmean_err,trange,trange_err,mean_autocorr, mean_tp_corr,    &
     Y_mean,Y_std, y_std_all,y_min, y_max, error, &
    PCP_2, POP_2, PCPERR_2,tmean_2,tmean_err_2,trange_2,trange_err_2)

  use strings
  use utim
  use type
  implicit none

  interface
     subroutine read_transform_exp(ntimes,file_name,texp)
       use type

       integer(I4B), intent(in)		:: ntimes
       character (len = *), intent(in)	:: file_name
       real(DP),allocatable,intent(out)	:: texp(:)

     end subroutine read_transform_exp

     subroutine read_station(stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss,&
                             vals_miss_t, error)
       use type
       character(len=100), intent(in) :: stnvar
       character(len=100), intent(in) :: stnid
       character(len=100), intent(in) :: site_var
       character(len=100), intent(in) :: site_var_t
       character(len=500), intent(in) :: site_list
       real(DP), intent(in) :: Times(:)
       real(DP), allocatable, intent(out) :: vals(:),tair_vals(:,:)
       logical, allocatable, intent(out) :: vals_miss(:),vals_miss_t(:)
       integer, intent(out) :: error
     end subroutine read_station

     subroutine normalize_X(X)
       use type
       real(DP), intent(inout) :: X(:,:)
     end subroutine normalize_X

     subroutine normalize_Xv(X,weight,mean,stdev,stdev_all,smin,smax,Yp)
       use type
       real(DP), intent(inout) :: X(:)
       real(DP), intent(in)  :: weight(:)
       real(DP), intent(out) :: mean
       real(DP), intent(out) :: stdev
       real(DP), intent(out) :: stdev_all
       real(DP), intent(out) :: smin
       real(DP), intent(out) :: smax
       integer(I4B),intent(out) :: Yp(:)
     end subroutine normalize_Xv

     subroutine normalize_Y(texp,Y)
       use type
       real(DP), intent(in)    :: texp  !transform exponent
       real(DP), intent(inout) :: Y(:)
     end subroutine normalize_Y

     subroutine calc_weights(Times, tt, X, W)
       use type
       real(DP), intent(in) :: Times(:)
       integer(I4B), intent(in) :: tt
       real(DP), intent(in) :: X(:,:)
       real(DP), allocatable, intent(out) :: W(:,:)
     end subroutine calc_weights
     
     subroutine least_squares(X, Y, TX, B)
       use type
       real(DP), intent(in) :: X(:,:)
       real(DP), intent(in) :: Y(:)
       real(DP), intent(in) :: TX(:,:)
       real(DP), allocatable, intent(out) :: B(:)
     end subroutine least_squares

     subroutine logistic_regression(X, Y, TX, Yp, B)
       use type
       real(DP), intent(in) :: X(:,:)
       real(DP), intent(in) :: Y(:)
       real(DP), intent(in) :: TX(:,:)
       integer(I4B),intent(in) :: Yp(:)
       real(DP), allocatable, intent(out) :: B(:)
     end subroutine logistic_regression

    subroutine generic_corr(stn_data,tair_data,lag,window,auto_corr,t_p_corr)
      use type

      real(DP),intent(in)       :: stn_data(:)
      real(DP),intent(in)       :: tair_data(:,:)
      integer(I4B),intent(in)   :: lag
      integer(I4B),intent(in)   :: window
      real(DP),intent(out)   :: auto_corr
      real(DP),intent(out)   :: t_p_corr

    end subroutine generic_corr

    subroutine calc_distance(lat1,lon1,lat2,lon2,dist)
      use type
      implicit none

      real(DP), intent(in) :: lat1, lon1, lat2, lon2
      real(DP), intent(out) :: dist

    end subroutine calc_distance


    subroutine heapsort(n,ra,rn)
      use type
      implicit none

      integer(I4B),intent(in)		:: N
      integer(I4B),dimension(:),intent(inout)	:: RN
      real(DP),dimension(:),intent(inout)	:: RA

    end subroutine heapsort

  end interface

  real(DP), intent(in) :: X(:,:), Z(:,:)   !station and grid point description arrays
  real(DP), intent(in) :: maxDistance		!max distance for weight function
  integer(I4B), intent(in) :: nsta, ngrid	!nuber of input stations and grid points
  real(DP), intent(in) :: Times(:)		!time step array
  character (len = 100), intent(in) :: stnid(:)	!station id array
  character(len=100), intent(in) :: stnvar, site_var, site_var_t	!control file variables
  character(len=500), intent(in) :: site_list				!file name of station list
  real(SP), allocatable, intent(out) :: PCP(:,:), POP(:,:), PCPERR(:,:)	!output variables for precipitation
  real(SP), allocatable, intent(out) :: tmean(:,:),tmean_err(:,:)   !OLS tmean estimate and error
  real(SP), allocatable, intent(out) :: trange(:,:),trange_err(:,:) !OLS trange estimate and error

  real(SP), allocatable, intent(out) :: tmean_2(:,:),tmean_err_2(:,:)   !OLS tmean estimate and error
  real(SP), allocatable, intent(out) :: trange_2(:,:),trange_err_2(:,:) !OLS trange estimate and error
  real(SP), allocatable, intent(out) :: PCP_2(:,:), POP_2(:,:), PCPERR_2(:,:)


  integer, intent(out) :: error		!integer error flag
  real(DP),intent(out) :: mean_autocorr(:)  !mean autocorrelation from all stations over entire time period
  real(DP),intent(out) :: mean_tp_corr(:)  !mean correlation for mean temp and precip

   !vary at each grid point and time step
  real(DP),intent(out) :: Y_mean(:,:), Y_std(:,:)  !std and mean of time step precipitation
  real(DP),intent(out) :: Y_std_all(:,:)           !std of time step precip including stations with zero precip
  real(DP),intent(out) :: Y_min(:,:), Y_max(:,:)  !min & max  of normalized time step precipitation

  real(DP), allocatable :: Y(:), TWX(:,:), B(:), TX(:,:)

  real(DP), allocatable :: Y_red(:), TWX_red(:,:), TX_red(:,:) !reduced matricies
  real(DP), allocatable :: X_red(:,:) !reduced matricies

  real(DP), allocatable :: TWX_red_t(:,:), TX_red_t(:,:) !reduced matricies
  real(DP), allocatable :: X_red_t(:,:) !reduced matricies

  real(DP), allocatable :: TWX_red_tr(:,:), TX_red_tr(:,:) !reduced matricies
  real(DP), allocatable :: X_red_tr(:,:) !reduced matricies

  real(DP), allocatable :: w_base(:,:) !initial distance weight matrix
  real(DP), allocatable :: w_pcp_1d(:),w_temp_1d(:)
  integer(I4B),allocatable :: w_pcp_1d_loc(:),w_temp_1d_loc(:)
!  real(DP), allocatable :: w_pcp(:,:), w_temp(:,:) !distance weight matrices for a specific grid point
  real(DP), allocatable :: w_pcp_red(:,:), w_temp_red(:,:) !reduced distance weigth matricies 


  real(DP), allocatable :: Y_tmean(:),Y_trange(:)			!transformed station data arrays
  real(DP), allocatable :: Y_tmean_red(:),Y_trange_red(:)   		!transformed station data arrays
  real(DP), allocatable :: stn_vals(:), stn_data(:,:), tair_data(:,:,:), stn_tair(:,:)  !original station data arrays
  real(DP), allocatable :: auto_corr(:)  		!lag-1 autocorrelation for stations over entire time period used
  real(DP), allocatable :: t_p_corr(:)   		!correlation between temp and precip 
  integer(I4B),allocatable :: Yp(:)      		!binary for logistic regression
  integer(I4B),allocatable :: Yp_red(:)  		!reduced binary for logistic regression


  logical, allocatable :: stn_miss(:), stn_miss_t(:)	!missing value logical arrays

  real(DP)              :: m_tmean, m_trange          !used to fill missing values in stations

  real(DP) :: errsum, wgtsum, p,sta_pcp,sta_temp
  real(DP) :: auto_corr_sum,tp_corr_sum
  real(DP) :: step_mean, step_std, step_std_all, step_min, step_max  !timestep statistics

  real(DP) :: rsqr, ss_tot,ss_res,vif !r-squared and variance correction

  integer(I4B) :: xsize  !size of second dimension of input X array

  integer(I4B) :: ntimes, nstns
  integer(I4B) :: t, i, g, ndata, nodata
  integer(I4B) :: ndata_t,nodata_t
  integer(I4B) :: lag,window,tc,trc
  integer(I4B) :: auto_cnt,tp_cnt

  integer(I4B) :: stn_count


  !variables for tracking closest N stations for precipitation
  integer(I4B) :: out_loc
  integer(I4B),parameter :: sta_limit = 30
  integer(I4B),allocatable :: close_loc(:,:)
  integer(I4B),allocatable :: close_count(:)

  real(DP),allocatable :: close_weights(:,:)
  real(DP),allocatable :: close_meta(:,:,:)
  real(DP)             :: min_weight
  real(DP)             :: max_distance
  real(DP),parameter   :: search_distance=1000.0

  !variables for tracking closest N stations for temperature
  integer(I4B) :: out_loc_t
  integer(I4B),allocatable :: close_loc_t(:,:)
  integer(I4B),allocatable :: close_count_t(:)

  real(DP),allocatable	:: close_weights_t(:,:)
  real(DP),allocatable	:: close_meta_t(:,:,:)
  real(DP)     		:: min_weight_t
  real(DP)     		:: max_distance_t
  

  real(DP)		:: tmp_pcp
  real(DP)		:: tmp_weight

  integer(I4B)		:: slope_flag_pcp
  integer(I4B)		:: slope_flag_temp

!variables to check for singular matrix
  real(DP),allocatable	:: tmp(:,:)
  real(DP),allocatable	:: vv(:)

!variables for timing code
  integer(I4B) :: t1, t2,count_rate
  integer(I4B) :: tg1, tg2


!variables for keeping track of max_distance modifications
  integer(I4B),allocatable	:: expand_flag(:), expand_flag_t(:)
  real(DP),allocatable 		:: expand_dist(:), expand_dist_t(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! code starts below here
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nstns = size(stnid)
!  nstns = 100
  ntimes = size(Times)

  xsize = size(X,2)

!allocate variables
  allocate(Y(nstns))
  allocate(Y_tmean(nstns),Y_trange(nstns))
  allocate(stn_data(nstns, ntimes))
  allocate(tair_data(2,nstns,ntimes))

!original allocations
!  allocate(w_pcp(nstns,nstns))
!  allocate(w_temp(nstns,nstns))
!  allocate(TWX(4,nstns))
!  allocate(TX(4,nstns))

  allocate(w_pcp_red(sta_limit,sta_limit))
  allocate(w_temp_red(sta_limit,sta_limit))
  allocate(Y_red(sta_limit))
  allocate(Y_tmean_red(sta_limit),Y_trange_red(sta_limit))
  allocate(X_red(sta_limit,xsize))
  allocate(X_red_t(sta_limit,xsize))

  allocate(w_pcp_1d(sta_limit))
  allocate(w_temp_1d(sta_limit))
  allocate(w_pcp_1d_loc(sta_limit))
  allocate(w_temp_1d_loc(sta_limit))

  allocate(tmp(6,6))
  allocate(vv(6))

  allocate(PCP(ngrid, ntimes))
  allocate(POP(ngrid, ntimes))
  allocate(PCPERR(ngrid, ntimes))

  allocate(tmean(ngrid,ntimes))
  allocate(tmean_err(ngrid,ntimes))
  allocate(trange(ngrid,ntimes))
  allocate(trange_err(ngrid,ntimes))


  allocate(PCP_2(ngrid, ntimes))
  allocate(POP_2(ngrid, ntimes))
  allocate(PCPERR_2(ngrid, ntimes))

  allocate(tmean_2(ngrid,ntimes))
  allocate(tmean_err_2(ngrid,ntimes))
  allocate(trange_2(ngrid,ntimes))
  allocate(trange_err_2(ngrid,ntimes))


  allocate(auto_corr(nstns))
  allocate(t_p_corr(nstns))

  allocate(Yp(nstns))
  allocate(Yp_red(sta_limit))

  !station limit arrays
  allocate(close_weights(ngrid,sta_limit))
  allocate(close_loc(ngrid,sta_limit))
  allocate(close_meta(5,ngrid,sta_limit))
  allocate(close_count(ngrid))

  allocate(close_weights_t(ngrid,sta_limit))
  allocate(close_loc_t(ngrid,sta_limit))
  allocate(close_meta_t(5,ngrid,sta_limit))
  allocate(close_count_t(ngrid))


  !base weight array
  allocate(w_base(ngrid,nstns))

  !max_dist tracking variables
  allocate(expand_dist(ngrid),expand_flag(ngrid))
  allocate(expand_dist_t(ngrid),expand_flag_t(ngrid))


  PCP = 0.0d0
  POP = 0.0d0
  PCPERR = 0.0d0
  auto_corr = 0.0d0
  stn_count = 1

  tmean = 0.0d0
  trange= 0.0d0
  tmean_err = 0.0d0
  trange_err = 0.0d0

  auto_corr_sum = 0.0d0
  auto_cnt      = 0
  tp_corr_sum   = 0.0d0
  tp_cnt         = 0

  w_base = 0.0d0

  expand_dist = 0.0d0
  expand_flag = 0
  expand_dist_t = 0.0d0
  expand_flag_t = 0

!  do i = 1, 151, 1
  do i = 1, nstns, 1
!print *,'station read'
!print *,stnvar,stnid(i),site_var,site_var_t,site_list,Times
     call read_station(stnvar, stnid(i), site_var, site_var_t, site_list, Times, stn_vals, stn_tair, stn_miss,&
                       stn_miss_t, error)


     stn_count = stn_count + 1


     stn_data(i,:) = stn_vals
     tair_data(1,i,:) = stn_tair(1,:)
     tair_data(2,i,:) = stn_tair(2,:)

!print *, 'precip',stn_data(i,:),stn_miss
! print *,tair_data(1,i,:),stn_miss_t

    !call subroutine that does various  correlation calculations
    !can do autocorrelations and correlation between temperature and precipitation
    !uses an n-day moving average (window) to remove "monthly" cycle from temp
    !and computes autocorrelation on the anomalies
    
    lag = 1
    window = 31
    call generic_corr(stn_data(i,:),tair_data(:,i,:),lag,window,auto_corr(i),t_p_corr(i))
!print *,'here'
!    print *,auto_corr(i)

    !compute mean autocorrelation for all stations and all times
    !check for values outside of -1 to 1
    !stations with incomplete data are set to -999
    if(auto_corr(i) .ge. -1.0 .and. auto_corr(i) .le. 1.0) then
      auto_corr_sum = auto_corr_sum + auto_corr(i)
      auto_cnt = auto_cnt + 1
    endif
    if(t_p_corr(i) .ge. -1.0 .and. t_p_corr(i) .le. 1.0) then
      tp_corr_sum = tp_corr_sum + t_p_corr(i)
      tp_cnt = tp_cnt + 1
    endif
!print *,'here',i

    deallocate(stn_miss_t)
!print *,'here1b'
    deallocate(stn_miss)
!print *,'here2'
    deallocate(stn_vals)
!print *,'here3'
    deallocate(stn_tair)
!print *,'here4'
  enddo   !end station read loop


  error  = 0

  mean_autocorr = auto_corr_sum/real(auto_cnt,kind(dp))
  mean_tp_corr  = tp_corr_sum/real(tp_cnt,kind(dp))

  print *,'Temp lag-1 autocorrelation: ',mean_autocorr(1)
  print *,'Temp-precip correlation: ', mean_tp_corr(1)


  !pull weight generation outside of time loop


  print *,'Generating base weight matrix & '
  print *,'finding nearest stations for each gridpoint'
  call system_clock (t1, count_rate)


  do g = 1, ngrid, 1


    close_count(g) = 1
    min_weight = 0.0d0
    close_weights(g,:) = 0.0d0

    close_count_t(g) = 1
    min_weight_t = 0.0d0
    close_weights_t(g,:) = 0.0d0

    do i = 1, nstns, 1
      !setup distinct weight matrices for precip and temperature
      call calc_distance_weight(search_distance, X(i,2), X(i,3), Z(g,2), Z(g,3), w_base(g,i))

!original call
!      call calc_distance_weight(maxDistance, X(i,2), X(i,3), Z(g,2), Z(g,3), w_temp(i,i))
      
  !also set some logic to limit the number of stations to the N closest
	min_weight = 0.0d0
	
	if(w_base(g,i) .gt. min_weight .and. stn_data(i,1) .gt. -100.0d0) then
	  if(close_count(g) .le. sta_limit) then

	    close_weights(g,close_count(g)) = w_base(g,i)
	    close_loc(g,close_count(g))  = i

	    close_meta(1,g,close_count(g)) = X(i,2)
	    close_meta(2,g,close_count(g)) = X(i,3)
	    close_meta(3,g,close_count(g)) = Z(g,2)
	    close_meta(4,g,close_count(g)) = Z(g,3)
	    call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta(5,g,close_count(g)))

	    close_count(g) = close_count(g) + 1
	  else
	    min_weight = minval(close_weights(g,:),1)
	    if(w_base(g,i) .gt. min_weight) then
	      out_loc = minloc(close_weights(g,:),1)
	      close_weights(g,out_loc) = w_base(g,i)
	      close_loc(g,out_loc) = i
	      
	      close_meta(1,g,out_loc) = X(i,2)
	      close_meta(2,g,out_loc) = X(i,3)
	      close_meta(3,g,out_loc) = Z(g,2)
	      close_meta(4,g,out_loc) = Z(g,3)
	      call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta(5,g,out_loc))
	    endif
	    
	  endif
	endif

  !need to repeat above for temperature since that data is independent of precipitation
	min_weight_t = 0.0d0

	if(w_base(g,i) .gt. min_weight_t .and. tair_data(1,i,1) .gt. -200.0d0) then
	  if(close_count_t(g) .le. sta_limit) then
	    
	    close_weights_t(g,close_count_t(g)) = w_base(g,i)
	    close_loc_t(g,close_count_t(g))  = i

	    close_meta_t(1,g,close_count_t(g)) = X(i,2)
	    close_meta_t(2,g,close_count_t(g)) = X(i,3)
	    close_meta_t(3,g,close_count_t(g)) = Z(g,2)
	    close_meta_t(4,g,close_count_t(g)) = Z(g,3)
	    call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta_t(5,g,close_count_t(g)))

	    close_count_t(g) = close_count_t(g) + 1
	  else
	    min_weight_t = minval(close_weights(g,:),1)
	    if(w_base(g,i) .gt. min_weight) then
	      out_loc_t = minloc(close_weights_t(g,:),1)
	      close_weights_t(g,out_loc_t) = w_base(g,i)

	      close_loc_t(g,out_loc_t) = i

	      close_meta_t(1,g,out_loc_t) = X(i,2)
	      close_meta_t(2,g,out_loc_t) = X(i,3)
	      close_meta_t(3,g,out_loc_t) = Z(g,2)
	      close_meta_t(4,g,out_loc_t) = Z(g,3)
	      call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta(5,g,out_loc_t))
	    endif
	  endif
	endif

    enddo  !end station loop


  enddo  !end grid point loop
call system_clock (t2, count_rate)
  print *,'Elapsed time for weight generation: ',real( t2 - t1 )/real(count_rate)


!  do t = 4, 4, 4
  do t = 1, ntimes, 1

     call system_clock(tg1,count_rate)
     print *, "TIME STEP = ", Times(t), " (",t,"/",ntimes,")"


     do i = 1, nstns, 1

      Y(i) = stn_data(i,t)

       Y_tmean(i) = (tair_data(1,i,t)+tair_data(2,i,t))/2.0d0 
       Y_trange(i) = abs(tair_data(2,i,t)-tair_data(1,i,t))

     enddo

!do power transformation on input vector
     call normalize_Y(4.0d0,Y)


     do g = 1, ngrid, 1

! call system_clock(tg1,count_rate)

	if(Z(g,4) .gt. -200) then

	  allocate(TWX_red(6,sta_limit))
	  allocate(TX_red(6,sta_limit))
	  allocate(TWX_red_t(6,sta_limit))
	  allocate(TX_red_t(6,sta_limit))

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
	  do i = 1,(close_count(g)-1)
	    if(close_meta(5,g,i) .gt. max_distance) then
	      max_distance = close_meta(5,g,i)
	    endif
	  enddo


	  
	  if(max_distance .le. maxDistance) then
	    expand_dist(g) = max_distance
	    expand_flag(g) = 0
	    max_distance = maxDistance  
	  else
	    max_distance = max_distance+1.0d0
	    expand_flag(g) = 1
	    expand_dist(g) = max_distance
	  endif

    !reduced matrices for precip
          slope_flag_pcp = 0
          do i=1,(close_count(g)-1)
	    call calc_distance_weight(max_distance,close_meta(1,g,i),close_meta(2,g,i), &
                                 close_meta(3,g,i),close_meta(4,g,i),tmp_weight )

	     w_pcp_red(i,i) = tmp_weight
	     w_pcp_1d(i) = tmp_weight
	     w_pcp_1d_loc(i) = close_loc(g,i)
	     Y_red(i) = Y(close_loc(g,i))
	     X_red(i,:) = X(close_loc(g,i),:)



	     if(stn_data(close_loc(g,i),t) .gt. 0.0) then
		ndata = ndata + 1
		Yp_red(i) = 1
	     else
		nodata = nodata + 1
             endif


	  enddo


	  call normalize_Xv(Y_red,w_pcp_1d,step_mean,step_std,step_std_all,step_min,step_max,Yp_red)

	  Y_mean(g,t) = step_mean
	  Y_std(g,t) = step_std
	  Y_std_all(g,t) = step_std_all
	  Y_min(g,t) = step_min
	  Y_max(g,t) = step_max


	  ndata_t = 0
          nodata_t = 0
	  w_temp_red = 0.0
          Y_tmean_red = 0.0
	  Y_trange_red = 0.0
          X_red_t = 0.0


	 ! max_distance_t = maxval(close_meta_t(5,g,:))
	  max_distance_t = 0.0d0
	  do i = 1,(close_count_t(g)-1)
	    if(close_meta_t(5,g,i) .gt. max_distance_t) then
	      max_distance_t = close_meta_t(5,g,i)
	    endif
	  enddo

	  if(max_distance_t .le. maxDistance) then
	    max_distance_t = maxDistance
	  else
	    max_distance_t = max_distance_t+1.0d0
	    expand_flag_t(g) = 1
	    expand_dist_t(g) = max_distance_t
	  endif

    !reduced matrices for temperature
	  slope_flag_temp = 0
          do i=1,(close_count_t(g)-1)
	    call calc_distance_weight(max_distance_t,close_meta_t(1,g,i),close_meta_t(2,g,i), &
                                 close_meta_t(3,g,i),close_meta_t(4,g,i), tmp_weight )

	     w_temp_red(i,i) = tmp_weight
	     w_temp_1d(i)    = tmp_weight
	     Y_tmean_red(i) = Y_tmean(close_loc_t(g,i))
	     Y_trange_red(i) = Y_trange(close_loc_t(g,i))
	     X_red_t(i,:) = X(close_loc_t(g,i),:)

	     if(y_tmean(close_loc_t(g,i)) .gt. -100.0) then
		ndata_t = ndata_t + 1
             else
                nodata_t = nodata_t + 1
             endif

	    
	  enddo

	if(ndata == 0 .AND. nodata == 0) then
	    !print *, "No stations within max distance of grid cell!"
	    POP(g,t) = 0.0
	    PCP(g,t) = 0.0
	    PCPERR(g,t) = 0.0
	    
	    POP_2(g,t) = 0.0
	    PCP_2(g,t) = 0.0
	    PCPERR_2(g,t) = 0.0

	endif

	if(ndata_t == 0 .and. nodata_t == 0) then
	  if(t .gt. 1) then
	    tmean(g,t) = tmean(g,t-1)
	    trange(g,t) = trange(g,t-1)
	    tmean_err(g,t) = tmean_err(g,t-1)
	    trange_err(g,t) = trange_err(g,t-1)

	    tmean_2(g,t) = tmean_2(g,t-1)
	    trange_2(g,t) = trange_2(g,t-1)
	    tmean_err_2(g,t) = tmean_err_2(g,t-1)
	    trange_err_2(g,t) = trange_err_2(g,t-1)
	  else
	    tmean(g,t) = -999
	    trange(g,t) = -999
	    tmean_err(g,t) = 0.0
	    trange_err(g,t) = 0.0

	    tmean_2(g,t) = -999
	    trange_2(g,t) = -999
	    tmean_err_2(g,t) = 0.0
	    trange_err_2(g,t) = 0.0
	  endif
	endif

! print *,ndata

	if(ndata >= 1) then

	    !original call
!	    TWX = matmul(TX, w_pcp)

            !tmp needs to be matmul(TX,X) where TX = TWX_red and X = X_red
	    TWX_red = matmul(transpose(X_red),w_pcp_red)
	    tmp = matmul(TWX_red, X_red)
	    vv=maxval(abs(tmp),dim=2)

	    if (any(vv == 0.0)) then
	      slope_flag_pcp = 0
	    else
	      slope_flag_pcp = 1
	    endif

	    if(nodata == 0) then
	      !print *, "All stations have precip, POP = 1.0"
	      POP(g,t) = 1.0
	
	      POP_2(g,t) = 1.0

	    else



	    if (slope_flag_pcp .eq. 0) then
	      POP(g,t) = -999.
	    else
	      !regression with slope
	      TX_red = transpose(X_red)
	      TWX_red = matmul(TX_red,w_pcp_red)

	      call logistic_regression(X_red,Y_red,TWX_red, Yp_red, B) 
	      POP(g,t) = real(1.0 / (1.0 + exp(-dot_product(Z(g,:),B)) ),kind(sp))

	      deallocate(B)
            endif


	      !regression without slope
	      deallocate(TWX_red)
	      deallocate(TX_red)
	      allocate(TWX_red(4,sta_limit))
	      allocate(TX_red(4,sta_limit))
	      TX_red = transpose(X_red(:,1:4))

	      TWX_red = matmul(TX_red,w_pcp_red)

	      

	      call logistic_regression(X_red(:,1:4),Y_red,TWX_red, Yp_red, B)
	      POP_2(g,t) = real(1.0 / (1.0 + exp(-dot_product(Z(g,1:4),B)) ),kind(sp))


	      deallocate(B)

	    endif

	    !print *, "POP: ", POP(g,t)


	    deallocate(TWX_red)
	    deallocate(TX_red)
	    allocate(TWX_red(6,sta_limit))
	    allocate(TX_red(6,sta_limit))


	    if (slope_flag_pcp .eq. 0) then
	      PCP(g,t) = -999.
            else
	      !regression with slope
	      TX_red = transpose(X_red)
	      TWX_red = matmul(TX_red,w_pcp_red)

	      call least_squares(X_red, Y_red, TWX_red, B)
	      PCP(g,t) = real(dot_product(Z(g,:),B),kind(sp))
!	      

	      wgtsum = 0.0
	      errsum = 0.0
	      ss_tot = 0.0
	      ss_res = 0.0
	      do i = 1,(close_count(g)-1),1
		wgtsum = wgtsum + w_pcp_red(i,i)


		errsum = errsum + (w_pcp_red(i,i) * (PCP(g,t) - Y_red(i))**2 )

	      enddo
	      PCPERR(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))

	    endif


      !regression without slope
	      
	    deallocate(TWX_red)
	    deallocate(TX_red)
	    allocate(TWX_red(4,sta_limit))
	    allocate(TX_red(4,sta_limit))

	    TX_red = transpose(X_red(:,1:4))
	    TWX_red = matmul(TX_red,w_pcp_red)

	    call least_squares(X_red(:,1:4), Y_red, TWX_red, B)
	    PCP_2(g,t) = real(dot_product(Z(g,1:4),B),kind(sp))

	    wgtsum = 0.0
	    errsum = 0.0
	    ss_tot = 0.0
	    ss_res = 0.0
	    do i = 1,(close_count(g)-1),1
	      wgtsum = wgtsum + w_pcp_red(i,i)

	      errsum = errsum + (w_pcp_red(i,i) * (PCP_2(g,t) - Y_red(i))**2 )

	    enddo
	    PCPERR_2(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))


	    deallocate(B)

	else
	    !print *, "Not enough stations with data within max distance"
	    POP(g,t) = 0.0
	    PCP(g,t) = 0.0
	    PCPERR(g,t) = 0.0

	    POP_2(g,t) = 0.0
	    PCP_2(g,t) = 0.0
	    PCPERR_2(g,t) = 0.0
	endif  !precip if statement



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!  do temperature ols now
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(ndata_t .ge. 1) then 


	  !regression with slope
	  TX_red_t = transpose(X_red_t)
	  TWX_red_t = matmul(TX_red_t, w_temp_red)
	  call least_squares(X_red_t, y_tmean_red, TWX_red_t, B)
	  tmean(g,t) = real(dot_product(Z(g,:),B),kind(sp))
	  
	  errsum = 0.0
	  wgtsum = 0.0
	  do i = 1,(close_count_t(g)-1),1
	      wgtsum = wgtsum + w_temp_red(i,i)

	      errsum = errsum + (w_temp_red(i,i) * (tmean(g,t) - Y_tmean_red(i))**2)
	  enddo
	  tmean_err(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))


	  !regression without slope
	  deallocate(B)
	  deallocate(TWX_red_t)
	  deallocate(TX_red_t)
	  allocate(TWX_red_t(4,sta_limit))
	  allocate(TX_red_t(4,sta_limit))
	  TX_red_t = transpose(X_red_t(:,1:4))
	  TWX_red_t = matmul(TX_red_t, w_temp_red)

	  call least_squares(X_red_t(:,1:4), y_tmean_red, TWX_red_t, B)
	  tmean_2(g,t) = real(dot_product(Z(g,1:4),B),kind(sp))

	  errsum = 0.0
	  wgtsum = 0.0
	  !do i = 1, nstns, 1
	  do i = 1,(close_count_t(g)-1),1
	      wgtsum = wgtsum + w_temp_red(i,i)


	      errsum = errsum + (w_temp_red(i,i) * (tmean_2(g,t) - Y_tmean_red(i))**2)

	  enddo
	  tmean_err_2(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))

	  deallocate(B)

	  !then do trange


	  !  print *,'OLS FOR TRANGE'
	  !  print *,'**********************'



	  deallocate(TWX_red_t)
	  deallocate(TX_red_t)
	  allocate(TWX_red_t(6,sta_limit))
	  allocate(TX_red_t(6,sta_limit))
	  TX_red_t = transpose(X_red_t)
	  TWX_red_t = matmul(TX_red_t, w_temp_red)

	  call least_squares(X_red_t, y_trange_red, TWX_red_t, B)
	  trange(g,t) = real(dot_product(Z(g,:),B),kind(sp))
	  
	  errsum = 0.0
	  wgtsum = 0.0
	  !do i = 1, nstns, 1
	  do i = 1,(close_count_t(g)-1),1
	      wgtsum = wgtsum + w_temp_red(i,i)

	      errsum = errsum + (w_temp_red(i,i) * (trange(g,t) - Y_trange_red(i))**2)

	  enddo
	  trange_err(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))

!print *,'trange no slope'

	  !regression without slope
	  deallocate(B)
	  deallocate(TWX_red_t)
	  deallocate(TX_red_t)
	  allocate(TWX_red_t(4,sta_limit))
	  allocate(TX_red_t(4,sta_limit))
	  TX_red_t = transpose(X_red_t(:,1:4))
	  TWX_red_t = matmul(TX_red_t, w_temp_red)

	  call least_squares(X_red_t(:,1:4), y_trange_red, TWX_red_t, B)
	  trange_2(g,t) = real(dot_product(Z(g,1:4),B),kind(sp))

	  errsum = 0.0
	  wgtsum = 0.0
	  !do i = 1, nstns, 1
	  do i = 1,(close_count_t(g)-1),1
	      wgtsum = wgtsum + w_temp_red(i,i)

	      sta_temp = real(dot_product(X_red_t(i,1:4),B),kind(sp))

	      errsum = errsum + (w_temp_red(i,i) * (trange_2(g,t) - Y_trange_red(i))**2)

	  enddo
	  trange_err_2(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))

	else

	    !if not enough stations with data
	    !just use value from previous grid point for now
            if(g .gt. 1) then
	      trange(g,t) = trange(g-1,t)
	      trange_err(g,t) = trange_err(g-1,t)
	      tmean(g,t) = tmean(g-1,t)
	      tmean_err(g,t) = tmean_err(g-1,t)

	      trange_2(g,t) = trange_2(g-1,t)
	      trange_err_2(g,t) = trange_err_2(g-1,t)
	      tmean_2(g,t) = tmean_2(g-1,t)
	      tmean_err_2(g,t) = tmean_err_2(g-1,t)
	    else
	      trange(g,t) = trange(g,t-1)
	      trange_err(g,t) = trange_err(g-1,t-1)
	      tmean(g,t) = tmean(g,t-1)
	      tmean_err(g,t) = tmean_err(g,t-1)

	      trange_2(g,t) = trange_2(g,t-1)
	      trange_err_2(g,t) = trange_err_2(g-1,t-1)
	      tmean_2(g,t) = tmean_2(g,t-1)
	      tmean_err_2(g,t) = tmean_err_2(g,t-1)
	    endif
	endif !end data check if statement for temperature


      endif !end check for valid elevation

     enddo !end grid loop
call system_clock(tg2,count_rate)
  print *,'Elapsed time for one time step: ',real( tg2 - tg1 )/real(count_rate)
  enddo  !end time loop

end subroutine estimate_precip


subroutine normalize_X(X)
  use type
  implicit none

  real(DP), intent(inout) :: X(:,:)

  real(DP) :: mean, stdev, sum_x, sum_x2
  integer(I4B) :: v, t
  integer(I4B) :: nvars, ntimes

  nvars = size(X,2)-1
  ntimes = size(X,1)

  do v = 2, nvars+1, 1
    sum_x = 0.0d0
    sum_x2 = 0.0d0
    do t = 1, ntimes, 1
      sum_x  = sum_x + x(t,v)
      sum_x2 = sum_x2 + x(t,v)**2 
    enddo
    mean = sum_x / real(ntimes)
    stdev = sqrt(  (real(ntimes) * sum_x2 - sum_x**2) /  (real(ntimes)*real(ntimes-1)) )
    do t = 1, ntimes, 1
       if(stdev .eq. 0.0) then
          x(t,v) = x(t,v)
       else
          x(t,v) = (x(t,v) - mean) / stdev
       endif
    enddo
  enddo

end subroutine normalize_X


subroutine normalize_Xv(X,weight,mean,stdev,stdev_all,smin,smax,Yp)
  use type
  implicit none

  real(DP), intent(inout) :: X(:)
  real(DP), intent(in)	:: weight(:)
  real(DP), intent(out) :: mean
  real(DP), intent(out) :: stdev
  real(DP), intent(out) :: stdev_all
  real(DP), intent(out) :: smin
  real(DP), intent(out) :: smax
  integer(I4B),intent(out) :: Yp(:)

!  real(DP) :: mean, stdev, sum_x2, sum_x
  real(DP) :: sum_x2, sum_x
  real(DP) :: sum_stdev,sum_std

  real(DP) :: mean_all,sum_weight,sum_xw

  integer(I4B) :: t
  integer(I4B) :: ntimes

  ntimes = size(X)
  
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
  do t = 1, ntimes, 1
    if(x(t) > 0.0) then
      sum_x  = sum_x + x(t)
      sum_x2 = sum_x2 + x(t)**2 
      sum_xw = sum_xw + weight(t)*x(t)
      sum_weight = sum_weight + weight(t)
      Yp(t) = 1

      if(x(t) .le. smin) then
	smin = x(t)
      endif
      if(x(t) .ge. smax) then
        smax = x(t)
      endif

    endif
  enddo

  mean_all = sum_x / real(ntimes)
  mean = sum_x / real(sum(Yp))
  do t = 1,ntimes,1
    sum_stdev = sum_stdev+(x(t)-mean_all)**2
    if(Yp(t) .eq. 1) then
      sum_std = sum_std + (x(t)-mean)**2
    endif
  enddo
!  stdev_all = sqrt(sum_stdev/(ntimes-1))
  

  if(sum(Yp) .ge. 2) then
   stdev = sqrt(sum_std/real(sum(Yp)-1.0))     
   stdev_all = sqrt(sum_std/real(ntimes-1.0))  

  else
    mean = sum_x
    stdev = 0.01
!    stdev_all = sum_xw
    stdev_all = 0.01
  endif


  if(stdev .gt. 0.0) then
    do t = 1, ntimes, 1
      if(Yp(t) .gt. 0.0) then
	x(t) = (x(t) - mean) / stdev
      endif
    enddo
  endif


  return

end subroutine normalize_Xv


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
subroutine normalize_Y(texp,Y)
  use type
  implicit none

  real(DP), intent(in)	  :: texp  !transform exponent
  real(DP), intent(inout) :: Y(:)
  integer(I4B) :: t
  integer(I4B) :: ntimes

  ntimes = size(Y)

  do t = 1, ntimes, 1
!     Y(t) = Y(t) ** (1.0d0/2.5d0)
!    Y(t) = Y(t) ** (1.0d0/4d0)
    Y(t) = Y(t) ** (1.0d0/texp)
  end do

end subroutine normalize_Y

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
subroutine calc_weights(Times, tt, X, W)
  use type
  implicit none
  
  real(DP), intent(in) :: Times(:)
  integer(I4B), intent(in) :: tt
  real(DP), intent(in) :: X(:,:)
  real(DP), allocatable, intent(out) :: W(:,:)
  real(DP) :: sum
  integer(I4B) :: v, t
  integer(I4B) :: nvars, ntimes

  nvars = size(X,2) -1
  ntimes = size(X,1)

  allocate(W(ntimes,ntimes))

  do t = 1, ntimes, 1
     W(t,:) = 0.0d0
     if(t /= tt) then
        sum = 0.0d0
        do v = 2, nvars+1, 1
           !print *, "X(t,v) : ",  X(t,v)
           sum = sum + (X(t,v) - X(tt,v))**2
        enddo
        W(t,:) = 0.0d0

	W(t,t) = 1/sqrt(sum / nvars)
     endif
!     print *, "W : ",  t, W(t,t)
  enddo

end subroutine calc_weights

! Great circle distance calculation
! Output in nm
subroutine calc_distance_weight(maxd, lat1, lon1, lat2, lon2, weight)
  use type
  implicit none

  real(DP), intent(in) :: maxd, lat1, lon1, lat2, lon2
  real(DP), intent(out) :: weight

  real(DP) :: dist, lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927

  lat1r = lat1*Pi/180.0d0;
  lon1r = lon1*Pi/180.0d0;
  lat2r = lat2*Pi/180.0d0;
  lon2r = lon2*Pi/180.0d0;
  dist =  ((180*60)/Pi)*(2*asin(sqrt( (sin((lat1r-lat2r)/2))**2 +  &
		     cos(lat1r)*cos(lat2r) * (sin((lon1r-lon2r)/2))**2 )));
  if(dist .gt. maxd) then
     weight = 0.0d0
  else
     weight = (1.0d0- (dist / maxd)**3 )**3
!    weight = 1.0d0 - (dist/maxd)**0.5
  endif

end subroutine calc_distance_weight

! Great circle distance calculation
! Output in nm
subroutine calc_distance(lat1, lon1, lat2, lon2, dist)
  use type
  implicit none

  real(DP), intent(in) :: lat1, lon1, lat2, lon2
  real(DP), intent(out) :: dist

  real(DP) :: lat1r, lon1r, lat2r, lon2r
  !real(DP) :: Pi
  !Pi = 3.1415927

  lat1r = lat1*Pi/180.0d0;
  lon1r = lon1*Pi/180.0d0;
  lat2r = lat2*Pi/180.0d0;
  lon2r = lon2*Pi/180.0d0;                                                           
  dist =  ((180*60)/Pi)*(2*asin(sqrt( (sin((lat1r-lat2r)/2))**2 +  &
		     cos(lat1r)*cos(lat2r) * (sin((lon1r-lon2r)/2))**2 )));

end subroutine calc_distance

!
! Solve linear equation for x (Ax = b => x = bA^-1) using LU decomposition and back substitution.
! Input:
!   X  = An m by n array.
!   TX = Precalculated transpose array of X, size n by m
!   Y  = An m-element vector containing the right-hand side of the linear system Ax = b. 
! Output: 
!   B  = An n-element vector.
subroutine least_squares(X, Y, TX, B)
  use type
  implicit none

  interface
     subroutine ludcmp(a,indx,d)
       use type
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
       REAL(DP), INTENT(OUT) :: d
     end subroutine ludcmp

     subroutine lubksb(a,indx,b)
       USE type
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
     end subroutine lubksb
  end interface

  real(DP), intent(in) :: X(:,:)
  real(DP), intent(in) :: Y(:)
  real(DP), intent(in) :: TX(:,:)
  real(DP), allocatable, intent(out) :: B(:)

  real(DP), allocatable :: A(:,:)
  integer(I4B), allocatable :: indx(:)
  integer(I4B) :: nvars, ntimes
  real(DP) :: d

  nvars = size(X,2) -1
  ntimes = size(Y)

  allocate(B(nvars+1))
  allocate(A(nvars+1,nvars+1))
  allocate(indx(nvars+1))

  B = matmul(TX, Y)
  A = matmul(TX, X)

  call ludcmp(A, indx, d )
  if (ANY(ABS(A) < 9.99999968E-15)) then
     B(:) = 0.0d0
!     print *, "Warning, LUdcmp produced a zero."
     return
  endif

  call lubksb(A, indx, B)

  deallocate(A)
  deallocate(indx)

end subroutine least_squares


subroutine logistic_regression(X, Y, TX, Yp, B)
  use type
  implicit none

  interface
     subroutine least_squares(X, Y, TX, B)
       use type
       real(DP), intent(in) :: X(:,:)
       real(DP), intent(in) :: Y(:)
       real(DP), intent(in) :: TX(:,:)
       real(DP), allocatable, intent(out) :: B(:)
     end subroutine least_squares
  end interface

  real(DP), intent(in) :: X(:,:)
  real(DP), intent(in) :: Y(:)
  real(DP), intent(in) :: TX(:,:)
  integer(I4B),intent(in) :: Yp(:)
  real(DP), allocatable, intent(out) :: B(:)

  real(DP), allocatable :: Ypd(:), P(:), YN(:), BN(:), V(:,:), XV(:,:)
!  real(DP), allocatable :: P(:), YN(:), BN(:), V(:,:), XV(:,:)
  integer(I4B) :: nvars, ntimes, i, t, f, it
  real(DP) :: d

  nvars = size(X,2) -1
  ntimes = size(Y)

  allocate(B(nvars+1))
  allocate(Ypd(ntimes))
  allocate(YN(ntimes))
  allocate(P(ntimes))
  allocate(V(ntimes,ntimes))
  allocate(XV(ntimes, nvars+1))

  do t = 1, ntimes, 1
    if(Yp(t) .gt. 0.0) then
      Ypd(t) = 1.0d0
    else
      Ypd(t) = 0.0d0
    end if
  end do

  B = 0.0d0
  i = 0
  it = 0

  do while (f /= 1)
!     print *, "Iteration ", it
     P = 1.0d0 / (1.0d0 + exp(-matmul(X, B))) 
     if (ANY(P > 0.97)) then
        !print *, "WARNING: logistic regression diverging"
        f = 1
     else

        YN = Ypd - P
        V = 0.0d0
        do t = 1, ntimes, 1
           V(t,t) = P(t)*(1.0d0-P(t))
        enddo
        XV = matmul(V,X)
        call least_squares(XV, YN, TX, BN)
       
        f = 1
        do i = 1, nvars+1, 1
           if(BN(i) .GT. 1.0E-04 .OR. BN(i) .LT. -1.0E-04) then
              f = 0
           end if
        end do
        if(it > 8) then
           !print *, "WARNING: logistic regression failed to converge"
           f = 1
        endif
        
        B = B + BN
!        print *, "Bnew = ", B
        deallocate(BN)
        
     endif
     it = it+1
  end do
!  print *, "Iterations = ", it
  !print *, "Final B = ", B

end subroutine logistic_regression


subroutine logistic_regressionrf(X, Y, TX, B)
  use type
  implicit none

  interface
     subroutine least_squares(X, Y, TX, B)
       use type
       real(DP), intent(in) :: X(:,:)
       real(DP), intent(in) :: Y(:)
       real(DP), intent(in) :: TX(:,:)
       real(DP), allocatable, intent(out) :: B(:)
     end subroutine least_squares
  end interface

  real(DP), intent(in) :: X(:,:)
  real(DP), intent(in) :: Y(:)
  real(DP), intent(in) :: TX(:,:)
  real(DP), allocatable, intent(out) :: B(:)

  real(DP), allocatable :: Ypd(:), P(:), YN(:), BN(:), V(:,:), XV(:,:)
  integer(I4B) :: nvars, ntimes, i, t, f, it
  real(DP) :: d

  nvars = size(X,2) -1
  ntimes = size(Y)

  allocate(B(nvars+1))
  allocate(Ypd(ntimes))
  allocate(YN(ntimes))
  allocate(P(ntimes))
  allocate(V(ntimes,ntimes))
  allocate(XV(ntimes, nvars+1))

  do t = 1, ntimes, 1
    if(Y(t) .ne. 0.0) then
      Ypd(t) = 1.0d0
    else
      Ypd(t) = 0.0d0
    end if
  end do
 
  B = 0.0d0
  i = 0
  it = 0
  !print *, "B = ", B
  do while (f /= 1)
!     print *, "Iteration ", it
     P = 1.0d0 / (1.0d0 + exp(-matmul(X, B))) 
!     print *, "Pie = ", P
     if (ANY(P > 0.97)) then
        !print *, "WARNING: logistic regression diverging"
        f = 1
     else

        YN = Ypd - P
        V = 0.0d0
        do t = 1, ntimes, 1
           V(t,t) = P(t)*(1.0d0-P(t))
        enddo
        XV = matmul(V,X)

        call least_squares(XV, YN, TX, BN)
        
        f = 1
        do i = 1, nvars+1, 1
           if(BN(i) .GT. 1.0E-04 .OR. BN(i) .LT. -1.0E-04) then
              f = 0
           end if
        end do
        if(it > 8) then
           !print *, "WARNING: logistic regression failed to converge"
           f = 1
        endif
        
        B = B + BN
        deallocate(BN)
        
     endif
     it = it+1
  end do


end subroutine logistic_regressionrf

subroutine generic_corr(stn_data,tair_data,lag,window,auto_corr,t_p_corr)
  use type

  implicit none

!input
  real(DP),intent(in)       :: stn_data(:)
  real(DP),intent(in)       :: tair_data(:,:)
  integer(I4B),intent(in)   :: lag
  integer(I4B),intent(in)   :: window


!output
  real(DP),intent(out)   :: auto_corr
  real(DP),intent(out)   :: t_p_corr

!local variables
  real(DP),allocatable  :: tmean(:),trange(:)
  real(DP),allocatable  :: moving_avg(:,:)
  real(DP)              :: lag_0_mean
  real(DP)              :: lag_n_mean
  real(DP)              :: lag_0_var
  real(DP)              :: lag_n_var
  real(DP)              :: lag_0_sum
  real(DP)              :: lag_n_sum
  real(DP)              :: cov
  real(DP)              :: lag_0_pmean,lag_0_pvar,lag_0_psum
  real(DP)              :: trange_mean,trange_sum,trange_var

  real(DP)              :: tmp_tmean, tmp_trange

  integer(I4B)          :: i,j,tmp_cnt
  integer(I4B)          :: ntimes
  integer(I4B)          :: half_window
  integer(I4B)          :: cnt_sums
  integer(I4B)          :: data_cnt

!code
  ntimes = size(stn_data)

  allocate(tmean(ntimes))
  allocate(trange(ntimes))
  allocate(moving_avg(2,ntimes))

  data_cnt = 0
  do i = 1,ntimes,1
    if(tair_data(1,i) .gt. -100.0 .and. tair_data(2,i) .gt. -100.0) then
      tmean(i) = ((tair_data(1,i)+tair_data(2,i))/2.d0)+273.15d0
      trange(i) = (tair_data(2,i)-tair_data(1,i)/2.d0)
      data_cnt = data_cnt + 1
    else
      tmean(i)= -999.0d0
      trange(i) = -999.0d0
    endif
  enddo


  half_window = floor(window/2.0d0)

!do the lag correlation for temperature
!first compute the moving average for climo removal
!need to check for missing values....

  do i = 1,ntimes,1
    if(i .lt. half_window) then
      tmp_tmean = 0.0
      tmp_trange = 0.0
      tmp_cnt = 0
      do j = 1,window,1
	if(tmean(j) .gt. -100.0) then
	  tmp_tmean = tmp_tmean + tmean(j)
          tmp_trange = tmp_trange + trange(j)
          tmp_cnt = tmp_cnt + 1
	endif
	if(tmp_cnt .gt. 0) then
	  moving_avg(1,i) = tmp_tmean/real(tmp_cnt,kind(dp))
	  moving_avg(2,i) = tmp_trange/real(tmp_cnt,kind(dp))
	else
	  moving_avg(1,i) = -999.0
	  moving_avg(2,i) = -999.0
	endif
      enddo

    elseif(i .gt. ntimes-half_window) then
      tmp_tmean = 0.0
      tmp_trange = 0.0
      tmp_cnt = 0
      do j = 1,window,1
      	if(tmean(ntimes-j) .gt. -100.0) then
	  tmp_tmean = tmp_tmean + tmean(ntimes-j)
          tmp_trange = tmp_trange + trange(ntimes-j)
          tmp_cnt = tmp_cnt + 1
	endif
      enddo
	if(tmp_cnt .gt. 0) then
	  moving_avg(1,i) = tmp_tmean/real(tmp_cnt,kind(dp))
	  moving_avg(2,i) = tmp_trange/real(tmp_cnt,kind(dp))
	else
	  moving_avg(1,i) = -999.0
	  moving_avg(2,i) = -999.0
	endif
    else
      tmp_tmean = 0.0
      tmp_trange = 0.0
      tmp_cnt = 0
      do j = -half_window,half_window,1
      	if(tmean(i+j) .gt. -100.0) then
	  tmp_tmean = tmp_tmean + tmean(i+j)
          tmp_trange = tmp_trange + trange(i+j)
          tmp_cnt = tmp_cnt + 1
	endif
      enddo
	if(tmp_cnt .gt. 0) then
	  moving_avg(1,i) = tmp_tmean/real(tmp_cnt,kind(dp))
	  moving_avg(2,i) = tmp_trange/real(tmp_cnt,kind(dp))
	else
	  moving_avg(1,i) = -999.0
	  moving_avg(2,i) = -999.0
	endif
    endif
  enddo


!only use portions of timeseries to compute auto_corr
!need to go through and check to see if values exist for lag-0 and lag-n and moving_avg
!if values do not exist for any of the three, don't add to running sums

!compute means
  lag_0_sum  = 0.0d0
  lag_n_sum  = 0.0d0
  lag_0_var  = 0.0d0
  lag_n_var  = 0.0d0
  cov        = 0.0d0
  cnt_sums = 0

  do i = lag+1,ntimes,1
    if(tmean(i) .gt. -100.0 .and. tmean(i-lag) .gt. -100.0 .and. &
    moving_avg(1,i) .gt. -100.0 .and. moving_avg(1,i-lag) .gt. -100.0) then
      lag_n_sum = lag_n_sum + tmean(i-lag)-moving_avg(1,i-lag)
      lag_0_sum = lag_0_sum + tmean(i)-moving_avg(1,i)
      cnt_sums = cnt_sums + 1
    endif
  enddo

  lag_0_mean = lag_0_sum/real(cnt_sums,kind(dp))
  lag_n_mean = lag_n_sum/real(cnt_sums,kind(dp))

!compute variance,covariance
  do i = lag+1,ntimes,1
    if(tmean(i) .gt. -100.0 .and. tmean(i-lag) .gt. -100.0 .and. &
    moving_avg(1,i) .gt. -100.0 .and. moving_avg(1,i-lag) .gt. -100.0) then
      lag_n_var = lag_n_var + ((tmean(i-lag)-moving_avg(1,i-lag)) -lag_n_mean)**2
      lag_0_var = lag_0_var + ((tmean(i)-moving_avg(1,i)) -lag_0_mean)**2
      cov       = cov +       ((tmean(i-lag)-moving_avg(1,i-lag)) -lag_n_mean)*((tmean(i)-moving_avg(1,i))-lag_0_mean)
    endif
  enddo

!compute autocorrelation
  auto_corr = cov/(sqrt(lag_0_var)*sqrt(lag_n_var))


!!!!!!!!!!!!!!!!
!
! now do the t - p correlation
! do correlation on trange, not tmean
!
!!!!!!!!!!!!!!!!!!!

  lag_0_pmean = 0.0d0
  lag_0_pvar  = 0.0d0
  lag_0_psum  = 0.0d0
  trange_sum = 0.0d0
  trange_mean= 0.0d0
  trange_var = 0.0d0
  cov        = 0.0d0
  tmp_cnt    = 0


!again need to check for missing values....
  do i = 1,ntimes,1
    if(trange(i) .gt. -100 .and. stn_data(i) .gt. -100.0) then
      !compute for precip mean
      lag_0_psum = lag_0_psum + stn_data(i)
      !anomaly means of trange
      trange_sum = trange_sum + (trange(i)-moving_avg(2,i))
      tmp_cnt = tmp_cnt + 1
    endif
  enddo
  lag_0_pmean = lag_0_psum/real(tmp_cnt,kind(dp))
  trange_mean = trange_sum/real(tmp_cnt,kind(dp))

!compute variance and covariance
  do i = 1,ntimes,1
    if(trange(i) .gt. -100 .and. stn_data(i) .gt. -100.0) then
      lag_0_pvar = lag_0_pvar + (stn_data(i)-lag_0_pmean)**2
      trange_var = trange_var + ((trange(i)-moving_avg(2,i))-trange_mean)**2
      cov       = cov       + ((trange(i)-moving_avg(2,i))-trange_mean)*(stn_data(i)-lag_0_pmean)
    endif
  enddo

  t_p_corr = cov/(sqrt(lag_0_pvar)*sqrt(trange_var))
  if(sqrt(lag_0_pvar)*sqrt(trange_var) .le. 0.00001) then
    t_p_corr = 0.0
  endif

!print *,t_p_corr

!in some situations, there are very limited data used for calculations
!set those cases to missing value
  if(data_cnt .lt. real(ntimes)*0.25) then
    auto_corr = -999.0
    t_p_corr  = -999.0
  endif

end subroutine generic_corr
