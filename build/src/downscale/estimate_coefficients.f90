!subroutine estimate_coefficients (d, nvars, lats, lons, times, st_rec, end_rec, stnid, stnlat, &
!& stnlon, stnalt, stnvar, site_var, site_list, directory, c, poc, error) !AWW added directory

subroutine estimate_coefficients (d, nvars, lats, lons, times, st_rec, end_rec, stnid, stnlat, &
  & stnlon, stnvar,  directory, c, poc, error) !AWW added directory

  ! ==================================================================================================
  ! This routine called in MODE 1 usage:  downscaling gridded data to create station/point ensembles
  ! ==================================================================================================

  use type
  use utim ! AWW-add:  can figure out start & ends records with this for station read
  implicit none
 
  ! ==== interfaces ======
  interface
 
    subroutine read_station (stnvar, stnid, directory, st_rec, end_rec, vals, tair_vals, &
   & vals_miss, vals_miss_t, error)
      use type
      use utim ! AWW
      character (len=100), intent (in) :: stnvar
      character (len=100), intent (in) :: stnid
      character (len=500), intent (in) :: directory ! AWW
      integer (i4b), intent (in) :: st_rec, end_rec ! AWW
      real (dp), allocatable, intent (out) :: vals (:), tair_vals (:, :)
      logical, allocatable, intent (out) :: vals_miss (:), vals_miss_t (:)
      integer, intent (out) :: error
    end subroutine read_station
 
    subroutine normalize_x (x)
      use type
      real (dp), intent (inout) :: x (:, :)
    end subroutine normalize_x
 
    subroutine normalize_y (texp, y)
      use type
      real (dp), intent (in) :: texp !transform exponent
      real (dp), intent (inout) :: y (:)
    end subroutine normalize_y
 
    subroutine calc_weights (tt, x, w)
      use type
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
 
    subroutine logistic_regressionrf (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine logistic_regressionrf
 
  end interface
  ! ======   end of interfaces -- start code  =======
 
  real (dp), intent (in) :: d (:, :, :), lats (:), lons (:)
  real (dp), intent (in) :: times (:)
 
  integer (i4b), intent (in) :: st_rec, end_rec ! AWW
 
  integer (i4b), intent (in) :: nvars
  character (len=100), intent (in) :: stnid (:)
  real (dp), intent (in) :: stnlat (:), stnlon (:)
  character (len=100), intent (in) :: stnvar
  character (len=500), intent (in) :: directory ! AWW-feb2016 added for data location
  real (dp), allocatable, intent (out) :: c (:, :, :), poc (:, :, :)
  integer, intent (out) :: error
 
  !  real(DP), allocatable :: X(:,:), Y(:), XS(:,:), TWXS(:,:), YS(:), W(:,:), B(:),tair_vals(:,:)
  real (dp), allocatable :: x(:, :), y(:), xs(:, :), xp (:, :), twxs(:, :), ys(:), w(:, :)
  real (dp), allocatable :: b(:), tair_vals (:, :)
 
  real (dp), allocatable :: ts(:)
  logical, allocatable :: y_miss(:), y_miss_t(:) ! modified AJN Sept 2013
  integer (i4b) :: ntimes, nlats, nlons, nstns
  integer (i4b) :: i, j, k, t, v, tt
  real (dp) :: minlondis, minlatdis
  integer (i4b) :: minlatk, minlonj, gridindex
  real (dp) :: transform_exp  ! precip. transform variable, the exponent of 
                              !   transform norm_pcp = pcp^(1/transform_exp)
 
  print *, "allocated memory in estimate_coefficients"
  error = 0
  transform_exp = 4.0d0
 
  ntimes = size (times)
  nlats = size (lats)  ! grid lats & lons
  nlons = size (lons)
  nstns = size (stnlat)
  print *, "In estimate_coefficients: ", ntimes, nlats, nlons, nstns
 
  allocate (x(ntimes, nvars+1))
  if (trim(stnvar) .eq. "PRCP") then
    allocate (poc(nstns, ntimes, nvars+1))
    poc = 0.0d0
  end if
  allocate (c(nstns, ntimes, nvars+1))
  c = 0.0d0
 
  ! ==== loop through stations ===
  do i = 1, nstns, 1
 
    minlondis = 360.0
    minlonj = -1
    do j = 1, nlons, 1
      if (abs(stnlon(i)-lons(j)) < minlondis) then
        minlondis = abs (stnlon(i)-lons(j))
        minlonj = j
      end if
    end do
    minlatdis = 180.0
    minlatk = -1
    do k = 1, nlats, 1
      if (abs(stnlat(i)-lats(k)) < minlatdis) then
        minlatdis = abs (stnlat(i)-lats(k))
        minlatk = k
      end if
    end do
 
    if (minlonj == -1 .or. minlatk == -1) then
      print *, "Failed to find closest grid point for station: ", trim (stnid(i))
      error = 1
      return
    end if
 
    !AJN wrong (old version)
    !     gridindex = ((minlonj-1) * nlats) + minlatk
 
    gridindex = ((minlatk-1)*nlons) + minlonj
 
    print *, "Station: ", trim (stnid(i)), stnlat (i), stnlon (i)
    print *, "Closest Grid point: ", gridindex, lats (minlatk), lons (minlonj)
 
    x (:, 1) = 1.0
    do v = 1, nvars, 1
      x (:, v+1) = d (v, gridindex, :)
      do t = 1, ntimes, 1
        if (x(t, v+1) /= 0.0) exit
        if (t == ntimes) then
          print *, "ERROR: var ", v, " is all zero for station ", i
          x (:, v+1) = 1.0
        end if
      end do
        !        print *, "VAR NUM:", v, X(:,v+1)
    end do
 
    ! AWW added directory, st_rec, end_rec
    call read_station (stnvar, stnid(i), directory, st_rec, end_rec, y, tair_vals, y_miss, &
   & y_miss_t, error)
    !     print *, "Y:", Y
 
    if (count(y_miss) < ntimes) then
      allocate (ts(count(y_miss)))
      allocate (xs(count(y_miss), nvars+1))
      allocate (xp(count(y_miss), nvars+1))
      allocate (twxs(nvars+1, count(y_miss)))
      allocate (ys(count(y_miss)))
      ts (:) = pack (times, y_miss)
      do v = 1, nvars + 1, 1
        xs (:, v) = pack (x(:, v), y_miss)
      end do
      xp (:, :) = xs (:, :)
      ys (:) = pack (y, y_miss)
    else
      allocate (ts(ntimes))
      allocate (xs(ntimes, nvars+1))
      allocate (xp(ntimes, nvars+1))
      allocate (twxs(nvars+1, ntimes))
      allocate (ys(ntimes))
      ts (:) = times (:)
      xs (:, :) = x (:, :)
      xp (:, :) = xs (:, :)
      ys (:) = y (:)
    end if
 
    call normalize_x (xs)
    call normalize_y (transform_exp, ys)
 
    tt = 1
    do t = 1, ntimes, 1
 
      if (y_miss(t) .eqv. .true.) then
 
        print *, t, ntimes, y_miss (t)
        ! AJN original version
 
        call calc_weights (tt, xs, w)
        ! call calc_weights(TS, tt, XP, W)
 
        ! AJN original version
        twxs = matmul (transpose(xs), w)
        ! TWXS = matmul(transpose(XP), W)
 
        if (trim(stnvar) .eq. "PRCP") then
 
          ! call logistic_regressionrf(XS, YS, TWXS, B)
          call logistic_regressionrf (xp, ys, twxs, b)
 
          poc (i, t, :) = b (:)
          deallocate (b)
        end if
 
        ! call least_squares(XS, YS, TWXS, B)
        call least_squares (xp, ys, twxs, b)
        c (i, t, :) = b (:)
 
        deallocate (b)
 
        tt = tt + 1
      else
        if (trim(stnvar) .eq. "PRCP") then
          poc (i, t, :) = -999.99
        end if
        c (i, t, :) = -999.99
      end if
 
    end do
 
    deallocate (ys)
    deallocate (twxs)
    deallocate (xs)
    deallocate (ts)
    deallocate (xp)
 
  end do
 
end subroutine estimate_coefficients
