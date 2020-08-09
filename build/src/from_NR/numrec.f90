! This collection of routines from num_rec was put together as part of the
! downscaling code effort.  NCAR 2013

 
subroutine nrerror (string)
  character (len=*), intent (in) :: string
  write (*,*) 'nrerror: ', string
  stop 'program terminated by nrerror'
end subroutine nrerror
 
function assert_eq (n1, n2, n3, string)
  character (len=*), intent (in) :: string
  integer, intent (in) :: n1, n2, n3
  integer :: assert_eq
  if (n1 == n2 .and. n2 == n3) then
    assert_eq = n1
  else
    write (*,*) 'nrerror: an assert_eq failed with this tag:', string
    stop 'program terminated by assert_eq3'
  end if
end function assert_eq
 
subroutine swap_d (a, b)
  use type
  real (dp), dimension (:), intent (inout) :: a, b
  real (dp), dimension (size(a)) :: dum
  dum = a
  a = b
  b = dum
end subroutine swap_d
 
function imaxloc (arr)
  use type
  real (dp), dimension (:), intent (in) :: arr
  integer (i4b) :: imaxloc
  integer (i4b), dimension (1) :: imax
  imax = maxloc (arr(:))
  imaxloc = imax (1)
end function imaxloc
 
function outerprod (a, b)
  use type
  real (dp), dimension (:), intent (in) :: a, b
  real (dp), dimension (size(a), size(b)) :: outerprod
  outerprod = spread (a, dim=2, ncopies=size(b)) * spread (b, dim=1, ncopies=size(a))
end function outerprod
 
subroutine ludcmp (a, indx, d)
  use type
  implicit none
  interface
    subroutine nrerror (string)
      character (len=*), intent (in) :: string
    end subroutine nrerror
 
    function assert_eq (n1, n2, n3, string)
      character (len=*), intent (in) :: string
      integer, intent (in) :: n1, n2, n3
      integer :: assert_eq
    end function assert_eq
 
    function imaxloc (arr)
      use type
      real (dp), dimension (:), intent (in) :: arr
      integer (i4b) :: imaxloc
    end function imaxloc
 
    function outerprod (a, b)
      use type
      real (dp), dimension (:), intent (in) :: a, b
      real (dp), dimension (size(a), size(b)) :: outerprod
    end function outerprod
 
    subroutine swap_d (a, b)
      use type
      real (dp), dimension (:), intent (inout) :: a, b
    end subroutine swap_d
  end interface
 
  real (dp), dimension (:, :), intent (inout) :: a
  integer (i4b), dimension (:), intent (out) :: indx
  real (dp), intent (out) :: d
  real (dp), dimension (size(a, 1)) :: vv
!  REAL(DP), PARAMETER :: TINY=1.0e-20
  integer (i4b) :: j, n, imax
  n = assert_eq (size(a, 1), size(a, 2), size(indx), 'ludcmp')
  d = 1.0d0
  vv = maxval (abs(a), dim=2)
  if (any(vv == 0.0)) call nrerror ('singular matrix in ludcmp')
  vv = 1.0d0 / vv
  do j = 1, n
    imax = (j-1) + imaxloc (vv(j:n)*abs(a(j:n, j)))
    if (j /= imax) then
      call swap_d (a(imax, :), a(j, :))
      d = - d
      vv (imax) = vv (j)
    end if
    indx (j) = imax
    if (a(j, j) == 0.0) a (j, j) = tiny
    a (j+1:n, j) = a (j+1:n, j) / a (j, j)
    a (j+1:n, j+1:n) = a (j+1:n, j+1:n) - outerprod (a(j+1:n, j), a(j, j+1:n))
  end do
end subroutine ludcmp
 
subroutine lubksb (a, indx, b)
  use type
  implicit none
  interface
    function assert_eq (n1, n2, n3, string)
      character (len=*), intent (in) :: string
      integer, intent (in) :: n1, n2, n3
      integer :: assert_eq
    end function assert_eq
  end interface
 
  real (dp), dimension (:, :), intent (in) :: a
  integer (i4b), dimension (:), intent (in) :: indx
  real (dp), dimension (:), intent (inout) :: b
  integer (i4b) :: i, n, ii, ll
  real (dp) :: summ
  n = assert_eq (size(a, 1), size(a, 2), size(indx), 'lubksb')
  ii = 0
  do i = 1, n
    ll = indx (i)
    summ = b (ll)
    b (ll) = b (i)
    if (ii /= 0) then
      summ = summ - dot_product (a(i, ii:i-1), b(ii:i-1))
    else if (summ /= 0.0) then
      ii = i
    end if
    b (i) = summ
  end do
  do i = n, 1, - 1
    b (i) = (b(i)-dot_product(a(i, i+1:n), b(i+1:n))) / a (i, i)
  end do
end subroutine lubksb
