!
! Solve linear equation for x (Ax = b => x = bA^-1) using LU decomposition and back substitution.
! Input:
!   X  = An m by n array.
!   TX = Precalculated transpose array of X, size n by m
!   Y  = An m-element vector containing the right-hand side of the linear system Ax = b.
! Output:
!   B  = An n-element vector.

subroutine least_squares (x, y, tx, b)
  use type
  implicit none
 
  interface
    subroutine ludcmp (a, indx, d)
      use type
      real (dp), dimension (:, :), intent (inout) :: a
      integer (i4b), dimension (:), intent (out) :: indx
      real (dp), intent (out) :: d
    end subroutine ludcmp
 
    subroutine lubksb (a, indx, b)
      use type
      real (dp), dimension (:, :), intent (in) :: a
      integer (i4b), dimension (:), intent (in) :: indx
      real (dp), dimension (:), intent (inout) :: b
    end subroutine lubksb
  end interface
 
  real (dp), intent (in) :: x (:, :)
  real (dp), intent (in) :: y (:)
  real (dp), intent (in) :: tx (:, :)
  real (dp), allocatable, intent (out) :: b (:)
 
  real (dp), allocatable :: a (:, :)
  integer (i4b), allocatable :: indx (:)
  integer (i4b) :: nvars, ntimes
  real (dp) :: d
 
  nvars = size (x, 2) - 1
  ntimes = size (y)
 
!  print *,'least'
  allocate (b(nvars+1))
  allocate (a(nvars+1, nvars+1))
  allocate (indx(nvars+1))
!  print *,'allocated'
  !print *, "Y = ", Y
  !print *, "X = ", X
  b = matmul (tx, y)
  a = matmul (tx, x)
  !print *, "A = ", A
  !print *, "B = ", B
!print *,'lu start'
  call ludcmp (a, indx, d)
!print *,'ludcmp'
  if (any(abs(a) < 9.99999968E-15)) then
    b (:) = 0.0d0
!AJN
!     print *, "Warning, LUdcmp produced a zero."
    return
  end if
!print *,'lubksb'
  !print *, "LU A = ", A
  call lubksb (a, indx, b)
 
  !print *, matmul(matmul(TX, X), B)
!print *,'deallocate'
  deallocate (a)
  deallocate (indx)
!print *,'done'
end subroutine least_squares
 
 
! -----------------------------
subroutine logistic_regression (x, y, tx, yp, b)
  use type
  implicit none
 
  interface
    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares
  end interface
 
  real (dp), intent (in) :: x (:, :)
  real (dp), intent (in) :: y (:)
  real (dp), intent (in) :: tx (:, :)
  integer (i4b), intent (in) :: yp (:)
  real (dp), allocatable, intent (out) :: b (:)
 
  real (dp), allocatable :: ypd (:), p (:), yn (:), bn (:), v (:, :), xv (:, :)
  !  real(DP), allocatable :: P(:), YN(:), BN(:), V(:,:), XV(:,:)
  integer (i4b) :: nvars, ntimes, i, t, f, it
  !real (dp) :: d
 
  nvars = size (x, 2) - 1
  ntimes = size (y)
 
  allocate (b(nvars+1))
  allocate (ypd(ntimes))
  allocate (yn(ntimes))
  allocate (p(ntimes))
  allocate (v(ntimes, ntimes))
  allocate (xv(ntimes, nvars+1))
 
  do t = 1, ntimes, 1
    if (yp(t) .gt. 0.0) then
      ypd (t) = 1.0d0
    else
      ypd (t) = 0.0d0
    end if
  end do
 
  b = 0.0d0
  i = 0
  it = 0
  f = 0 ! AWW added initialization
  do while (f /=  1)
    !print*, 'matmul(x,b):'
    !print*, matmul(x, b)
    if(any(-matmul(X,B) > 50)) then
      f = 1
    else
      p = 1.0d0 / (1.0d0+exp(-matmul(x, b)))
    endif
    if (any(p > 0.99)) then
      ! print *, "WARNING: logistic regression diverging"
      f = 1
    else
 
      yn = ypd - p
      v = 0.0d0
      do t = 1, ntimes, 1
        v (t, t) = p (t) * (1.0d0-p(t))
      end do
      xv = matmul (v, x)

      call least_squares (xv, yn, tx, bn)
 
      f = 1
      do i = 1, nvars + 1, 1
        if (bn(i) .gt. 1.0E-04 .or. bn(i) .lt.-1.0E-04) then
          f = 0
        end if
      end do
      if (it > 10) then
        ! print *, "WARNING: logistic regression failed to converge"
        f = 1
      end if
 
      b = b + bn
      deallocate (bn)
 
    end if
    it = it + 1
  end do
  ! print *, "Iterations = ", it
  ! print *, "Final B = ", B
 
end subroutine logistic_regression
 
 
subroutine logistic_regressionrf (x, y, tx, b)
  use type
  implicit none
 
  interface
    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares
  end interface
 
  real (dp), intent (in) :: x (:, :)
  real (dp), intent (in) :: y (:)
  real (dp), intent (in) :: tx (:, :)
  real (dp), allocatable, intent (out) :: b (:)
 
  real (dp), allocatable :: ypd (:), p (:), yn (:), bn (:), v (:, :), xv (:, :)
  !  real(DP), allocatable :: P(:), YN(:), BN(:), V(:,:), XV(:,:)
  integer (i4b) :: nvars, ntimes, i, t, f, it
  !real (dp) :: d
 
  nvars = size (x, 2) - 1
  ntimes = size (y)
 
  !print *, 'logistic regression',ntimes,nvars,Yp
 
  allocate (b(nvars+1))
  allocate (ypd(ntimes))
  allocate (yn(ntimes))
  allocate (p(ntimes))
  allocate (v(ntimes, ntimes))
  allocate (xv(ntimes, nvars+1))
 
  do t = 1, ntimes, 1
    if (y(t) .ne. 0.0) then
      ypd (t) = 1.0d0
    else
      ypd (t) = 0.0d0
    end if
  end do
 
  b = 0.0d0
  i = 0
  it = 0
  f = 1 ! AWW added initialization
  !print *, "B = ", B
  do while (f /=  1)
    p = 1.0d0 / (1.0d0+exp(-matmul(x, b)))
    if (any(p > 0.97)) then
        !print *, "WARNING: logistic regression diverging"
      f = 1
    else
 
      yn = ypd - p
      v = 0.0d0
      do t = 1, ntimes, 1
        v (t, t) = p (t) * (1.0d0-p(t))
      end do
      xv = matmul (v, x)
      call least_squares (xv, yn, tx, bn)
 
      f = 1
      do i = 1, nvars + 1, 1
        if (bn(i) .gt. 1.0E-04 .or. bn(i) .lt.-1.0E-04) then
          f = 0
        end if
      end do
      if (it > 8) then
        ! print *, "WARNING: logistic regression failed to converge"
        f = 1
      end if
 
      b = b + bn
      deallocate (bn)
 
    end if
    it = it + 1
  end do
  ! print *, "Iterations = ", it
  ! print *, "Final B = ", B
 
end subroutine logistic_regressionrf
