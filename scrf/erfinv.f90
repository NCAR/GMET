function erfinv(x)
  use nrtype
  implicit none

  real(sp),intent(in)    :: x

! real(sp),dimension(4) :: a = (/0.886226899, -1.645349621, 0.914624893, -0.140543331/)
! real(dp),dimension(4) :: b = (/-2.118377725, 1.442710462, -0.329097515, 0.0122229801/)
! real(dp),dimension(4) :: c = (/-1.970840454, -1.624906493, 3.429567803, 1.641345311/)
! real(dp),dimension(2) :: d = (/3.543889200, 1.637067800/)
!  real(dp)              :: y_o

  real(sp)               :: erfinv
  real(sp)               :: tmp
  integer                :: neg
  
  neg = 0
!  print *,'erfinv ',neg,x
  if(x .lt. 0.0) then
    neg = 1
    erfinv = -x
  else
    erfinv = x
  endif

  if(erfinv .le. 0.7) then
    tmp = erfinv*erfinv
    erfinv = erfinv * (((-0.140543331*tmp+0.914624893)*tmp-1.645349621)*tmp+0.886226899)/((((0.012229801*tmp-0.329097515)*tmp+1.442710462)*tmp-2.118377725)*tmp+1.0)
  else
    tmp = sqrt(-log(0.5*(1.0-erfinv)))
    erfinv = (((1.641345311*tmp+3.429567803)*tmp-1.624906493)*tmp-1.970840454)/((1.637067800*tmp+3.543889200)*tmp+1.0)
  endif

!  print *, 'erfinv ',tmp,erfinv,neg
  if(neg) then
    erfinv = -erfinv
  endif

  return
end function erfinv
