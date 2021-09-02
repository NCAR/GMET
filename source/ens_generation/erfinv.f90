function erfinv (x)
  use nrtype
  implicit none
!
  real (sp), intent (in) :: x
  real (sp) :: erfinv
  real (sp) :: tmp
  integer :: neg
!
  neg = 0
!
  if (x .lt. 0.0) then
    neg = 1
    erfinv = - x
  else
    erfinv = x
  end if
!
  if (erfinv .le. 0.7) then
    tmp = erfinv * erfinv
    erfinv = erfinv * (((-0.140543331*tmp+0.914624893)*tmp-1.645349621)*tmp+0.886226899) / &
   & ((((0.012229801*tmp-0.329097515)*tmp+1.442710462)*tmp-2.118377725)*tmp+1.0)
  else
    tmp = sqrt (-log(0.5*(1.0-erfinv)))
    erfinv = (((1.641345311*tmp+3.429567803)*tmp-1.624906493)*tmp-1.970840454) / &
   & ((1.637067800*tmp+3.543889200)*tmp+1.0)
  end if
!
  !if (neg) then   ! AW
  if (neg .eq. 1) then

    erfinv = - erfinv
  end if
!
  return
end function erfinv
