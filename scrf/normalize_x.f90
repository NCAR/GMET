subroutine normalize_x(x,mean,stdev)
  use nrtype
  implicit none

  real(DP), intent(in) :: x(:,:)
  real(DP), intent(out) :: mean
  real(DP), intent(out) :: stdev

!  real(DP) :: mean, stdev, sum_x2, sum_x
  real(DP) :: sum_x2, sum_x
  integer(I4B) :: i,j
  integer(I4B) :: nx,ny

  nx = size(x,1)
  ny = size(x,2)

  sum_x = 0.0d0
  sum_x2 = 0.0d0
  do i = 1, nx, 1
    do j = 1, ny, 1
      sum_x  = sum_x + x(i,j)
      sum_x2 = sum_x2 + x(i,j)**2 
    enddo
  enddo

  mean = sum_x / real(nx*ny)
  stdev = sqrt(  (real(nx*ny) * sum_x2 - sum_x**2) /  (real(nx*ny)*real(nx*ny-1)) )

!  do i = 1, nx, 1
!    do j = 1, ny, 1 
!      if(stdev .eq. 0.0) then
!	x(i,j) = x(i,j)
!      else
!	x(i,j) = (x(i,j) - mean) / stdev
!      endif
!    enddo
!  enddo

end subroutine normalize_x