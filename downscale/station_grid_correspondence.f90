subroutine station_grid_correspondence(close_weights,close_weights_t,close_loc,close_loc_t,nSta,nearestGridpoint)
  use type
  implicit none

  real(dp), intent(in)        :: close_weights(:,:)
  real(dp), intent(in)        :: close_weights_t(:,:)
  integer(I4B), intent(in)    :: close_loc(:,:) 
  integer(I4B), intent(in)    :: close_loc_t(:,:) 
  integer(I4B), intent(in)    :: nSta
  
  integer(I4B), intent(out)   :: nearestGridpoint(:)  

  !local variables 
  integer(I4B)                :: g, s                  ! counter variables
  integer(I4B)                :: vshape(2)             ! array shape 
  integer(I4B)                :: nStaGrid
  integer(I4B)                :: nGrid
  real(dp),allocatable        :: closestWeight(:)      ! weights for closest gridpoint to a station

  !code starts below
  vshape = shape(close_weights)


  allocate(closestWeight(nSta))

  closestWeight(:) = 0.0


  ngrid = vshape(1)
  nStaGrid = vshape(2)

  do g = 1, ngrid, 1
    do s = 1, nStaGrid, 1
      if( close_weights(g,s) >= closestWeight(close_loc(g,s)) ) then
        closestWeight(close_loc(g,s)) = close_weights(g,s)
        nearestGridPoint(close_loc(g,s)) = g;
      end if
      if( close_weights_t(g,s) >= closestWeight(close_loc_t(g,s)) ) then
        closestWeight(close_loc_t(g,s)) = close_weights_t(g,s)
        nearestGridPoint(close_loc_t(g,s)) = g;
      end if

    end do
  end do

end subroutine
