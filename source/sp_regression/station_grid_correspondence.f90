subroutine station_grid_correspondence(X,Z,close_weights,close_weights_t,close_loc,close_loc_t,nSta,nearestGridpoint)
  use type
  implicit none

  real(dp), intent(in)        :: X(:,:)
  real(dp), intent(in)        :: Z(:,:)
  real(dp), intent(in)        :: close_weights(:,:)
  real(dp), intent(in)        :: close_weights_t(:,:)
  integer(I4B), intent(in)    :: close_loc(:,:) 
  integer(I4B), intent(in)    :: close_loc_t(:,:) 
  integer(I4B), intent(in)    :: nSta
  integer(I4B), intent(out)   :: nearestGridpoint(:)  

  ! local variables 
  integer(I4B)                :: g, s                  ! counter variables
  integer(I4B)                :: vshape(2)             ! array shape 
  integer(I4B)                :: nStaGrid
  integer(I4B)                :: nGrid
  real(dp),allocatable        :: closestWeight(:)      ! weights for closest gridpoint to a station
  integer(I4B)                :: tmp_grid_pt
  real(dp)                    :: tmp_weight
  real(dp)                    :: max_tmp_weight

  ! ------- code starts below ----------
  
  ! Associate each grid point with the nearest station and station weight
  !   AW: to-do:  this can be vectorized with advanced indexing
  
  vshape = shape(close_weights)

  allocate(closestWeight(nSta))

  ! initialize
  closestWeight    = 0.0
  max_tmp_weight   = 0.0
  nearestGridpoint = -1
  
  ngrid    = vshape(1)
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

      ! catch errors
      if(nearestGridPoint(close_loc(g,s)) .lt. 1) then
        print *,'p cell',g,'closest weights not found', close_weights(g,:)
      end if
      if(nearestGridPoint(close_loc_t(g,s)) .lt. 1) then
        print *,'t cell',g,'closest weights not found', close_weights_t(g,:)
      end if
      
    end do
  end do

  !check for stations that do not have a grid correspondence yet
  do s = 1, nSta, 1
    if(nearestGridPoint(s) .eq. 0) then
      do g = 1, ngrid, 1
        ! x() are station lonlat; z() are grid lonlat; returns weight (w_base) for grd-to-stn
        call calc_distance_weight(real(500.,dp), X(s,2), X(s,3), Z(g,2), Z(g,3), tmp_weight)      
        if(tmp_weight > max_tmp_weight) then
          max_tmp_weight = tmp_weight
          tmp_grid_pt = g
        end if
      end do
      nearestGridPoint(s) = tmp_grid_pt
      !reset max_tmp_weight
      max_tmp_weight = 0
    end if
  end do

!  print *,nearestGridPoint
end subroutine
