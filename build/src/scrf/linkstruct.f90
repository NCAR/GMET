module linkstruct
!
!	Modified:
!	Andy Newman, Aug 2013
!
  use nrtype
  implicit none
  save
 ! --------------------------------------------------------------------------------------
 ! structure to hold indices of desired data in the input files
  type des_ix
    integer (i4b) :: spl1_start ! start index (1st spatial dimension)
    integer (i4b) :: spl1_count ! start index (1st spatial dimension)
    integer (i4b) :: spl2_start ! count (2nd spatial dimension)
    integer (i4b) :: spl2_count ! count (2nd spatial dimension)
  end type des_ix
 ! --------------------------------------------------------------------------------------
 ! coordinates of the "current" grid
  type coords
    type (des_ix) :: idx ! indices of data in input files
    real (dp), dimension (:, :), pointer :: lat ! latitude
    real (dp), dimension (:, :), pointer :: lon ! longitude
    real (dp), dimension (:, :), pointer :: elv ! elevation
  end type coords
 ! --------------------------------------------------------------------------------------
 ! (i,j) and (lat,lon) of individual grid points that are assigned to a given sub-basin
  type basgrd
    integer (i4b) :: ipos ! i-position
    integer (i4b) :: jpos ! j-position
    real (dp) :: wght ! weight assigned to grid
    real (dp) :: crat ! mean annual precip (bas/sta)
    real (dp) :: corr ! correction (precip/(1+corr))
    real (dp) :: alat ! latitude of i-j position
    real (dp) :: alon ! longitude of i-j position
  end type basgrd
 ! --------------------------------------------------------------------------------------
 ! collection of grid points that are assigned to a given sub-basin
  type interp
    type (basgrd), dimension (:), pointer :: link
  end type interp
 ! --------------------------------------------------------------------------------------
 ! structure to hold the grid-to-basin information
  type (interp), dimension (:), pointer :: grid_close ! save temporary linkage information
 ! --------------------------------------------------------------------------------------
 ! structure to hold the current grid
  type (coords), pointer :: holding ! "current" grid of lat-lon
 ! --------------------------------------------------------------------------------------
end module linkstruct
