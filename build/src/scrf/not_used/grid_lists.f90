module grid_lists
  use nrtype
  implicit none
  save
 ! --------------------------------------------------------------------------------------
 ! This structure is used in gridsubset (and routines that call gridsubset) to
 ! provide a convenient structure to hold output information.
 ! --------------------------------------------------------------------------------------
 ! structure of indices in a 2d array
  type indx2d
    integer (i4b) :: i, j ! indices into a 2d array
  end type indx2d
 ! structure for each link in the list
  type grdlst
    type (indx2d) :: indx ! pair of indices
    real (dp) :: dist ! distance to grid point
    type (grdlst), pointer :: next ! pointer to next link in list
  end type grdlst
 ! structure for the main array holding information about number of links in a
 ! list as well as a pointer to the first link
  type grdarray
    integer (i4b) :: lnknr ! number of links in list
    type (grdlst), pointer :: first ! pointer to first link in list
  end type grdarray
 ! --------------------------------------------------------------------------------------
 ! structure used in disaggrhum.f90 holding linkage information between
 ! temperature data grids and relative humidity grids
  type (grdarray), dimension (:, :), pointer :: gl_rhtp
 ! --------------------------------------------------------------------------------------
end module grid_lists
