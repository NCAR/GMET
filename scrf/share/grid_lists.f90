MODULE grid_lists
 USE nrtype
 IMPLICIT NONE
 SAVE
 ! --------------------------------------------------------------------------------------
 ! This structure is used in gridsubset (and routines that call gridsubset) to
 ! provide a convenient structure to hold output information.
 ! --------------------------------------------------------------------------------------
 ! structure of indices in a 2d array
 TYPE INDX2D
   INTEGER(I4B)                          :: I,J         ! indices into a 2d array
 END TYPE INDX2D
 ! structure for each link in the list
 TYPE GRDLST
  TYPE(INDX2D)                           :: INDX        ! pair of indices
  REAL(DP)                               :: DIST        ! distance to grid point
  TYPE(GRDLST),POINTER                   :: NEXT        ! pointer to next link in list
 END TYPE GRDLST
 ! structure for the main array holding information about number of links in a
 ! list as well as a pointer to the first link
 TYPE GRDARRAY
  INTEGER(I4B)                           :: LNKNR       ! number of links in list
  TYPE(GRDLST),POINTER                   :: FIRST       ! pointer to first link in list
 END TYPE GRDARRAY
 ! --------------------------------------------------------------------------------------
 ! structure used in disaggrhum.f90 holding linkage information between
 ! temperature data grids and relative humidity grids
 TYPE(GRDARRAY),DIMENSION(:,:),POINTER   :: GL_RHTP
 ! --------------------------------------------------------------------------------------
END MODULE grid_lists
