Module grid_lists
  Use nrtype
  Implicit None
  Save
 ! --------------------------------------------------------------------------------------
 ! This structure is used in gridsubset (and routines that call gridsubset) to
 ! provide a convenient structure to hold output information.
 ! --------------------------------------------------------------------------------------
 ! structure of indices in a 2d array
  Type INDX2D
    Integer (I4B) :: I, J ! indices into a 2d array
  End Type INDX2D
 ! structure for each link in the list
  Type GRDLST
    Type (INDX2D) :: INDX ! pair of indices
    Real (DP) :: DIST ! distance to grid point
    Type (GRDLST), Pointer :: NEXT ! pointer to next link in list
  End Type GRDLST
 ! structure for the main array holding information about number of links in a
 ! list as well as a pointer to the first link
  Type GRDARRAY
    Integer (I4B) :: LNKNR ! number of links in list
    Type (GRDLST), Pointer :: FIRST ! pointer to first link in list
  End Type GRDARRAY
 ! --------------------------------------------------------------------------------------
 ! structure used in disaggrhum.f90 holding linkage information between
 ! temperature data grids and relative humidity grids
  Type (GRDARRAY), Dimension (:, :), Pointer :: GL_RHTP
 ! --------------------------------------------------------------------------------------
End Module grid_lists
