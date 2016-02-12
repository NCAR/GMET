Module linkstruct
!
!	Modified:
!	Andy Newman, Aug 2013
!
  Use nrtype
  Implicit None
  Save
 ! --------------------------------------------------------------------------------------
 ! structure to hold indices of desired data in the input files
  Type DES_IX
    Integer (I4B) :: SPL1_START ! start index (1st spatial dimension)
    Integer (I4B) :: SPL1_COUNT ! start index (1st spatial dimension)
    Integer (I4B) :: SPL2_START ! count (2nd spatial dimension)
    Integer (I4B) :: SPL2_COUNT ! count (2nd spatial dimension)
  End Type DES_IX
 ! --------------------------------------------------------------------------------------
 ! coordinates of the "current" grid
  Type COORDS
    Type (DES_IX) :: IDX ! indices of data in input files
    Real (DP), Dimension (:, :), Pointer :: LAT ! latitude
    Real (DP), Dimension (:, :), Pointer :: LON ! longitude
    Real (DP), Dimension (:, :), Pointer :: ELV ! elevation
  End Type COORDS
 ! --------------------------------------------------------------------------------------
 ! (i,j) and (lat,lon) of individual grid points that are assigned to a given sub-basin
  Type BASGRD
    Integer (I4B) :: IPOS ! i-position
    Integer (I4B) :: JPOS ! j-position
    Real (DP) :: WGHT ! weight assigned to grid
    Real (DP) :: CRAT ! mean annual precip (bas/sta)
    Real (DP) :: CORR ! correction (precip/(1+corr))
    Real (DP) :: ALAT ! latitude of i-j position
    Real (DP) :: ALON ! longitude of i-j position
  End Type BASGRD
 ! --------------------------------------------------------------------------------------
 ! collection of grid points that are assigned to a given sub-basin
  Type INTERP
    Type (BASGRD), Dimension (:), Pointer :: LINK
  End Type INTERP
 ! --------------------------------------------------------------------------------------
 ! structure to hold the grid-to-basin information
  Type (INTERP), Dimension (:), Pointer :: GRID_CLOSE ! save temporary linkage information
 ! --------------------------------------------------------------------------------------
 ! structure to hold the current grid
  Type (COORDS), Pointer :: HOLDING ! "current" grid of lat-lon
 ! --------------------------------------------------------------------------------------
End Module linkstruct
