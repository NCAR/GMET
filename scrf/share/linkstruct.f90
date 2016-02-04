MODULE linkstruct
!
!	Modified:
!	Andy Newman, Aug 2013
!
  USE nrtype
  IMPLICIT NONE
  SAVE
 ! --------------------------------------------------------------------------------------
 ! structure to hold indices of desired data in the input files
  TYPE DES_IX
   INTEGER (I4B) :: SPL1_START ! start index (1st spatial dimension)
   INTEGER (I4B) :: SPL1_COUNT ! start index (1st spatial dimension)
   INTEGER (I4B) :: SPL2_START ! count (2nd spatial dimension)
   INTEGER (I4B) :: SPL2_COUNT ! count (2nd spatial dimension)
  END TYPE DES_IX
 ! --------------------------------------------------------------------------------------
 ! coordinates of the "current" grid
  TYPE COORDS
   TYPE (DES_IX) :: IDX ! indices of data in input files
   REAL (DP), DIMENSION (:, :), POINTER :: LAT ! latitude
   REAL (DP), DIMENSION (:, :), POINTER :: LON ! longitude
   REAL (DP), DIMENSION (:, :), POINTER :: ELV ! elevation
  END TYPE COORDS
 ! --------------------------------------------------------------------------------------
 ! (i,j) and (lat,lon) of individual grid points that are assigned to a given sub-basin
  TYPE BASGRD
   INTEGER (I4B) :: IPOS ! i-position
   INTEGER (I4B) :: JPOS ! j-position
   REAL (DP) :: WGHT ! weight assigned to grid
   REAL (DP) :: CRAT ! mean annual precip (bas/sta)
   REAL (DP) :: CORR ! correction (precip/(1+corr))
   REAL (DP) :: ALAT ! latitude of i-j position
   REAL (DP) :: ALON ! longitude of i-j position
  END TYPE BASGRD
 ! --------------------------------------------------------------------------------------
 ! collection of grid points that are assigned to a given sub-basin
  TYPE INTERP
   TYPE (BASGRD), DIMENSION (:), POINTER :: LINK
  END TYPE INTERP
 ! --------------------------------------------------------------------------------------
 ! structure to hold the grid-to-basin information
  TYPE (INTERP), DIMENSION (:), POINTER :: GRID_CLOSE ! save temporary linkage information
 ! --------------------------------------------------------------------------------------
 ! structure to hold the current grid
  TYPE (COORDS), POINTER :: HOLDING ! "current" grid of lat-lon
 ! --------------------------------------------------------------------------------------
END MODULE linkstruct
