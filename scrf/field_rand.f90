SUBROUTINE FIELD_RAND (NSPL1, NSPL2, CFIELD)
! ----------------------------------------------------------------------------------------
! Purpose:
!
!   Used to estimate a correlated field of random numbers for the basin gridpoints
!
! ----------------------------------------------------------------------------------------
!
! Modified:
!	Andy Newman: Aug 2013
!
! ----------------------------------------------------------------------------------------
! I/O:
!
!   Input(s):
!   ---------
!    NSPL1: number of location points (1st spatial dimension)
!    NSPL2: number of location points (2nd spatial dimension)
!
! ----------------------------------------------------------------------------------------
! Structures Used:
!
!   (1) INPUTDAT2D
!   --------------
!   Correspondence between basins and gridpoints
 
!   (2) GRIDWEIGHT
!   --------------
!   Weights applied to each grid point to generate spatially correlated random fields
!
! ----------------------------------------------------------------------------------------
 
  USE nrtype ! variable types (DP, I4B, etc.)
  USE nr, ONLY: gasdev ! Num. Recipies
  USE inputdat2d ! use to relate basins to gridpoints
 
  USE gridweight ! correlation structure
  IMPLICIT NONE
 
! input
  INTEGER (I4B), INTENT (IN) :: NSPL1 ! # points (1st spatial dimension)
  INTEGER (I4B), INTENT (IN) :: NSPL2 ! # points (2nd spatial dimension)
 
  ! Output
  REAL (DP), DIMENSION (NSPL1, NSPL2) :: CFIELD ! correlated random field
 
 
  ! Internal variables
  INTEGER (I4B) :: NLON ! number of x points in the 2-d grid
  INTEGER (I4B) :: NLAT ! number of y points in the 2-d grid
  INTEGER (I4B) :: IERR ! error code for alolocate statement
  INTEGER (I4B) :: IGRD ! loop through the gridpoints
  INTEGER (I4B) :: ILON, ILAT ! (i,j) position of the igrd-th point
  INTEGER (I4B), POINTER :: JLON, JLAT ! (i,j) position of prev generated point
  INTEGER (I4B) :: IPREV ! loop thru previously generated points
  INTEGER (I4B) :: NPRV ! number of previously generated points
  INTEGER (I4B) :: ILNK ! index to the nearest grid point
  REAL (DP), DIMENSION (:), ALLOCATABLE :: VPRV ! vector of previously generated points
  REAL (DP) :: XBAR ! conditional mean at each grid
  REAL (SP) :: ARAN ! a single random number
  REAL (DP), DIMENSION (:, :), ALLOCATABLE :: CRAN ! grid of correlated random numbers
  INTEGER (I4B) :: IRCH ! loop through the stream segments
 
  ! ----------------------------------------------------------------------------------------
  ! (1) GET THE NUMBER OF X AND Y POINTS AND ALLOCATE SPACE FOR THE RANDOM GRID
  ! ----------------------------------------------------------------------------------------
  NLON = SIZE (SPCORR, 1)
  NLAT = SIZE (SPCORR, 2)
  ALLOCATE (CRAN(NLON, NLAT), STAT=IERR)
  IF (IERR .NE. 0) CALL EXIT_SCRF (1, 'CRAN: PROBLEM ALLOCATING SPACE')
  ! ----------------------------------------------------------------------------------------
  ! (1) GENERATE GRID OF RANDOM NUMBERS
  ! ----------------------------------------------------------------------------------------
  ! loop through the grid points
  DO IGRD = 1, NLON * NLAT
  ! identify the (i,j) position of the igrd-th point
   ILON = IORDER (IGRD)
   ILAT = JORDER (IGRD)
  ! assign a random number to the first grid-point
   IF (IGRD .EQ. 1) THEN
    CALL gasdev (ARAN)
    CRAN (ILON, ILAT) = ARAN
  ! process gridpoints 2,...,n
   ELSE
    ! get the number of "previously generated points"
    NPRV = SIZE (SPCORR(ILON, ILAT)%WGHT)
    ALLOCATE (VPRV(NPRV), STAT=IERR)
    IF (IERR .NE. 0) CALL EXIT_SCRF (1, 'problem allocating array [field_rand.f90]')
    ! build a vector of previously generated points
    DO IPREV = 1, NPRV
     JLON => SPCORR(ILON, ILAT)%IPOS(IPREV)! i-position of previously generated point
     JLAT => SPCORR(ILON, ILAT)%JPOS(IPREV)! j-position of previously generated point
     VPRV (IPREV) = CRAN (JLON, JLAT)! (previously generated point)
    END DO ! iprev
    ! and generate the "current" point
    CALL gasdev (ARAN)
    XBAR = DOT_PRODUCT (VPRV(1:NPRV), SPCORR(ILON, ILAT)%WGHT(1:NPRV))
    CRAN (ILON, ILAT) = XBAR + SPCORR(ILON, ILAT)%SDEV * ARAN
    ! free up VPRV so that we can use it again for the next grid
    DEALLOCATE (VPRV, STAT=IERR)
    IF (IERR .NE. 0) CALL EXIT_SCRF (1, 'problem deallocating array [field_rand.f90]')
   END IF ! (if not the first point)
  END DO ! igrd
 
  CFIELD = CRAN
 
  RETURN
 
END SUBROUTINE FIELD_RAND
