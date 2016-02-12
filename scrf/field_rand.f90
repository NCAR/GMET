Subroutine FIELD_RAND (NSPL1, NSPL2, CFIELD)
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
 
  Use nrtype ! variable types (DP, I4B, etc.)
  Use nr, Only: gasdev ! Num. Recipies
  Use inputdat2d ! use to relate basins to gridpoints
 
  Use gridweight ! correlation structure
  Implicit None
 
! input
  Integer (I4B), Intent (In) :: NSPL1 ! # points (1st spatial dimension)
  Integer (I4B), Intent (In) :: NSPL2 ! # points (2nd spatial dimension)
 
  ! Output
  Real (DP), Dimension (NSPL1, NSPL2) :: CFIELD ! correlated random field
 
 
  ! Internal variables
  Integer (I4B) :: NLON ! number of x points in the 2-d grid
  Integer (I4B) :: NLAT ! number of y points in the 2-d grid
  Integer (I4B) :: IERR ! error code for alolocate statement
  Integer (I4B) :: IGRD ! loop through the gridpoints
  Integer (I4B) :: ILON, ILAT ! (i,j) position of the igrd-th point
  Integer (I4B), Pointer :: JLON, JLAT ! (i,j) position of prev generated point
  Integer (I4B) :: IPREV ! loop thru previously generated points
  Integer (I4B) :: NPRV ! number of previously generated points
  Integer (I4B) :: ILNK ! index to the nearest grid point
  Real (DP), Dimension (:), Allocatable :: VPRV ! vector of previously generated points
  Real (DP) :: XBAR ! conditional mean at each grid
  Real (SP) :: ARAN ! a single random number
  Real (DP), Dimension (:, :), Allocatable :: CRAN ! grid of correlated random numbers
  Integer (I4B) :: IRCH ! loop through the stream segments
 
  ! ----------------------------------------------------------------------------------------
  ! (1) GET THE NUMBER OF X AND Y POINTS AND ALLOCATE SPACE FOR THE RANDOM GRID
  ! ----------------------------------------------------------------------------------------
  NLON = SIZE (SPCORR, 1)
  NLAT = SIZE (SPCORR, 2)
  Allocate (CRAN(NLON, NLAT), Stat=IERR)
  If (IERR .Ne. 0) Call EXIT_SCRF (1, 'CRAN: PROBLEM ALLOCATING SPACE')
  ! ----------------------------------------------------------------------------------------
  ! (1) GENERATE GRID OF RANDOM NUMBERS
  ! ----------------------------------------------------------------------------------------
  ! loop through the grid points
  Do IGRD = 1, NLON * NLAT
  ! identify the (i,j) position of the igrd-th point
    ILON = IORDER (IGRD)
    ILAT = JORDER (IGRD)
  ! assign a random number to the first grid-point
    If (IGRD .Eq. 1) Then
      Call gasdev (ARAN)
      CRAN (ILON, ILAT) = ARAN
  ! process gridpoints 2,...,n
    Else
    ! get the number of "previously generated points"
      NPRV = SIZE (SPCORR(ILON, ILAT)%WGHT)
      Allocate (VPRV(NPRV), Stat=IERR)
      If (IERR .Ne. 0) Call EXIT_SCRF (1, 'problem allocating array [field_rand.f90]')
    ! build a vector of previously generated points
      Do IPREV = 1, NPRV
        JLON => SPCORR(ILON, ILAT)%IPOS(IPREV)! i-position of previously generated point
        JLAT => SPCORR(ILON, ILAT)%JPOS(IPREV)! j-position of previously generated point
        VPRV (IPREV) = CRAN (JLON, JLAT)! (previously generated point)
      End Do ! iprev
    ! and generate the "current" point
      Call gasdev (ARAN)
      XBAR = DOT_PRODUCT (VPRV(1:NPRV), SPCORR(ILON, ILAT)%WGHT(1:NPRV))
      CRAN (ILON, ILAT) = XBAR + SPCORR(ILON, ILAT)%SDEV * ARAN
    ! free up VPRV so that we can use it again for the next grid
      Deallocate (VPRV, Stat=IERR)
      If (IERR .Ne. 0) Call EXIT_SCRF (1, 'problem deallocating array [field_rand.f90]')
    End If ! (if not the first point)
  End Do ! igrd
 
  CFIELD = CRAN
 
  Return
 
End Subroutine FIELD_RAND
