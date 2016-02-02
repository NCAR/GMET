MODULE gridweight
 USE nrtype
 IMPLICIT NONE
 SAVE
 TYPE SPLNUM
  INTEGER(I4B), DIMENSION(:), POINTER      :: IPOS     ! i position of previously generated points
  INTEGER(I4B), DIMENSION(:), POINTER      :: JPOS     ! j position of previously generated points
  REAL(DP), DIMENSION(:), POINTER          :: WGHT     ! weights for previously generated points
  REAL(DP)                                 :: SDEV     ! standard deviation of the estimate
 ENDTYPE SPLNUM
 TYPE(SPLNUM), DIMENSION(:,:), POINTER     :: SPCORR   ! spatially-correlated random numbers
 INTEGER(I4B), DIMENSION(:), POINTER       :: IORDER   ! i-position, in processing order
 INTEGER(I4B), DIMENSION(:), POINTER       :: JORDER   ! j-position, in processing order
END MODULE gridweight
