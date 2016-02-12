Module gridweight
  Use nrtype
  Implicit None
  Save
  Type SPLNUM
    Integer (I4B), Dimension (:), Pointer :: IPOS ! i position of previously generated points
    Integer (I4B), Dimension (:), Pointer :: JPOS ! j position of previously generated points
    Real (DP), Dimension (:), Pointer :: WGHT ! weights for previously generated points
    Real (DP) :: SDEV ! standard deviation of the estimate
  End Type SPLNUM
  Type (SPLNUM), Dimension (:, :), Pointer :: SPCORR ! spatially-correlated random numbers
  Integer (I4B), Dimension (:), Pointer :: IORDER ! i-position, in processing order
  Integer (I4B), Dimension (:), Pointer :: JORDER ! j-position, in processing order
End Module gridweight
