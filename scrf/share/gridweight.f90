module gridweight
  use nrtype
  implicit none
  save
  type splnum
    integer (i4b), dimension (:), pointer :: ipos ! i position of previously generated points
    integer (i4b), dimension (:), pointer :: jpos ! j position of previously generated points
    real (dp), dimension (:), pointer :: wght ! weights for previously generated points
    real (dp) :: sdev ! standard deviation of the estimate
  end type splnum
  type (splnum), dimension (:, :), pointer :: spcorr ! spatially-correlated random numbers
  integer (i4b), dimension (:), pointer :: iorder ! i-position, in processing order
  integer (i4b), dimension (:), pointer :: jorder ! j-position, in processing order
end module gridweight
