module type
  integer, parameter :: i4b = selected_int_kind (9)
  integer, parameter :: sp = kind (1.0)
  integer, parameter :: dp = kind (1.0D0)
  real (dp), parameter :: pi = 3.141592653589793238462643383279502884197D0
  real (dp), parameter :: tiny = 0.00000000000000000000000000000001D0
end module type
