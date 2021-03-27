module trig_degrees
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Einar Örn Hreinsson, Nov 2007 (einarhre@gmail.com)
  ! The family of functions sind, cosd and tand
  ! are not intrinsic to the gfortran compiler. These functions
  ! are implemented here in terms of the instrinsic functions sin,
  ! cos and tan.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use nrtype
  implicit none
  interface sind
    module procedure sind_r, sind_d
  end interface
  interface cosd
    module procedure cosd_r, cosd_d
  end interface
  interface tand
    module procedure tand_r, tand_d
  end interface
contains
  function sind_r (degree)
  ! Input argument to sin in degrees
    real (sp), intent (in) :: degree ! argument to sind in degrees
  ! Internal variable
    real (sp) :: crad ! conversion from degrees to radians
    real (sp) :: degree_rad ! argument to sin in radians
  ! Output
    real (sp) :: sind_r ! the value of sin function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    crad = pi / 180.0_sp
    degree_rad = degree * crad
   ! Value of SIND is equal to value of SIN when the argument is given
   ! in radians
    sind_r = sin (degree_rad)
  end function sind_r
  function sind_d (degree)
  ! Input argument to sin in degrees
    real (dp), intent (in) :: degree ! argument to sind in degrees
  ! Internal variable
    real (dp) :: crad ! conversion from degrees to radians
    real (dp) :: degree_rad ! argument to sin in radians
  ! Output
    real (dp) :: sind_d ! the value of sin function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    crad = pi_d / 180.0_dp
    degree_rad = degree * crad
   ! Value of SIND is equal to value of SIN when the argument is given
   ! in radians
    sind_d = sin (degree_rad)
  end function sind_d
  function cosd_r (degree)
  ! Input argument to cosd in degrees
    real (sp), intent (in) :: degree ! argument to cosd in degrees
  ! Internal variable
    real (sp) :: crad ! conversion from degrees to radians
    real (sp) :: degree_rad ! argument to cos in radians
  ! Output
    real (sp) :: cosd_r ! the value of cos function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    crad = pi / 180.0_sp
    degree_rad = degree * crad
   ! Value of COSD is equal to value of COS when the argument is given
   ! in radians
    cosd_r = cos (degree_rad)
  end function cosd_r
  function cosd_d (degree)
  ! Input argument to cosd in degrees
    real (dp), intent (in) :: degree ! argument to cosd in degrees
  ! Internal variable
    real (dp) :: crad ! conversion from degrees to radians
    real (dp) :: degree_rad ! argument to cos in radians
  ! Output
    real (dp) :: cosd_d ! the value of cos function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    crad = pi_d / 180.0_dp
    degree_rad = degree * crad
   ! Value of COSD is equal to value of COS when the argument is given
   ! in radians
    cosd_d = cos (degree_rad)
  end function cosd_d
  function tand_r (degree)
  ! Input argument to tand in degrees
    real (sp), intent (in) :: degree ! argument to tand in degrees
  ! Internal variable
    real (sp) :: crad ! conversion from degrees to radians
    real (sp) :: degree_rad ! argument to tan in radians
  ! Output
    real (sp) :: tand_r ! the value of tan function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    crad = pi / 180.0_sp
    degree_rad = degree * crad
   ! Value of tand is equal to value of tan when the argument is given
   ! in radians
    tand_r = tan (degree_rad)
  end function tand_r
  function tand_d (degree)
  ! Input argument to tand in degrees
    real (dp), intent (in) :: degree ! argument to tand in degrees
  ! Internal variable
    real (dp) :: crad ! conversion from degrees to radians
    real (dp) :: degree_rad ! argument to tan in radians
  ! Output
    real (dp) :: tand_d ! the value of tan function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    crad = pi_d / 180.0_dp
    degree_rad = degree * crad
   ! Value of tand is equal to value of tan when the argument is given
   ! in radians
    tand_d = tan (degree_rad)
  end function tand_d
end module trig_degrees
