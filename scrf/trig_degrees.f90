Module trig_degrees
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Einar Örn Hreinsson, Nov 2007 (einarhre@gmail.com)
  ! The family of functions sind, cosd and tand
  ! are not intrinsic to the gfortran compiler. These functions
  ! are implemented here in terms of the instrinsic functions sin,
  ! cos and tan.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Use nrtype
  Implicit None
  Interface SIND
    Module Procedure SIND_R, SIND_D
  End Interface
  Interface COSD
    Module Procedure COSD_R, COSD_D
  End Interface
  Interface TAND
    Module Procedure TAND_R, TAND_D
  End Interface
Contains
  Function SIND_R (DEGREE)
  ! Input argument to sin in degrees
    Real (SP), Intent (In) :: DEGREE ! argument to sind in degrees
  ! Internal variable
    Real (SP) :: CRAD ! conversion from degrees to radians
    Real (SP) :: DEGREE_RAD ! argument to sin in radians
  ! Output
    Real (SP) :: SIND_R ! the value of sin function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    CRAD = PI / 180.0_SP
    DEGREE_RAD = DEGREE * CRAD
   ! Value of SIND is equal to value of SIN when the argument is given
   ! in radians
    SIND_R = Sin (DEGREE_RAD)
  End Function SIND_R
  Function SIND_D (DEGREE)
  ! Input argument to sin in degrees
    Real (DP), Intent (In) :: DEGREE ! argument to sind in degrees
  ! Internal variable
    Real (DP) :: CRAD ! conversion from degrees to radians
    Real (DP) :: DEGREE_RAD ! argument to sin in radians
  ! Output
    Real (DP) :: SIND_D ! the value of sin function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    CRAD = PI_D / 180.0_DP
    DEGREE_RAD = DEGREE * CRAD
   ! Value of SIND is equal to value of SIN when the argument is given
   ! in radians
    SIND_D = Sin (DEGREE_RAD)
  End Function SIND_D
  Function COSD_R (DEGREE)
  ! Input argument to cosd in degrees
    Real (SP), Intent (In) :: DEGREE ! argument to cosd in degrees
  ! Internal variable
    Real (SP) :: CRAD ! conversion from degrees to radians
    Real (SP) :: DEGREE_RAD ! argument to cos in radians
  ! Output
    Real (SP) :: COSD_R ! the value of cos function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    CRAD = PI / 180.0_SP
    DEGREE_RAD = DEGREE * CRAD
   ! Value of COSD is equal to value of COS when the argument is given
   ! in radians
    COSD_R = Cos (DEGREE_RAD)
  End Function COSD_R
  Function COSD_D (DEGREE)
  ! Input argument to cosd in degrees
    Real (DP), Intent (In) :: DEGREE ! argument to cosd in degrees
  ! Internal variable
    Real (DP) :: CRAD ! conversion from degrees to radians
    Real (DP) :: DEGREE_RAD ! argument to cos in radians
  ! Output
    Real (DP) :: COSD_D ! the value of cos function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    CRAD = PI_D / 180.0_DP
    DEGREE_RAD = DEGREE * CRAD
   ! Value of COSD is equal to value of COS when the argument is given
   ! in radians
    COSD_D = Cos (DEGREE_RAD)
  End Function COSD_D
  Function TAND_R (DEGREE)
  ! Input argument to tand in degrees
    Real (SP), Intent (In) :: DEGREE ! argument to tand in degrees
  ! Internal variable
    Real (SP) :: CRAD ! conversion from degrees to radians
    Real (SP) :: DEGREE_RAD ! argument to tan in radians
  ! Output
    Real (SP) :: TAND_R ! the value of tan function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    CRAD = PI / 180.0_SP
    DEGREE_RAD = DEGREE * CRAD
   ! Value of tand is equal to value of tan when the argument is given
   ! in radians
    TAND_R = Tan (DEGREE_RAD)
  End Function TAND_R
  Function TAND_D (DEGREE)
  ! Input argument to tand in degrees
    Real (DP), Intent (In) :: DEGREE ! argument to tand in degrees
  ! Internal variable
    Real (DP) :: CRAD ! conversion from degrees to radians
    Real (DP) :: DEGREE_RAD ! argument to tan in radians
  ! Output
    Real (DP) :: TAND_D ! the value of tan function
   ! Convert degrees to radians (2*pi rad = 360 deg)
    CRAD = PI_D / 180.0_DP
    DEGREE_RAD = DEGREE * CRAD
   ! Value of tand is equal to value of tan when the argument is given
   ! in radians
    TAND_D = Tan (DEGREE_RAD)
  End Function TAND_D
End Module trig_degrees
