MODULE trig_degrees
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Einar Örn Hreinsson, Nov 2007 (einarhre@gmail.com)
  ! The family of functions sind, cosd and tand
  ! are not intrinsic to the gfortran compiler. These functions 
  ! are implemented here in terms of the instrinsic functions sin, 
  ! cos and tan.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE nrtype
  IMPLICIT NONE
  INTERFACE SIND
    MODULE PROCEDURE SIND_R, SIND_D
  END INTERFACE
  INTERFACE COSD
    MODULE PROCEDURE COSD_R, COSD_D
  END INTERFACE
  INTERFACE TAND
    MODULE PROCEDURE TAND_R, TAND_D
  END INTERFACE
CONTAINS 
  FUNCTION SIND_R(DEGREE)
  ! Input argument to sin in degrees
  REAL(SP), INTENT(IN)                 :: DEGREE     ! argument to sind in degrees
  ! Internal variable
  REAL(SP)                             :: CRAD       ! conversion from degrees to radians
  REAL(SP)                             :: DEGREE_RAD ! argument to sin in radians
  ! Output
  REAL(SP)                             :: SIND_R    ! the value of sin function
   ! Convert degrees to radians (2*pi rad = 360 deg)
   CRAD = PI/180.0_sp
   DEGREE_RAD = DEGREE*CRAD
   ! Value of SIND is equal to value of SIN when the argument is given
   ! in radians
   SIND_R = SIN(DEGREE_RAD) 
  END FUNCTION SIND_R
  FUNCTION SIND_D(DEGREE)
  ! Input argument to sin in degrees
  REAL(DP), INTENT(IN)                 :: DEGREE     ! argument to sind in degrees
  ! Internal variable
  REAL(DP)                             :: CRAD       ! conversion from degrees to radians
  REAL(DP)                             :: DEGREE_RAD ! argument to sin in radians
  ! Output
  REAL(DP)                             :: SIND_D    ! the value of sin function
   ! Convert degrees to radians (2*pi rad = 360 deg)
   CRAD = PI_D/180.0_dp
   DEGREE_RAD = DEGREE*CRAD
   ! Value of SIND is equal to value of SIN when the argument is given
   ! in radians
   SIND_D = SIN(DEGREE_RAD) 
  END FUNCTION SIND_D
  FUNCTION COSD_R(DEGREE)
  ! Input argument to cosd in degrees
  REAL(SP), INTENT(IN)                 :: DEGREE     ! argument to cosd in degrees
  ! Internal variable
  REAL(SP)                             :: CRAD       ! conversion from degrees to radians
  REAL(SP)                             :: DEGREE_RAD ! argument to cos in radians
  ! Output
  REAL(SP)                             :: COSD_R       ! the value of cos function
   ! Convert degrees to radians (2*pi rad = 360 deg)
   CRAD = PI/180.0_sp
   DEGREE_RAD = DEGREE*CRAD
   ! Value of COSD is equal to value of COS when the argument is given
   ! in radians
   COSD_R = COS(DEGREE_RAD) 
  END FUNCTION COSD_R
  FUNCTION COSD_D(DEGREE)
  ! Input argument to cosd in degrees
  REAL(DP), INTENT(IN)                 :: DEGREE     ! argument to cosd in degrees
  ! Internal variable
  REAL(DP)                             :: CRAD       ! conversion from degrees to radians
  REAL(DP)                             :: DEGREE_RAD ! argument to cos in radians
  ! Output
  REAL(DP)                             :: COSD_D       ! the value of cos function
   ! Convert degrees to radians (2*pi rad = 360 deg)
   CRAD = PI_D/180.0_dp
   DEGREE_RAD = DEGREE*CRAD
   ! Value of COSD is equal to value of COS when the argument is given
   ! in radians
   COSD_D = COS(DEGREE_RAD) 
  END FUNCTION COSD_D
  FUNCTION TAND_R(DEGREE)
  ! Input argument to tand in degrees
  REAL(SP), INTENT(IN)                 :: DEGREE     ! argument to tand in degrees
  ! Internal variable
  REAL(SP)                             :: CRAD       ! conversion from degrees to radians
  REAL(SP)                             :: DEGREE_RAD ! argument to tan in radians
  ! Output
  REAL(SP)                             :: TAND_R       ! the value of tan function
   ! Convert degrees to radians (2*pi rad = 360 deg)
   CRAD = PI/180.0_sp
   DEGREE_RAD = DEGREE*CRAD
   ! Value of tand is equal to value of tan when the argument is given
   ! in radians
   TAND_R = TAN(DEGREE_RAD) 
  END FUNCTION TAND_R
  FUNCTION TAND_D(DEGREE)
  ! Input argument to tand in degrees
  REAL(DP), INTENT(IN)                 :: DEGREE     ! argument to tand in degrees
  ! Internal variable
  REAL(DP)                             :: CRAD       ! conversion from degrees to radians
  REAL(DP)                             :: DEGREE_RAD ! argument to tan in radians
  ! Output
  REAL(DP)                             :: TAND_D       ! the value of tan function
   ! Convert degrees to radians (2*pi rad = 360 deg)
   CRAD = PI_D/180.0_dp
   DEGREE_RAD = DEGREE*CRAD
   ! Value of tand is equal to value of tan when the argument is given
   ! in radians
   TAND_D = TAN(DEGREE_RAD) 
  END FUNCTION TAND_D
END MODULE trig_degrees
