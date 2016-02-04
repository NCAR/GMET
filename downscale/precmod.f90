MODULE precision
 
! Real kinds
 
  INTEGER, PARAMETER :: kr4 = Selected_Real_Kind (6, 37)! single precision real
  INTEGER, PARAMETER :: kr8 = Selected_Real_Kind (15, 307)! double precision real
 
! Integer kinds
 
  INTEGER, PARAMETER :: ki4 = Selected_Int_Kind (9)! single precision integer
  INTEGER, PARAMETER :: ki8 = Selected_Int_Kind (18)! double precision integer
 
!Complex kinds
 
  INTEGER, PARAMETER :: kc4 = kr4 ! single precision complex
  INTEGER, PARAMETER :: kc8 = kr8 ! double precision complex
 
END MODULE precision
