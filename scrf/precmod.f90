Module precision
 
! Real kinds
 
  Integer, Parameter :: kr4 = Selected_Real_Kind (6, 37)! single precision real
  Integer, Parameter :: kr8 = Selected_Real_Kind (15, 307)! double precision real
 
! Integer kinds
 
  Integer, Parameter :: ki4 = Selected_Int_Kind (9)! single precision integer
  Integer, Parameter :: ki8 = Selected_Int_Kind (18)! double precision integer
 
!Complex kinds
 
  Integer, Parameter :: kc4 = kr4 ! single precision complex
  Integer, Parameter :: kc8 = kr8 ! double precision complex
 
End Module precision
