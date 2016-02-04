FUNCTION erfinv (x)
  USE nrtype
  IMPLICIT NONE
 
  REAL (sp), INTENT (IN) :: x
  REAL (sp) :: erfinv
  REAL (sp) :: tmp
  INTEGER :: neg
 
  neg = 0
 
  IF (x .LT. 0.0) THEN
   neg = 1
   erfinv = - x
  ELSE
   erfinv = x
  END IF
 
  IF (erfinv .LE. 0.7) THEN
   tmp = erfinv * erfinv
   erfinv = erfinv * (((-0.140543331*tmp+0.914624893)*tmp-1.645349621)*tmp+0.886226899) / &
  & ((((0.012229801*tmp-0.329097515)*tmp+1.442710462)*tmp-2.118377725)*tmp+1.0)
  ELSE
   tmp = Sqrt (-Log(0.5*(1.0-erfinv)))
   erfinv = (((1.641345311*tmp+3.429567803)*tmp-1.624906493)*tmp-1.970840454) / ((1.637067800*tmp+3.543889200)*tmp+1.0)
  END IF
 
  IF (neg) THEN
   erfinv = - erfinv
  END IF
 
  RETURN
END FUNCTION erfinv
