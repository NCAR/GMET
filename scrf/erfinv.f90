Function erfinv (x)
  Use nrtype
  Implicit None
 
  Real (sp), Intent (In) :: x
  Real (sp) :: erfinv
  Real (sp) :: tmp
  Integer :: neg
 
  neg = 0
 
  If (x .Lt. 0.0) Then
    neg = 1
    erfinv = - x
  Else
    erfinv = x
  End If
 
  If (erfinv .Le. 0.7) Then
    tmp = erfinv * erfinv
    erfinv = erfinv * (((-0.140543331*tmp+0.914624893)*tmp-1.645349621)*tmp+0.886226899) / &
   & ((((0.012229801*tmp-0.329097515)*tmp+1.442710462)*tmp-2.118377725)*tmp+1.0)
  Else
    tmp = Sqrt (-Log(0.5*(1.0-erfinv)))
    erfinv = (((1.641345311*tmp+3.429567803)*tmp-1.624906493)*tmp-1.970840454) / &
   & ((1.637067800*tmp+3.543889200)*tmp+1.0)
  End If
 
  If (neg) Then
    erfinv = - erfinv
  End If
 
  Return
End Function erfinv
