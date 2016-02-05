Subroutine EXIT_SCRF (EXIT_STATUS, OUTTXT)
! ----------------------------------------------------------------------------------------
! Creator:
!   Einar Örn Hreinsson, 2009
!
! Modified:
!   Andy Newman, Aug 2012
!
! ----------------------------------------------------------------------------------------
! Purpose:
!   To exit from spatially correlated random fields with an exit status, which can then be examined by
!   looking at the $? variable in the shell that the code was called from
!
! ----------------------------------------------------------------------------------------
  Use nrtype ! data types
  Implicit None
  Character (Len=*), Intent (In) :: OUTTXT ! exit text
  Integer (I4B), Intent (In) :: EXIT_STATUS ! exit status
! ----------------------------------------------------------------------------------------
! print out the exit text
  Print *, OUTTXT
! also put the exit text in the log file
  Write (99,*) OUTTXT
! exit with exit status
  Call EXIT (EXIT_STATUS)
! ----------------------------------------------------------------------------------------
End Subroutine EXIT_SCRF
