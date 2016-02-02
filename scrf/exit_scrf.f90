SUBROUTINE EXIT_SCRF(EXIT_STATUS,OUTTXT)
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
USE nrtype                                                  ! data types
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN)                  :: OUTTXT      ! exit text
INTEGER(I4B),INTENT(IN)                      :: EXIT_STATUS ! exit status
! ----------------------------------------------------------------------------------------
! print out the exit text
print *,OUTTXT
! also put the exit text in the log file
write (99,*) OUTTXT
! exit with exit status
CALL EXIT(EXIT_STATUS)
! ----------------------------------------------------------------------------------------
END SUBROUTINE EXIT_SCRF
