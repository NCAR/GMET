subroutine exit_scrf (exit_status, outtxt)
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
  use nrtype ! data types
  implicit none
  character (len=*), intent (in) :: outtxt ! exit text
  integer (i4b), intent (in) :: exit_status ! exit status
! ----------------------------------------------------------------------------------------
! print out the exit text
  print *, outtxt
! also put the exit text in the log file
  write (99,*) outtxt
! exit with exit status
  call exit (exit_status)
! ----------------------------------------------------------------------------------------
end subroutine exit_scrf
