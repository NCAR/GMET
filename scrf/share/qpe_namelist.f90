module qpe_namelist
  use nrtype

  implicit none

  integer(I4B)		:: nens			!number of ensemble members to generate
  integer(I4B)		:: ntimes		!number of times in timeseries to generate ensemble fields for
  integer(I4B)		:: start_time		!time step to start at

  real(dp)		:: clen			!correlation length for scrf

  character(len=1024)	:: grid_name		!name of ascii grid file used in jason's code
  character(len=1024)	:: out_name_base	!base output name for netcdf files
  character(len=1024)	:: qpe_nc_name		!name of netcdf output file from jason's code


  !define namelist required variables
  namelist / PARAMS / nens,ntimes,grid_name,out_name_base,qpe_nc_name,clen,start_time


  save

  contains


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc                                                               ccc
!cc                      SUBROUTINES!!!!!!!!!                     ccc
!cc 								  ccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccc
subroutine read_namelist
  
  implicit none

!local variables
  integer :: ierr

  open(UNIT=30, file="namelist.params",form="FORMATTED")

  read(UNIT=30, NML=PARAMS, iostat=ierr)
  if (ierr /= 0) then
    write(*,'(/," ***** ERROR: Problem reading namelist PARAMS",/)')
    rewind(UNIT=30)
    read(UNIT=30, NML=PARAMS)
    stop " ***** ERROR: Problem reading namelist PARAMS"
  endif

  return
end subroutine



end module qpe_namelist