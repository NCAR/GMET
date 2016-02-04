MODULE qpe_namelist
  USE nrtype
 
  IMPLICIT NONE
 
  INTEGER (I4B) :: nens !number of ensemble members to generate
  INTEGER (I4B) :: ntimes !number of times in timeseries to generate ensemble fields for
  INTEGER (I4B) :: start_time !time step to start at
 
  REAL (dp) :: clen !correlation length for scrf
 
  CHARACTER (LEN=1024) :: grid_name !name of ascii grid file used in jason's code
  CHARACTER (LEN=1024) :: out_name_base !base output name for netcdf files
  CHARACTER (LEN=1024) :: qpe_nc_name !name of netcdf output file from jason's code
 
 
  !define namelist required variables
  NAMELIST / PARAMS / nens, ntimes, grid_name, out_name_base, qpe_nc_name, clen, start_time
 
 
  SAVE
 
CONTAINS
 
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc                                                               ccc
!cc                      SUBROUTINES!!!!!!!!!                     ccc
!cc                                                               ccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!ccccccccccccccccccccccc
  SUBROUTINE read_namelist
 
   IMPLICIT NONE
 
!local variables
   INTEGER :: ierr
 
   OPEN (UNIT=30, FILE="namelist.params", FORM="FORMATTED")
 
   READ (UNIT=30, NML=PARAMS, IOSTAT=ierr)
   IF (ierr /= 0) THEN
    WRITE (*, '(/," ***** ERROR: Problem reading namelist PARAMS",/)')
    REWIND (UNIT=30)
    READ (UNIT=30, NML=PARAMS)
    STOP " ***** ERROR: Problem reading namelist PARAMS"
   END IF
 
   RETURN
  END SUBROUTINE
 
 
 
END MODULE qpe_namelist
