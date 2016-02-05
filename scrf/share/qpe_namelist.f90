Module qpe_namelist
  Use nrtype
 
  Implicit None
 
  Integer (I4B) :: nens !number of ensemble members to generate
  Integer (I4B) :: ntimes !number of times in timeseries to generate ensemble fields for
  Integer (I4B) :: start_time !time step to start at
 
  Real (dp) :: clen !correlation length for scrf
 
  Character (Len=1024) :: grid_name !name of ascii grid file used in jason's code
  Character (Len=1024) :: out_name_base !base output name for netcdf files
  Character (Len=1024) :: qpe_nc_name !name of netcdf output file from jason's code
 
 
  !define namelist required variables
  Namelist / PARAMS / nens, ntimes, grid_name, out_name_base, qpe_nc_name, clen, start_time
 
 
  Save
 
Contains
 
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc                                                               ccc
!cc                      SUBROUTINES!!!!!!!!!                     ccc
!cc                                                               ccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!ccccccccccccccccccccccc
  Subroutine read_namelist
 
    Implicit None
 
!local variables
    Integer :: ierr
 
    Open (Unit=30, File="namelist.params", Form="FORMATTED")
 
    Read (Unit=30, Nml=PARAMS, IoStat=ierr)
    If (ierr /= 0) Then
      Write (*, '(/," ***** ERROR: Problem reading namelist PARAMS",/)')
      Rewind (UNIT=30)
      Read (Unit=30, Nml=PARAMS)
      Stop " ***** ERROR: Problem reading namelist PARAMS"
    End If
 
    Return
  End Subroutine
 
 
 
End Module qpe_namelist
