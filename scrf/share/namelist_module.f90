module namelist_module
  use nrtype
  implicit none
 
  ! integer (i4b) :: nens !number of ensemble members to generate
  integer (i4b) :: start_ens, stop_ens ! start and stop numbers of ens. members to make
  integer (i4b) :: ntimes !number of times in timeseries to generate ensemble fields for
  integer (i4b) :: start_time !time step to start at

  real (dp) :: clen !correlation length for scrf

  character(len=1024)   :: climo_path           !base path to climatological netcdf fields
  character (len=1024) :: grid_name !name of grid file
  character (len=1024) :: out_forc_name_base !base output name for netcdf forcing ens file
  character (len=1024) :: in_regr_name ! name of netcdf regression file -- input

  character (len=1024) :: time_mode     ! mode of ensemble generation (climo, daily climo, daily only)
 
  ! define namelist required variables
  ! namelist / params / nens, ntimes, grid_name, out_name_base, qpe_nc_name, clen, start_time
  namelist / params / start_ens, stop_ens, ntimes, grid_name, out_forc_name_base, in_regr_name, clen, &
                      start_time, climo_path, time_mode
 
  save
contains
 
    ! AWW-16 - updated to process namelist file given as argument
  subroutine read_namelist (namelist_filename)
    implicit none
 
      !input variables
    character (len=200), intent (in) :: namelist_filename !AWW-2016
 
      !local variables
    integer :: ierr
 
    open (unit=30, file=namelist_filename, form="FORMATTED")
 
    read (unit=30, nml=params, iostat=ierr)
    if (ierr /= 0) then
      write (*, '(/," ***** ERROR: Problem reading namelist PARAMS",/)')
      rewind (unit=30)
      read (unit=30, nml=params)
      stop " ***** ERROR: Problem reading namelist PARAMS"
    end if
 
    return
  end subroutine
 
end module namelist_module
