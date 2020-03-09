subroutine read_climo_temp(climo_path,current_month,climo_field,climo_out,error)
  use netcdf
  use nrtype
  implicit none

  character (len = *), intent(in)     :: climo_path      !path to climo netcdf files
  integer, intent(in)                 :: current_month   !integer month
  character (len = *), intent(in)     :: climo_field     !name of climo variable to read
  real(SP), intent(out)               :: climo_out(:,:)    !climo variable grid
  integer, intent(out)                :: error           !error code


  !local variables
  character (len = 5000)                     :: climo_file     ! name of climo netcdf file
  character (len = 2)                        :: mnth_str       ! string of current month
  integer(I4B)                               :: ncid           ! id of netcdf file
  integer(I4B)                               :: varid          ! variable id from netcdf file
  integer(I4B)                               :: dimid          ! variable id from netcdf file
  integer(I4B)                               :: ndims          ! number of dimensions
  integer(I4B)                               :: nx             ! x-dimension
  integer(I4B)                               :: ny             ! y-dimension
  
  !code starts below

  !create string of month
  write(mnth_str,"(I0.2)") current_month
  
  !generate file name
  climo_file=trim(climo_path)//'/ENS_CLIMO_'//mnth_str//'.nc'

  
  !open file  
  call check(nf90_open(trim(climo_file),nf90_nowrite,ncid),"Climo file open error: "//trim(climo_file), error)
  if(error /= 0) return

  !inquire about variables
  call check(nf90_inq_varid(ncid,climo_field,varid),"Climo variable name error",error)
  if(error /= 0) return  

  !get dimensions
  call check(nf90_inquire_variable(ncid,varid,ndims = ndims),"Grid dimension inq error",error)
  if(error /= 0) return
  if(ndims/=2) then; print *,'Wrong number of dimensions',ndims,'should be 2'; return; endif

  !get x,y dimensions
  call check(nf90_inq_dimid(ncid,'y',dimid),"y dim inquiry error",error)
  call check(nf90_inquire_dimension(ncid,dimid,len = nx),"y dim error",error)

  if(error .ne. 0) then !check for lat
    error = 0
    call check(nf90_inq_dimid(ncid,'lat',dimid),"lat dim inquiry error",error)
    call check(nf90_inquire_dimension(ncid,dimid,len = ny),"lat dim error",error)  
  endif
  
  error = 0
  call check(nf90_inq_dimid(ncid,'x',dimid),"x dim inquiry error",error)
  call check(nf90_inquire_dimension(ncid,dimid,len = ny),"x dim error",error)
  if(error .ne. 0) then !check for lon
    error = 0
    call check(nf90_inq_dimid(ncid,'lon',dimid),"lon dim inquiry error",error)
    call check(nf90_inquire_dimension(ncid,dimid,len = nx),"lon dim error",error)
  endif
  
  !get it
  call check(nf90_get_var(ncid,varid,climo_out),"Climo variable read error",error)
  if(error /= 0) return
  
  !Close the netcdf file
  error = nf90_close(ncid)

 
contains
  subroutine check(status, info, error)
    integer, intent (in) :: status
    character (len = *), intent(in) :: info
    integer, intent(out) :: error

    if(status /= nf90_noerr) then
       print *, trim(info)//": "// trim(nf90_strerror(status))
       error = 1
    end if
  end subroutine check  

end subroutine read_climo_temp
