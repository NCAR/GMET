
subroutine read_config(fname, n, names, values)

  use strings
  implicit none

  character(len=*) :: fname
  integer :: n
  character(len=500) :: names(n)
  character(len=500) :: values(n)


  character(len=500) :: peices(2)
  character (1000) line
  integer ipos, stat, npeices, i
  logical fexist

  inquire (file=fname,exist=fexist)
  if(.not.fexist) then
     print *, "Cannot find config file: ", trim(fname)
     return
  endif

  do i = 1, n, 1
     values(i) = ""
  enddo

  open (11,file=fname,status='old')

  do
     read(11, "(A)", iostat=stat) line
     if(stat < 0) exit
     line = adjustl(line)
     if(line(1:1) .NE. "!" .AND. len(trim(line)) /= 0) then
        ipos=index(line,'=')
        if(ipos > 0) then
           call parse(line, "=", peices, npeices)
           do i = 1, n, 1
              if(index(peices(1), names(i)) == 1) values(i) = peices(2)
           enddo
        endif

     endif

  end do

end subroutine read_config
