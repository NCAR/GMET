! modified AWW Dec 2015
subroutine read_config (fname, n, names, values)
 
  use strings
  implicit none
 
  character (len=*) :: fname
  integer :: n
  character (len=500) :: names (n)
  character (len=500) :: values (n)
  character (len=500) :: settings (2)
  character (1000) line
  integer ipos, stat, nsettings, i
  logical fexist
 
  inquire (file=fname, exist=fexist)
  if ( .not. fexist) then
    print *, "Cannot find config file: ", trim (fname)
    return
  end if
 
  do i = 1, n, 1
    values (i) = ""
  end do
 
  open (11, file=fname, status='old', access='sequential')
 
  do
    read (11, "(A)", iostat=stat) line
    if (stat < 0) exit
    line = adjustl (line)
    if (line(1:1) .ne. "!" .and. len(trim(line)) /= 0) then
      ipos = index (line, '=')
      if (ipos > 0) then
        call parse (line, "=", settings, nsettings)
        do i = 1, n, 1
          if (index(settings(1), names(i)) == 1) values (i) = settings (2)
        end do
      end if
    end if
  end do
 
end subroutine read_config
