subroutine read_grid_list(file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
  use strings
  use nrtype
  implicit none

  character(len=500), intent(in) :: file_name
  real(DP), allocatable, intent(out) :: lats(:), lons(:), alts(:), slp_n(:), slp_e(:)
  integer(I4B), intent(out) :: nx, ny
  integer, intent(out) :: error

  character(len=100) :: peices(5)
  integer :: ngrid, i, npeices, stat, err, ipos
  real(DP) :: dx, dy, startx, starty
  character (300) :: line
  logical fexist

  ngrid = 0
  nx = 0
  ny = 0
  dx = 0.0
  dy = 0.0
  startx = 0.0
  starty = 0.0

  error = 0
  print *, "Reading grid file: ", trim(file_name)

  inquire (file=file_name,exist=fexist)
  if(.not.fexist) then
     print *, "Cannot find file: ", trim(file_name)
     error = 1
     return
  endif

  open (13,file=file_name,status='old')

  i = 1
  do
     read(13, "(A)", iostat=stat) line
     if(stat < 0) exit
     line = adjustl(line)
     if(line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) then
        
        ipos=index(line,"NX")
        if(ipos > 0) then
           ipos = ipos + 2
           call value(line(ipos:), nx, err)
           if(nx > 0 .AND. ny > 0) then
              ngrid = nx*ny
              allocate(lats(ngrid))
              allocate(lons(ngrid))
              allocate(alts(ngrid))
              allocate(slp_n(ngrid))
	      allocate(slp_e(ngrid))
           endif
        endif
        ipos=index(line,"NY")
        if(ipos > 0) then
           ipos = ipos + 2
           call value(line(ipos:), ny, err)
           if(nx > 0 .AND. ny > 0) then
              ngrid = nx*ny
              allocate(lats(ngrid))
              allocate(lons(ngrid))
              allocate(alts(ngrid))
	      allocate(slp_n(ngrid))
	      allocate(slp_e(ngrid))
           endif
        endif
        ipos=index(line,"DX")
        if(ipos > 0) then
           ipos = ipos + 2
           call value(line(ipos:), dx, err)
        endif
        ipos=index(line,"DY")
        if(ipos > 0) then
           ipos = ipos + 2
           call value(line(ipos:), dy, err)
        endif
        ipos=index(line,"STARTX")
        if(ipos > 0) then
           ipos = ipos + 6
           call value(line(ipos:), startx, err)
        endif
        ipos=index(line,"STARTY")
        if(ipos > 0) then
           ipos = ipos + 6
           call value(line(ipos:), starty, err)
        endif
     endif

     ipos=index(line,',')
     if(ipos > 0 .AND. ngrid > 0) then
        call parse(line, ",", peices, npeices)
        if(npeices == 5) then
           call value(peices(1), lats(i), err)
           if(err /= 0) lats(i) = -999.99
           call value(peices(2), lons(i), err)
           if(err /= 0) lons(i) = -999.99
           call value(peices(3), alts(i), err)
           if(err /= 0) alts(i) = -999.99
	   call value(peices(4), slp_n(i), err)
           if(err /= 0) slp_n(i) = -999.99
	   call value(peices(5), slp_e(i), err)
           if(err /= 0) slp_e(i) = -999.99
           i = i + 1
        endif
     endif

  end do
  if(ngrid == 0) then
     print *, "Failed to find NX, NY in grid list: ", trim(file_name)
     error = 1
  else
     if(i /= ngrid+1) then
        print *, "Found only ", i, " out of ", ngrid, " points from: ", trim(file_name)
        error = 1
     endif
  endif

end subroutine read_grid_list
