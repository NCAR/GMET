
subroutine read_refcst(startdate, enddate, file_var, perturbation, var_name, forecast, V, X, Y, T, error)
  use strings
  use utim
  use type
  implicit none

  interface
    subroutine read_nc_file(file_name, valid_time, var_name, var, lats, lons, error)
      use type
      character (len = *), intent(in) :: file_name
      character (len = *), intent(in) :: var_name
      real(DP), intent(in) :: valid_time
      real(DP), allocatable, intent(out) :: var(:, :, :)
      real(DP), allocatable, intent(out) :: lats(:), lons(:)
      integer, intent(out) :: error
    end subroutine read_nc_file
  end interface

  character (len = 100), intent(in) :: startdate, enddate, file_var, perturbation
  character (len = *), intent(in) :: var_name
  integer(I4B), intent(in) :: forecast
  real(DP), allocatable, intent(out) :: V(:,:), X(:), Y(:)
  real(DP), allocatable, intent(out) :: T(:)
  integer, intent(out) :: error

  real(DP), allocatable :: lats(:), lons(:), var(:, :, :)
  real(DP) :: valid_time, utime  
  character (len = 500) :: file
  character (len = 100) :: date
  integer(I4B) :: vshape(3)
  integer(I4B) :: sday, eday, ntimes, ngrids
  integer(I4B) :: sec, min, hour, day, month, year
  integer(I4B) :: i, j, k

  error = 0

  call parse_date(startdate, year, month, day, hour, min, sec, error)
  sday = julian_date(day, month, year)
  call parse_date(enddate, year, month, day, hour, min, sec, error)
  eday = julian_date(day, month, year)
  ntimes = eday - sday +1

  allocate(T(ntimes))

  i = 1
  utime= date_to_unix(startdate)
  do 
     if(utime > date_to_unix(enddate)) exit
     call unix_to_date(utime, year, month, day, hour, min, sec)

     write(date,"(I4.4I2.2I2.2)"), year, month, day
     if(forecast >= 190) then
        file = trim(date)//"/"//trim(file_var)//"_"//trim(date)//"00_"//trim(perturbation)//"_t190.nc" 
     else
        file = trim(date)//"/"//trim(file_var)//"_"//trim(date)//"00_"//trim(perturbation)//".nc" 
     endif


     valid_time = utime + (forecast*3600.0)

     call read_nc_file(file, valid_time, var_name, var, lats, lons, error)

     if(trim(var_name) == 'APCP_ens_mean_surface') then
	var = var**(1.0/4.0)
	print *,'normalizing GEFS precip'
     endif

     if(error == 1) then 
        print *, "Failed reading file ", trim(file)
	print *,'var ',trim(var_name)
        exit
     else
        print ("(AAAF11.0)"), "Success reading file ", trim(file), " Valid at ", valid_time
     endif

     vshape = shape(var)
     if(vshape(1) /= size(lons) .OR. vshape(2) /= size(lats) .OR. vshape(3) /= 1) then
        print *, "Dimensions from file do not match"
        error = 1
        exit
     endif

     if(i == 1) then
        ngrids = vshape(1)*vshape(2)
        allocate(V(ngrids, ntimes))
        allocate(X(vshape(2)))
        allocate(Y(vshape(1)))
        X(:) = lats(:)
        do j = 1, vshape(1), 1
          Y(j) = mod(lons(j),360.0)-360.0
       enddo

     else
        if(vshape(2) /= size(X) .OR. vshape(1) /= size(Y)) then
           print *, "Dimensions from file do not match previous files"
           error = 1
           exit
        endif
     endif

     T(i) = valid_time
     V(1:ngrids,i) = reshape(var, (/ngrids/))

     i = i + 1
     utime = utime+86400

  enddo

end subroutine read_refcst


subroutine read_station_list(file_name, id, name, lat, lon, alt, sslp_n, sslp_e, n_stations, error)
  use strings
  use type
  implicit none

  character(len=500), intent(in) :: file_name
  character(len=100), allocatable, intent(out) :: id(:), name(:)
  real(DP), allocatable, intent(out) :: lat(:), lon(:), alt(:), sslp_n(:), sslp_e(:)
  integer(I4B), intent(out) :: n_stations
  integer, intent(out) :: error

  character(len=100) :: peices(7)
  integer(I4B) i, npeices, stat, err, ipos
  character (500) line
  logical fexist

  error = 0
  print *, "Reading stations list: ", trim(file_name)

  inquire (file=file_name,exist=fexist)
  if(.not.fexist) then
     print *, "Cannot find file: ", trim(file_name)
     error = 1
     return
  endif

  open (11,file=file_name,status='old')

  n_stations = 0
  i = 1
  do
     read(11, "(A)", iostat=stat) line
     if(stat < 0) exit
     line = adjustl(line)

     if(line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) then
        ipos=index(line,"NSITES")
        if(ipos > 0) then
           ipos = ipos + 6
           call value(line(ipos:), n_stations, err)
           if(err /= 0) then
              n_stations = 0
           else
              print *, "Stations: ", n_stations
              allocate(id(n_stations))
              allocate(name(n_stations))
              allocate(lat(n_stations))
              allocate(lon(n_stations))
              allocate(alt(n_stations))
              allocate(sslp_n(n_stations))
              allocate(sslp_e(n_stations))
           endif
        endif
        ipos=index(line,',')
        if(ipos > 0 .AND. n_stations > 0) then
           call parse(line, ",", peices, npeices)
           if(npeices == 7) then
              id(i) = peices(1)
              name(i) = peices(7)
              call delall(name(i), '"')
              call value(peices(2), lat(i), err)
              if(err /= 0) lat(i) = -999.99
              call value(peices(3), lon(i), err)
              if(err /= 0) lon(i) = -999.99
              call value(peices(4), alt(i), err)
              if(err /= 0) alt(i) = -999.99
              call value(peices(5), sslp_n(i), err)
              if(err /= 0) sslp_n(i) = -999.99
              call value(peices(6), sslp_e(i), err)
              if(err /= 0) sslp_e(i) = -999.99

              i = i + 1
           endif
        endif

     endif

  end do
  if(n_stations == 0) then
     print *, "Failed to find NSITES in station list: ", trim(file_name)
     error = 1
  else
     if(i /= n_stations+1) then
        print *, "Found only ", i, " out of ", n_stations, " stations from: ", trim(file_name)
        error = 1
     endif
  endif

end subroutine read_station_list


subroutine read_station(stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss, vals_miss_t, error)
  use strings
  use utim
  use type
  implicit none

  character(len=100), intent(in) :: stnvar
  character(len=100), intent(in) :: stnid
  character(len=100), intent(in) :: site_var, site_var_t
  character(len=500), intent(in) :: site_list
  real(DP), intent(in) :: Times(:)
  real(DP), allocatable, intent(out) :: vals(:), tair_vals(:,:)
  logical, allocatable, intent(out) :: vals_miss(:),vals_miss_t(:)
  integer, intent(out) :: error

  character (len = 500) :: file_name
  character (len = 500) :: directory
  character (len = 500) :: line
  character (len = 100) :: peices(10)
  character (len = 100) :: stnvar_t
  real(DP) :: utime  
  integer(I4B) :: i, npeices, ntimes, ipos, stat, var_index, val_index
  integer(I4B) :: sec, min, hour, day, month, year
  logical :: fexist


!!!! first read in station precipitation data
  error = 0

  call parse(site_list, "/", peices, npeices)

  if(npeices == 0) then 
     directory = "."
  else
     directory = peices(1)
     if(npeices > 2) then
        do i = 2, npeices-1, 1
           directory = trim(directory)//"/"//peices(i)
        enddo
     endif
  endif
  file_name = trim(directory)//"/"//trim(stnid)//"_"//trim(site_var)//".txt"
  print *, "Reading station file: ", trim(file_name)

  ntimes = size(Times)
  allocate(vals(ntimes))
  allocate(vals_miss(ntimes))
  vals_miss = .FALSE.
  vals = -999.0

  inquire (file=file_name,exist=fexist)

  if(.not.fexist) then
     print *, "Cannot find file: ", trim(file_name)
     print *, "Setting station precipitation timeseries to missing"

  else
    open (12,file=file_name,status='old')

    var_index = -1
    i = 1
    do
      read(12, "(A)", iostat=stat) line
      if(stat < 0) exit
      line = adjustl(line)

      if(line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) then

	  if(var_index == -1) then
	    ipos=index(line,"DATE")
	    if(ipos > 0) then
		call parse(line, " ", peices, npeices)
		if(npeices < 3 .OR. trim(peices(1)) .NE. "DATE" .OR. trim(peices(2)) .NE. "HHMMSS") then
		  print *, "Failed to read header from file ", trim(file_name)
		  error = 1
		  return
		endif
		do i = 3, npeices, 1
		  if(trim(peices(i)) .EQ. trim(stnvar)) then
		      var_index = i
		  endif
		enddo
	    endif

	  else

	    call parse(line, " ", peices, npeices)

	    if(npeices < 3) then
		print *, "Confusing line from file ", trim(file_name), " line: ", trim(line)
		error = 1
		return
	    endif

	    utime = date_to_unix(trim(peices(1))//trim(peices(2)))
	    if(utime > Times(1) - 43200.0 .AND. utime <= Times(ntimes) + 43200.0) then
		val_index = FLOOR((utime - Times(1)) / 86400.0) + 1
		vals_miss(val_index) = .TRUE.
		call value(peices(var_index), vals(val_index), error)
		if(error /= 0) print *, "Problem converting string to float: ", peices(var_index)
		error = 0
		if(vals(val_index) == -999.0) then
		  vals_miss(val_index) = .FALSE.
		endif
	    endif
	  endif

      endif
    end do

    if(var_index == -1) then
      print *, "Failed to find header from file ", trim(file_name)
      error = 1
    endif
  endif  !end file exist if statement

!!! read in station temperature data now
  error = 0
  
  stnvar_t = 'Tmin'
  ntimes = size(Times)
  allocate(tair_vals(2,ntimes))
  allocate(vals_miss_t(ntimes))
  vals_miss_t = .FALSE.
  tair_vals = -999.0

  call parse(site_list, "/", peices, npeices)

  if(npeices == 0) then 
     directory = "."
  else
     directory = peices(1)
     if(npeices > 2) then
        do i = 2, npeices-1, 1
           directory = trim(directory)//"/"//peices(i)
        enddo
     endif
  endif

  file_name = trim(directory)//"/"//trim(stnid)//"_"//trim(site_var_t)//".txt"

  print *, "Reading station file: ", trim(file_name)

  inquire (file=file_name,exist=fexist)
  if(.not.fexist) then
     print *, "Cannot find file: ", trim(file_name)
     print *, "Setting station temperature timeseries to missing"

  else
    open (12,file=file_name,status='old')

    var_index = -1
    i = 1
    do
      read(12, "(A)", iostat=stat) line
      if(stat < 0) exit
      line = adjustl(line)

      if(line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) then

	  if(var_index == -1) then
	    ipos=index(line,"DATE")
	    if(ipos > 0) then
		call parse(line, " ", peices, npeices)
		if(npeices < 3 .OR. trim(peices(1)) .NE. "DATE" .OR. trim(peices(2)) .NE. "HHMMSS") then
		  print *, "Failed to read header from file ", trim(file_name)
		  error = 1
		  return
		endif
		do i = 3, npeices, 1
		  if(trim(peices(i)) .EQ. trim(stnvar_t)) then
		      var_index = i
		  endif
		enddo
	    endif

	  else

	    call parse(line, " ", peices, npeices)

	    if(npeices < 3) then
		print *, "Confusing line from file ", trim(file_name), " line: ", trim(line)
		error = 1
		return
	    endif

	    utime = date_to_unix(trim(peices(1))//trim(peices(2)))
	    if(utime > Times(1) - 43200.0 .AND. utime <= Times(ntimes) + 43200.0) then
!		val_index = NINT((utime - Times(1)) / 86400.0) + 1
		val_index = FLOOR((utime - Times(1)) / 86400.0) + 1
		vals_miss_t(val_index) = .TRUE.
		call value(peices(var_index), tair_vals(1,val_index), error)
		call value(peices(var_index+1),tair_vals(2,val_index),error)

		if(error /= 0) print *, "Problem converting string to float: ", peices(var_index)
		error = 0

		if(tair_vals(1,val_index) == -999.0 .or. tair_vals(2,val_index) == -999.0) then
		  vals_miss_t(val_index) = .FALSE.
		endif
	    endif
	  endif

      endif
    end do

    if(var_index == -1) then
      print *, "Failed to find header from file ", trim(file_name)
      error = 1
    endif
  endif  !end file exist if statement

end subroutine read_station



subroutine read_grid_list(file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
  use strings
  use type
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
