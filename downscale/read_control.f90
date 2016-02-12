 
Subroutine read_refcst (startdate, enddate, file_var, perturbation, var_name, forecast, V, X, Y, T, error)
  Use strings
  Use utim
  Use type
  Implicit None
 
  Interface
    Subroutine read_nc_file (file_name, valid_time, var_name, var, lats, lons, error)
      Use type
      Character (Len=*), Intent (In) :: file_name
      Character (Len=*), Intent (In) :: var_name
      Real (DP), Intent (In) :: valid_time
      Real (DP), Allocatable, Intent (Out) :: var (:, :, :)
      Real (DP), Allocatable, Intent (Out) :: lats (:), lons (:)
      Integer, Intent (Out) :: error
    End Subroutine read_nc_file
  End Interface
 
  Character (Len=100), Intent (In) :: startdate, enddate, file_var, perturbation
  Character (Len=*), Intent (In) :: var_name
  Integer (I4B), Intent (In) :: forecast
  Real (DP), Allocatable, Intent (Out) :: V (:, :), X (:), Y (:)
  Real (DP), Allocatable, Intent (Out) :: T (:)
  Integer, Intent (Out) :: error
 
  Real (DP), Allocatable :: lats (:), lons (:), var (:, :, :)
  Real (DP) :: valid_time, utime
  Character (Len=500) :: file
  Character (Len=100) :: date
  Integer (I4B) :: vshape (3)
  Integer (I4B) :: sday, eday, ntimes, ngrids
  Integer (I4B) :: sec, Min, hour, day, month, year
  Integer (I4B) :: i, j, k
 
  error = 0
 
  Call parse_date (startdate, year, month, day, hour, Min, sec, error)
  sday = julian_date (day, month, year)
  Call parse_date (enddate, year, month, day, hour, Min, sec, error)
  eday = julian_date (day, month, year)
  ntimes = eday - sday + 1
 
  Allocate (T(ntimes))
 
  i = 1
  utime = date_to_unix (startdate)
  Do
    If (utime > date_to_unix(enddate)) Exit
    Call unix_to_date (utime, year, month, day, hour, Min, sec)
 
    Write (date, "(I4.4I2.2I2.2)"), year, month, day
    If (forecast >= 190) Then
      file = trim (date) // "/" // trim (file_var) // "_" // trim (date) // "00_" // trim (perturbation) // "_t190.nc"
    Else
      file = trim (date) // "/" // trim (file_var) // "_" // trim (date) // "00_" // trim (perturbation) // ".nc"
    End If
 
 
    valid_time = utime + (forecast*3600.0)
 
    Call read_nc_file (file, valid_time, var_name, var, lats, lons, error)
 
    If (trim(var_name) == 'APCP_ens_mean_surface') Then
      var = var ** (1.0/4.0)
      Print *, 'normalizing GEFS precip'
    End If
 
    If (error == 1) Then
      Print *, "Failed reading file ", trim (file)
      Print *, 'var ', trim (var_name)
      Exit
    Else
      Print ("(AAAF11.0)"), "Success reading file ", trim (file), " Valid at ", valid_time
    End If
 
    vshape = shape (var)
    If (vshape(1) /= size(lons) .Or. vshape(2) /= size(lats) .Or. vshape(3) /= 1) Then
      Print *, "Dimensions from file do not match"
      error = 1
      Exit
    End If
 
    If (i == 1) Then
      ngrids = vshape (1) * vshape (2)
      Allocate (V(ngrids, ntimes))
      Allocate (X(vshape(2)))
      Allocate (Y(vshape(1)))
      X (:) = lats (:)
      Do j = 1, vshape (1), 1
        Y (j) = Mod (lons(j), 360.0) - 360.0
      End Do
 
    Else
      If (vshape(2) /= size(X) .Or. vshape(1) /= size(Y)) Then
        Print *, "Dimensions from file do not match previous files"
        error = 1
        Exit
      End If
    End If
 
    T (i) = valid_time
    V (1:ngrids, i) = reshape (var, (/ ngrids /))
 
    i = i + 1
    utime = utime + 86400
 
  End Do
 
End Subroutine read_refcst
 
 
Subroutine read_station_list (file_name, id, name, lat, lon, alt, sslp_n, sslp_e, n_stations, error)
  Use strings
  Use type
  Implicit None
 
  Character (Len=500), Intent (In) :: file_name
  Character (Len=100), Allocatable, Intent (Out) :: id (:), name (:)
  Real (DP), Allocatable, Intent (Out) :: lat (:), lon (:), alt (:), sslp_n (:), sslp_e (:)
  Integer (I4B), Intent (Out) :: n_stations
  Integer, Intent (Out) :: error
 
  Character (Len=100) :: peices (7)
  Integer (I4B) i, npeices, stat, err, ipos
  Character (500) line
  Logical fexist
 
  error = 0
  Print *, "Reading stations list: ", trim (file_name)
 
  Inquire (File=file_name, Exist=fexist)
  If ( .Not. fexist) Then
    Print *, "Cannot find file: ", trim (file_name)
    error = 1
    Return
  End If
 
  Open (11, File=file_name, Status='old')
 
  n_stations = 0
  i = 1
  Do
    Read (11, "(A)", IoStat=stat) line
    If (stat < 0) Exit
    line = adjustl (line)
 
    If (line(1:1) .Ne. "#" .And. len(trim(line)) /= 0) Then
      ipos = index (line, "NSITES")
      If (ipos > 0) Then
        ipos = ipos + 6
        Call value (line(ipos:), n_stations, err)
        If (err /= 0) Then
          n_stations = 0
        Else
          Print *, "Stations: ", n_stations
          Allocate (id(n_stations))
          Allocate (name(n_stations))
          Allocate (lat(n_stations))
          Allocate (lon(n_stations))
          Allocate (alt(n_stations))
          Allocate (sslp_n(n_stations))
          Allocate (sslp_e(n_stations))
        End If
      End If
      ipos = index (line, ',')
      If (ipos > 0 .And. n_stations > 0) Then
        Call parse (line, ",", peices, npeices)
        If (npeices == 7) Then
          id (i) = peices (1)
          name (i) = peices (7)
          Call delall (name(i), '"')
          Call value (peices(2), lat(i), err)
          If (err /= 0) lat (i) = - 999.99
          Call value (peices(3), lon(i), err)
          If (err /= 0) lon (i) = - 999.99
          Call value (peices(4), alt(i), err)
          If (err /= 0) alt (i) = - 999.99
          Call value (peices(5), sslp_n(i), err)
          If (err /= 0) sslp_n (i) = - 999.99
          Call value (peices(6), sslp_e(i), err)
          If (err /= 0) sslp_e (i) = - 999.99
 
          i = i + 1
        End If
      End If
 
    End If
 
  End Do
  If (n_stations == 0) Then
    Print *, "Failed to find NSITES in station list: ", trim (file_name)
    error = 1
  Else
    If (i /= n_stations+1) Then
      Print *, "Found only ", i, " out of ", n_stations, " stations from: ", trim (file_name)
      error = 1
    End If
  End If
 
End Subroutine read_station_list
 
 
Subroutine read_station (stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss, &
& vals_miss_t, error)
  Use strings
  Use utim
  Use type
  Implicit None
 
  Character (Len=100), Intent (In) :: stnvar
  Character (Len=100), Intent (In) :: stnid
  Character (Len=100), Intent (In) :: site_var, site_var_t
  Character (Len=500), Intent (In) :: site_list
  Real (DP), Intent (In) :: Times (:)
  Real (DP), Allocatable, Intent (Out) :: vals (:), tair_vals (:, :)
  Logical, Allocatable, Intent (Out) :: vals_miss (:), vals_miss_t (:)
  Integer, Intent (Out) :: error
 
  Character (Len=500) :: file_name
  Character (Len=500) :: directory
  Character (Len=500) :: line
  Character (Len=100) :: peices (10)
  Character (Len=100) :: stnvar_t
  Real (DP) :: utime
  Integer (I4B) :: i, npeices, ntimes, ipos, stat, var_index, val_index
  Integer (I4B) :: sec, Min, hour, day, month, year
  Logical :: fexist
 
 
!!!! first read in station precipitation data
  error = 0
 
  Call parse (site_list, "/", peices, npeices)
 
  If (npeices == 0) Then
    directory = "."
  Else
    directory = peices (1)
    If (npeices > 2) Then
      Do i = 2, npeices - 1, 1
        directory = trim (directory) // "/" // peices (i)
      End Do
    End If
  End If
  file_name = trim (directory) // "/" // trim (stnid) // "_" // trim (site_var) // ".txt"
  Print *, "Reading station file: ", trim (file_name)
 
  ntimes = size (Times)
  Allocate (vals(ntimes))
  Allocate (vals_miss(ntimes))
  vals_miss = .False.
  vals = - 999.0
 
  Inquire (File=file_name, Exist=fexist)
 
  If ( .Not. fexist) Then
    Print *, "Cannot find file: ", trim (file_name)
    Print *, "Setting station precipitation timeseries to missing"
 
  Else
    Open (12, File=file_name, Status='old')
 
    var_index = - 1
    i = 1
    Do
      Read (12, "(A)", IoStat=stat) line
      If (stat < 0) Exit
      line = adjustl (line)
 
      If (line(1:1) .Ne. "#" .And. len(trim(line)) /= 0) Then
 
        If (var_index ==-1) Then
          ipos = index (line, "DATE")
          If (ipos > 0) Then
            Call parse (line, " ", peices, npeices)
            If (npeices < 3 .Or. trim(peices(1)) .Ne. "DATE" .Or. trim(peices(2)) .Ne. "HHMMSS") Then
              Print *, "Failed to read header from file ", trim (file_name)
              error = 1
              Return
            End If
            Do i = 3, npeices, 1
              If (trim(peices(i)) .Eq. trim(stnvar)) Then
                var_index = i
              End If
            End Do
          End If
 
        Else
 
          Call parse (line, " ", peices, npeices)
 
          If (npeices < 3) Then
            Print *, "Confusing line from file ", trim (file_name), " line: ", trim (line)
            error = 1
            Return
          End If
 
          utime = date_to_unix (trim(peices(1))//trim(peices(2)))
          If (utime > Times(1)-43200.0 .And. utime <= Times(ntimes)+43200.0) Then
            val_index = FLOOR ((utime-Times(1))/86400.0) + 1
            vals_miss (val_index) = .True.
            Call value (peices(var_index), vals(val_index), error)
            If (error /= 0) Print *, "Problem converting string to float: ", peices (var_index)
            error = 0
            If (vals(val_index) ==-999.0) Then
              vals_miss (val_index) = .False.
            End If
          End If
        End If
 
      End If
    End Do
 
    If (var_index ==-1) Then
      Print *, "Failed to find header from file ", trim (file_name)
      error = 1
    End If
  End If !end file exist if statement
 
!!! read in station temperature data now
  error = 0
 
  stnvar_t = 'Tmin'
  ntimes = size (Times)
  Allocate (tair_vals(2, ntimes))
  Allocate (vals_miss_t(ntimes))
  vals_miss_t = .False.
  tair_vals = - 999.0
 
  Call parse (site_list, "/", peices, npeices)
 
  If (npeices == 0) Then
    directory = "."
  Else
    directory = peices (1)
    If (npeices > 2) Then
      Do i = 2, npeices - 1, 1
        directory = trim (directory) // "/" // peices (i)
      End Do
    End If
  End If
 
  file_name = trim (directory) // "/" // trim (stnid) // "_" // trim (site_var_t) // ".txt"
 
  Print *, "Reading station file: ", trim (file_name)
 
  Inquire (File=file_name, Exist=fexist)
  If ( .Not. fexist) Then
    Print *, "Cannot find file: ", trim (file_name)
    Print *, "Setting station temperature timeseries to missing"
 
  Else
    Open (12, File=file_name, Status='old')
 
    var_index = - 1
    i = 1
    Do
      Read (12, "(A)", IoStat=stat) line
      If (stat < 0) Exit
      line = adjustl (line)
 
      If (line(1:1) .Ne. "#" .And. len(trim(line)) /= 0) Then
 
        If (var_index ==-1) Then
          ipos = index (line, "DATE")
          If (ipos > 0) Then
            Call parse (line, " ", peices, npeices)
            If (npeices < 3 .Or. trim(peices(1)) .Ne. "DATE" .Or. trim(peices(2)) .Ne. "HHMMSS") Then
              Print *, "Failed to read header from file ", trim (file_name)
              error = 1
              Return
            End If
            Do i = 3, npeices, 1
              If (trim(peices(i)) .Eq. trim(stnvar_t)) Then
                var_index = i
              End If
            End Do
          End If
 
        Else
 
          Call parse (line, " ", peices, npeices)
 
          If (npeices < 3) Then
            Print *, "Confusing line from file ", trim (file_name), " line: ", trim (line)
            error = 1
            Return
          End If
 
          utime = date_to_unix (trim(peices(1))//trim(peices(2)))
          If (utime > Times(1)-43200.0 .And. utime <= Times(ntimes)+43200.0) Then
            val_index = FLOOR ((utime-Times(1))/86400.0) + 1
            vals_miss_t (val_index) = .True.
            Call value (peices(var_index), tair_vals(1, val_index), error)
            Call value (peices(var_index+1), tair_vals(2, val_index), error)
 
            If (error /= 0) Print *, "Problem converting string to float: ", peices (var_index)
            error = 0
 
            If (tair_vals(1, val_index) ==-999.0 .Or. tair_vals(2, val_index) ==-999.0) Then
              vals_miss_t (val_index) = .False.
            End If
          End If
        End If
 
      End If
    End Do
 
    If (var_index ==-1) Then
      Print *, "Failed to find header from file ", trim (file_name)
      error = 1
    End If
  End If !end file exist if statement
 
End Subroutine read_station
 
 
 
Subroutine read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
  Use strings
  Use type
  Implicit None
 
  Character (Len=500), Intent (In) :: file_name
  Real (DP), Allocatable, Intent (Out) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
  Integer (I4B), Intent (Out) :: nx, ny
  Integer, Intent (Out) :: error
 
  Character (Len=100) :: peices (5)
  Integer :: ngrid, i, npeices, stat, err, ipos
  Real (DP) :: dx, dy, startx, starty
  Character (300) :: line
  Logical fexist
 
  ngrid = 0
  nx = 0
  ny = 0
  dx = 0.0
  dy = 0.0
  startx = 0.0
  starty = 0.0
 
  error = 0
  Print *, "Reading grid file: ", trim (file_name)
 
  Inquire (File=file_name, Exist=fexist)
  If ( .Not. fexist) Then
    Print *, "Cannot find file: ", trim (file_name)
    error = 1
    Return
  End If
 
  Open (13, File=file_name, Status='old')
 
  i = 1
  Do
    Read (13, "(A)", IoStat=stat) line
    If (stat < 0) Exit
    line = adjustl (line)
    If (line(1:1) .Ne. "#" .And. len(trim(line)) /= 0) Then
 
      ipos = index (line, "NX")
      If (ipos > 0) Then
        ipos = ipos + 2
        Call value (line(ipos:), nx, err)
        If (nx > 0 .And. ny > 0) Then
          ngrid = nx * ny
          Allocate (lats(ngrid))
          Allocate (lons(ngrid))
          Allocate (alts(ngrid))
          Allocate (slp_n(ngrid))
          Allocate (slp_e(ngrid))
        End If
      End If
      ipos = index (line, "NY")
      If (ipos > 0) Then
        ipos = ipos + 2
        Call value (line(ipos:), ny, err)
        If (nx > 0 .And. ny > 0) Then
          ngrid = nx * ny
          Allocate (lats(ngrid))
          Allocate (lons(ngrid))
          Allocate (alts(ngrid))
          Allocate (slp_n(ngrid))
          Allocate (slp_e(ngrid))
        End If
      End If
      ipos = index (line, "DX")
      If (ipos > 0) Then
        ipos = ipos + 2
        Call value (line(ipos:), dx, err)
      End If
      ipos = index (line, "DY")
      If (ipos > 0) Then
        ipos = ipos + 2
        Call value (line(ipos:), dy, err)
      End If
      ipos = index (line, "STARTX")
      If (ipos > 0) Then
        ipos = ipos + 6
        Call value (line(ipos:), startx, err)
      End If
      ipos = index (line, "STARTY")
      If (ipos > 0) Then
        ipos = ipos + 6
        Call value (line(ipos:), starty, err)
      End If
    End If
 
    ipos = index (line, ',')
    If (ipos > 0 .And. ngrid > 0) Then
      Call parse (line, ",", peices, npeices)
      If (npeices == 5) Then
        Call value (peices(1), lats(i), err)
        If (err /= 0) lats (i) = - 999.99
        Call value (peices(2), lons(i), err)
        If (err /= 0) lons (i) = - 999.99
        Call value (peices(3), alts(i), err)
        If (err /= 0) alts (i) = - 999.99
        Call value (peices(4), slp_n(i), err)
        If (err /= 0) slp_n (i) = - 999.99
        Call value (peices(5), slp_e(i), err)
        If (err /= 0) slp_e (i) = - 999.99
        i = i + 1
      End If
    End If
 
  End Do
  If (ngrid == 0) Then
    Print *, "Failed to find NX, NY in grid list: ", trim (file_name)
    error = 1
  Else
    If (i /= ngrid+1) Then
      Print *, "Found only ", i, " out of ", ngrid, " points from: ", trim (file_name)
      error = 1
    End If
  End If
 
End Subroutine read_grid_list
