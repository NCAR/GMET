Subroutine read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
  Use strings
  Use nrtype
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
