SUBROUTINE read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
  USE strings
  USE nrtype
  IMPLICIT NONE
 
  CHARACTER (LEN=500), INTENT (IN) :: file_name
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
  INTEGER (I4B), INTENT (OUT) :: nx, ny
  INTEGER, INTENT (OUT) :: error
 
  CHARACTER (LEN=100) :: peices (5)
  INTEGER :: ngrid, i, npeices, stat, err, ipos
  REAL (DP) :: dx, dy, startx, starty
  CHARACTER (300) :: line
  LOGICAL fexist
 
  ngrid = 0
  nx = 0
  ny = 0
  dx = 0.0
  dy = 0.0
  startx = 0.0
  starty = 0.0
 
  error = 0
  PRINT *, "Reading grid file: ", trim (file_name)
 
  INQUIRE (FILE=file_name, EXIST=fexist)
  IF ( .NOT. fexist) THEN
   PRINT *, "Cannot find file: ", trim (file_name)
   error = 1
   RETURN
  END IF
 
  OPEN (13, FILE=file_name, STATUS='old')
 
  i = 1
  DO
   READ (13, "(A)", IOSTAT=stat) line
   IF (stat < 0) EXIT
   line = adjustl (line)
   IF (line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) THEN
 
    ipos = index (line, "NX")
    IF (ipos > 0) THEN
     ipos = ipos + 2
     CALL value (line(ipos:), nx, err)
     IF (nx > 0 .AND. ny > 0) THEN
      ngrid = nx * ny
      ALLOCATE (lats(ngrid))
      ALLOCATE (lons(ngrid))
      ALLOCATE (alts(ngrid))
      ALLOCATE (slp_n(ngrid))
      ALLOCATE (slp_e(ngrid))
     END IF
    END IF
    ipos = index (line, "NY")
    IF (ipos > 0) THEN
     ipos = ipos + 2
     CALL value (line(ipos:), ny, err)
     IF (nx > 0 .AND. ny > 0) THEN
      ngrid = nx * ny
      ALLOCATE (lats(ngrid))
      ALLOCATE (lons(ngrid))
      ALLOCATE (alts(ngrid))
      ALLOCATE (slp_n(ngrid))
      ALLOCATE (slp_e(ngrid))
     END IF
    END IF
    ipos = index (line, "DX")
    IF (ipos > 0) THEN
     ipos = ipos + 2
     CALL value (line(ipos:), dx, err)
    END IF
    ipos = index (line, "DY")
    IF (ipos > 0) THEN
     ipos = ipos + 2
     CALL value (line(ipos:), dy, err)
    END IF
    ipos = index (line, "STARTX")
    IF (ipos > 0) THEN
     ipos = ipos + 6
     CALL value (line(ipos:), startx, err)
    END IF
    ipos = index (line, "STARTY")
    IF (ipos > 0) THEN
     ipos = ipos + 6
     CALL value (line(ipos:), starty, err)
    END IF
   END IF
 
   ipos = index (line, ',')
   IF (ipos > 0 .AND. ngrid > 0) THEN
    CALL parse (line, ",", peices, npeices)
    IF (npeices == 5) THEN
     CALL value (peices(1), lats(i), err)
     IF (err /= 0) lats (i) = - 999.99
     CALL value (peices(2), lons(i), err)
     IF (err /= 0) lons (i) = - 999.99
     CALL value (peices(3), alts(i), err)
     IF (err /= 0) alts (i) = - 999.99
     CALL value (peices(4), slp_n(i), err)
     IF (err /= 0) slp_n (i) = - 999.99
     CALL value (peices(5), slp_e(i), err)
     IF (err /= 0) slp_e (i) = - 999.99
     i = i + 1
    END IF
   END IF
 
  END DO
  IF (ngrid == 0) THEN
   PRINT *, "Failed to find NX, NY in grid list: ", trim (file_name)
   error = 1
  ELSE
   IF (i /= ngrid+1) THEN
    PRINT *, "Found only ", i, " out of ", ngrid, " points from: ", trim (file_name)
    error = 1
   END IF
  END IF
 
END SUBROUTINE read_grid_list
