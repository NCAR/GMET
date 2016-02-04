 
SUBROUTINE read_refcst (startdate, enddate, file_var, perturbation, var_name, forecast, V, X, Y, T, error)
  USE strings
  USE utim
  USE type
  IMPLICIT NONE
 
  INTERFACE
   SUBROUTINE read_nc_file (file_name, valid_time, var_name, var, lats, lons, error)
    USE type
    CHARACTER (LEN=*), INTENT (IN) :: file_name
    CHARACTER (LEN=*), INTENT (IN) :: var_name
    REAL (DP), INTENT (IN) :: valid_time
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: var (:, :, :)
    REAL (DP), ALLOCATABLE, INTENT (OUT) :: lats (:), lons (:)
    INTEGER, INTENT (OUT) :: error
   END SUBROUTINE read_nc_file
  END INTERFACE
 
  CHARACTER (LEN=100), INTENT (IN) :: startdate, enddate, file_var, perturbation
  CHARACTER (LEN=*), INTENT (IN) :: var_name
  INTEGER (I4B), INTENT (IN) :: forecast
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: V (:, :), X (:), Y (:)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: T (:)
  INTEGER, INTENT (OUT) :: error
 
  REAL (DP), ALLOCATABLE :: lats (:), lons (:), var (:, :, :)
  REAL (DP) :: valid_time, utime
  CHARACTER (LEN=500) :: file
  CHARACTER (LEN=100) :: date
  INTEGER (I4B) :: vshape (3)
  INTEGER (I4B) :: sday, eday, ntimes, ngrids
  INTEGER (I4B) :: sec, Min, hour, day, month, year
  INTEGER (I4B) :: i, j, k
 
  error = 0
 
  CALL parse_date (startdate, year, month, day, hour, Min, sec, error)
  sday = julian_date (day, month, year)
  CALL parse_date (enddate, year, month, day, hour, Min, sec, error)
  eday = julian_date (day, month, year)
  ntimes = eday - sday + 1
 
  ALLOCATE (T(ntimes))
 
  i = 1
  utime = date_to_unix (startdate)
  DO
   IF (utime > date_to_unix(enddate)) EXIT
   CALL unix_to_date (utime, year, month, day, hour, Min, sec)
 
   WRITE (date, "(I4.4I2.2I2.2)"), year, month, day
   IF (forecast >= 190) THEN
    file = trim (date) // "/" // trim (file_var) // "_" // trim (date) // "00_" // trim (perturbation) // "_t190.nc"
   ELSE
    file = trim (date) // "/" // trim (file_var) // "_" // trim (date) // "00_" // trim (perturbation) // ".nc"
   END IF
 
 
   valid_time = utime + (forecast*3600.0)
 
   CALL read_nc_file (file, valid_time, var_name, var, lats, lons, error)
 
   IF (trim(var_name) == 'APCP_ens_mean_surface') THEN
    var = var ** (1.0/4.0)
    PRINT *, 'normalizing GEFS precip'
   END IF
 
   IF (error == 1) THEN
    PRINT *, "Failed reading file ", trim (file)
    PRINT *, 'var ', trim (var_name)
    EXIT
   ELSE
    PRINT ("(AAAF11.0)"), "Success reading file ", trim (file), " Valid at ", valid_time
   END IF
 
   vshape = shape (var)
   IF (vshape(1) /= size(lons) .OR. vshape(2) /= size(lats) .OR. vshape(3) /= 1) THEN
    PRINT *, "Dimensions from file do not match"
    error = 1
    EXIT
   END IF
 
   IF (i == 1) THEN
    ngrids = vshape (1) * vshape (2)
    ALLOCATE (V(ngrids, ntimes))
    ALLOCATE (X(vshape(2)))
    ALLOCATE (Y(vshape(1)))
    X (:) = lats (:)
    DO j = 1, vshape (1), 1
     Y (j) = Mod (lons(j), 360.0) - 360.0
    END DO
 
   ELSE
    IF (vshape(2) /= size(X) .OR. vshape(1) /= size(Y)) THEN
     PRINT *, "Dimensions from file do not match previous files"
     error = 1
     EXIT
    END IF
   END IF
 
   T (i) = valid_time
   V (1:ngrids, i) = reshape (var, (/ ngrids /))
 
   i = i + 1
   utime = utime + 86400
 
  END DO
 
END SUBROUTINE read_refcst
 
 
SUBROUTINE read_station_list (file_name, id, name, lat, lon, alt, sslp_n, sslp_e, n_stations, error)
  USE strings
  USE type
  IMPLICIT NONE
 
  CHARACTER (LEN=500), INTENT (IN) :: file_name
  CHARACTER (LEN=100), ALLOCATABLE, INTENT (OUT) :: id (:), name (:)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: lat (:), lon (:), alt (:), sslp_n (:), sslp_e (:)
  INTEGER (I4B), INTENT (OUT) :: n_stations
  INTEGER, INTENT (OUT) :: error
 
  CHARACTER (LEN=100) :: peices (7)
  INTEGER (I4B) i, npeices, stat, err, ipos
  CHARACTER (500) line
  LOGICAL fexist
 
  error = 0
  PRINT *, "Reading stations list: ", trim (file_name)
 
  INQUIRE (FILE=file_name, EXIST=fexist)
  IF ( .NOT. fexist) THEN
   PRINT *, "Cannot find file: ", trim (file_name)
   error = 1
   RETURN
  END IF
 
  OPEN (11, FILE=file_name, STATUS='old')
 
  n_stations = 0
  i = 1
  DO
   READ (11, "(A)", IOSTAT=stat) line
   IF (stat < 0) EXIT
   line = adjustl (line)
 
   IF (line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) THEN
    ipos = index (line, "NSITES")
    IF (ipos > 0) THEN
     ipos = ipos + 6
     CALL value (line(ipos:), n_stations, err)
     IF (err /= 0) THEN
      n_stations = 0
     ELSE
      PRINT *, "Stations: ", n_stations
      ALLOCATE (id(n_stations))
      ALLOCATE (name(n_stations))
      ALLOCATE (lat(n_stations))
      ALLOCATE (lon(n_stations))
      ALLOCATE (alt(n_stations))
      ALLOCATE (sslp_n(n_stations))
      ALLOCATE (sslp_e(n_stations))
     END IF
    END IF
    ipos = index (line, ',')
    IF (ipos > 0 .AND. n_stations > 0) THEN
     CALL parse (line, ",", peices, npeices)
     IF (npeices == 7) THEN
      id (i) = peices (1)
      name (i) = peices (7)
      CALL delall (name(i), '"')
      CALL value (peices(2), lat(i), err)
      IF (err /= 0) lat (i) = - 999.99
      CALL value (peices(3), lon(i), err)
      IF (err /= 0) lon (i) = - 999.99
      CALL value (peices(4), alt(i), err)
      IF (err /= 0) alt (i) = - 999.99
      CALL value (peices(5), sslp_n(i), err)
      IF (err /= 0) sslp_n (i) = - 999.99
      CALL value (peices(6), sslp_e(i), err)
      IF (err /= 0) sslp_e (i) = - 999.99
 
      i = i + 1
     END IF
    END IF
 
   END IF
 
  END DO
  IF (n_stations == 0) THEN
   PRINT *, "Failed to find NSITES in station list: ", trim (file_name)
   error = 1
  ELSE
   IF (i /= n_stations+1) THEN
    PRINT *, "Found only ", i, " out of ", n_stations, " stations from: ", trim (file_name)
    error = 1
   END IF
  END IF
 
END SUBROUTINE read_station_list
 
 
SUBROUTINE read_station (stnvar, stnid, site_var, site_var_t, site_list, Times, vals, tair_vals, vals_miss, &
& vals_miss_t, error)
  USE strings
  USE utim
  USE type
  IMPLICIT NONE
 
  CHARACTER (LEN=100), INTENT (IN) :: stnvar
  CHARACTER (LEN=100), INTENT (IN) :: stnid
  CHARACTER (LEN=100), INTENT (IN) :: site_var, site_var_t
  CHARACTER (LEN=500), INTENT (IN) :: site_list
  REAL (DP), INTENT (IN) :: Times (:)
  REAL (DP), ALLOCATABLE, INTENT (OUT) :: vals (:), tair_vals (:, :)
  LOGICAL, ALLOCATABLE, INTENT (OUT) :: vals_miss (:), vals_miss_t (:)
  INTEGER, INTENT (OUT) :: error
 
  CHARACTER (LEN=500) :: file_name
  CHARACTER (LEN=500) :: directory
  CHARACTER (LEN=500) :: line
  CHARACTER (LEN=100) :: peices (10)
  CHARACTER (LEN=100) :: stnvar_t
  REAL (DP) :: utime
  INTEGER (I4B) :: i, npeices, ntimes, ipos, stat, var_index, val_index
  INTEGER (I4B) :: sec, Min, hour, day, month, year
  LOGICAL :: fexist
 
 
!!!! first read in station precipitation data
  error = 0
 
  CALL parse (site_list, "/", peices, npeices)
 
  IF (npeices == 0) THEN
   directory = "."
  ELSE
   directory = peices (1)
   IF (npeices > 2) THEN
    DO i = 2, npeices - 1, 1
     directory = trim (directory) // "/" // peices (i)
    END DO
   END IF
  END IF
  file_name = trim (directory) // "/" // trim (stnid) // "_" // trim (site_var) // ".txt"
  PRINT *, "Reading station file: ", trim (file_name)
 
  ntimes = size (Times)
  ALLOCATE (vals(ntimes))
  ALLOCATE (vals_miss(ntimes))
  vals_miss = .FALSE.
  vals = - 999.0
 
  INQUIRE (FILE=file_name, EXIST=fexist)
 
  IF ( .NOT. fexist) THEN
   PRINT *, "Cannot find file: ", trim (file_name)
   PRINT *, "Setting station precipitation timeseries to missing"
 
  ELSE
   OPEN (12, FILE=file_name, STATUS='old')
 
   var_index = - 1
   i = 1
   DO
    READ (12, "(A)", IOSTAT=stat) line
    IF (stat < 0) EXIT
    line = adjustl (line)
 
    IF (line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) THEN
 
     IF (var_index ==-1) THEN
      ipos = index (line, "DATE")
      IF (ipos > 0) THEN
       CALL parse (line, " ", peices, npeices)
       IF (npeices < 3 .OR. trim(peices(1)) .NE. "DATE" .OR. trim(peices(2)) .NE. "HHMMSS") THEN
        PRINT *, "Failed to read header from file ", trim (file_name)
        error = 1
        RETURN
       END IF
       DO i = 3, npeices, 1
        IF (trim(peices(i)) .EQ. trim(stnvar)) THEN
         var_index = i
        END IF
       END DO
      END IF
 
     ELSE
 
      CALL parse (line, " ", peices, npeices)
 
      IF (npeices < 3) THEN
       PRINT *, "Confusing line from file ", trim (file_name), " line: ", trim (line)
       error = 1
       RETURN
      END IF
 
      utime = date_to_unix (trim(peices(1))//trim(peices(2)))
      IF (utime > Times(1)-43200.0 .AND. utime <= Times(ntimes)+43200.0) THEN
       val_index = FLOOR ((utime-Times(1))/86400.0) + 1
       vals_miss (val_index) = .TRUE.
       CALL value (peices(var_index), vals(val_index), error)
       IF (error /= 0) PRINT *, "Problem converting string to float: ", peices (var_index)
       error = 0
       IF (vals(val_index) ==-999.0) THEN
        vals_miss (val_index) = .FALSE.
       END IF
      END IF
     END IF
 
    END IF
   END DO
 
   IF (var_index ==-1) THEN
    PRINT *, "Failed to find header from file ", trim (file_name)
    error = 1
   END IF
  END IF !end file exist if statement
 
!!! read in station temperature data now
  error = 0
 
  stnvar_t = 'Tmin'
  ntimes = size (Times)
  ALLOCATE (tair_vals(2, ntimes))
  ALLOCATE (vals_miss_t(ntimes))
  vals_miss_t = .FALSE.
  tair_vals = - 999.0
 
  CALL parse (site_list, "/", peices, npeices)
 
  IF (npeices == 0) THEN
   directory = "."
  ELSE
   directory = peices (1)
   IF (npeices > 2) THEN
    DO i = 2, npeices - 1, 1
     directory = trim (directory) // "/" // peices (i)
    END DO
   END IF
  END IF
 
  file_name = trim (directory) // "/" // trim (stnid) // "_" // trim (site_var_t) // ".txt"
 
  PRINT *, "Reading station file: ", trim (file_name)
 
  INQUIRE (FILE=file_name, EXIST=fexist)
  IF ( .NOT. fexist) THEN
   PRINT *, "Cannot find file: ", trim (file_name)
   PRINT *, "Setting station temperature timeseries to missing"
 
  ELSE
   OPEN (12, FILE=file_name, STATUS='old')
 
   var_index = - 1
   i = 1
   DO
    READ (12, "(A)", IOSTAT=stat) line
    IF (stat < 0) EXIT
    line = adjustl (line)
 
    IF (line(1:1) .NE. "#" .AND. len(trim(line)) /= 0) THEN
 
     IF (var_index ==-1) THEN
      ipos = index (line, "DATE")
      IF (ipos > 0) THEN
       CALL parse (line, " ", peices, npeices)
       IF (npeices < 3 .OR. trim(peices(1)) .NE. "DATE" .OR. trim(peices(2)) .NE. "HHMMSS") THEN
        PRINT *, "Failed to read header from file ", trim (file_name)
        error = 1
        RETURN
       END IF
       DO i = 3, npeices, 1
        IF (trim(peices(i)) .EQ. trim(stnvar_t)) THEN
         var_index = i
        END IF
       END DO
      END IF
 
     ELSE
 
      CALL parse (line, " ", peices, npeices)
 
      IF (npeices < 3) THEN
       PRINT *, "Confusing line from file ", trim (file_name), " line: ", trim (line)
       error = 1
       RETURN
      END IF
 
      utime = date_to_unix (trim(peices(1))//trim(peices(2)))
      IF (utime > Times(1)-43200.0 .AND. utime <= Times(ntimes)+43200.0) THEN
       val_index = FLOOR ((utime-Times(1))/86400.0) + 1
       vals_miss_t (val_index) = .TRUE.
       CALL value (peices(var_index), tair_vals(1, val_index), error)
       CALL value (peices(var_index+1), tair_vals(2, val_index), error)
 
       IF (error /= 0) PRINT *, "Problem converting string to float: ", peices (var_index)
       error = 0
 
       IF (tair_vals(1, val_index) ==-999.0 .OR. tair_vals(2, val_index) ==-999.0) THEN
        vals_miss_t (val_index) = .FALSE.
       END IF
      END IF
     END IF
 
    END IF
   END DO
 
   IF (var_index ==-1) THEN
    PRINT *, "Failed to find header from file ", trim (file_name)
    error = 1
   END IF
  END IF !end file exist if statement
 
END SUBROUTINE read_station
 
 
 
SUBROUTINE read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
  USE strings
  USE type
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
