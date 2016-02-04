 
SUBROUTINE read_config (fname, n, names, values)
 
  USE strings
  IMPLICIT NONE
 
  CHARACTER (LEN=*) :: fname
  INTEGER :: n
  CHARACTER (LEN=500) :: names (n)
  CHARACTER (LEN=500) :: values (n)
 
 
  CHARACTER (LEN=500) :: peices (2)
  CHARACTER (1000) line
  INTEGER ipos, stat, npeices, i
  LOGICAL fexist
 
  INQUIRE (FILE=fname, EXIST=fexist)
  IF ( .NOT. fexist) THEN
   PRINT *, "Cannot find config file: ", trim (fname)
   RETURN
  END IF
 
  DO i = 1, n, 1
   values (i) = ""
  END DO
 
  OPEN (11, FILE=fname, STATUS='old', FORM='formatted')
 
  DO
   READ (11, "(A)", IOSTAT=stat) line
   IF (stat < 0) EXIT
   line = adjustl (line)
   IF (line(1:1) .NE. "!" .AND. len(trim(line)) /= 0) THEN
    ipos = index (line, '=')
    IF (ipos > 0) THEN
     CALL parse (line, "=", peices, npeices)
     DO i = 1, n, 1
      IF (index(peices(1), names(i)) == 1) values (i) = peices (2)
     END DO
    END IF
 
   END IF
 
  END DO
 
END SUBROUTINE read_config
