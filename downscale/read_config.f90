 
Subroutine read_config (fname, n, names, values)
 
  Use strings
  Implicit None
 
  Character (Len=*) :: fname
  Integer :: n
  Character (Len=500) :: names (n)
  Character (Len=500) :: values (n)
 
 
  Character (Len=500) :: peices (2)
  Character (1000) line
  Integer ipos, stat, npeices, i
  Logical fexist
 
  Inquire (File=fname, Exist=fexist)
  If ( .Not. fexist) Then
    Print *, "Cannot find config file: ", trim (fname)
    Return
  End If
 
  Do i = 1, n, 1
    values (i) = ""
  End Do
 
  Open (11, File=fname, Status='old', Form='formatted')
 
  Do
    Read (11, "(A)", IoStat=stat) line
    If (stat < 0) Exit
    line = adjustl (line)
    If (line(1:1) .Ne. "!" .And. len(trim(line)) /= 0) Then
      ipos = index (line, '=')
      If (ipos > 0) Then
        Call parse (line, "=", peices, npeices)
        Do i = 1, n, 1
          If (index(peices(1), names(i)) == 1) values (i) = peices (2)
        End Do
      End If
 
    End If
 
  End Do
 
End Subroutine read_config
