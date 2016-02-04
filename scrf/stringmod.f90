MODULE strings
 
  USE precision
 
  PRIVATE :: value_dr, value_sr, value_di, value_si
  PRIVATE :: write_dr, write_sr, write_di, write_si
  PRIVATE :: writeq_dr, writeq_sr, writeq_di, writeq_si
 
  INTERFACE value ! Generic operator for converting a number string to a
                 ! number. Calling syntax is 'call value(numstring,number,ios)'
                 ! where 'numstring' is a number string and 'number' is a
                 ! real number or an integer (single or double precision).
   MODULE PROCEDURE value_dr
   MODULE PROCEDURE value_sr
   MODULE PROCEDURE value_di
   MODULE PROCEDURE value_si
  END INTERFACE
 
  INTERFACE writenum ! Generic  interface for writing a number to a string. The
                    ! number is left justified in the string. The calling syntax
                    ! is 'call writenum(number,string,format)' where 'number' is
                    ! a real number or an integer, 'string' is a character string
                    ! containing the result, and 'format' is the format desired,
                    ! e.g., 'e15.6' or 'i5'.
   MODULE PROCEDURE write_dr
   MODULE PROCEDURE write_sr
   MODULE PROCEDURE write_di
   MODULE PROCEDURE write_si
  END INTERFACE
 
  INTERFACE writeq ! Generic interface equating a name to a numerical value. The
                  ! calling syntax is 'call writeq(unit,name,value,format)' where
                  ! unit is the integer output unit number, 'name' is the variable
                  ! name, 'value' is the real or integer value of the variable,
                  ! and 'format' is the format of the value. The result written to
                  ! the output unit has the form <name> = <value>.
   MODULE PROCEDURE writeq_dr
   MODULE PROCEDURE writeq_sr
   MODULE PROCEDURE writeq_di
   MODULE PROCEDURE writeq_si
  END INTERFACE
 
 
!**********************************************************************
 
CONTAINS
 
!**********************************************************************
 
  SUBROUTINE parse (str, delims, args, nargs)
 
! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.
 
   CHARACTER (LEN=*) :: str, delims
   CHARACTER (LEN=len_trim(str)) :: strsav
   CHARACTER (LEN=*), DIMENSION (:) :: args
 
   strsav = str
   CALL compact (str)
   na = size (args)
   DO i = 1, na
    args (i) = ' '
   END DO
   nargs = 0
   lenstr = len_trim (str)
   IF (lenstr == 0) RETURN
   k = 0
 
   DO
    IF (len_trim(str) == 0) EXIT
    nargs = nargs + 1
    CALL split (str, delims, args(nargs))
    CALL removebksl (args(nargs))
   END DO
   str = strsav
 
  END SUBROUTINE parse
 
!**********************************************************************
 
  SUBROUTINE compact (str)
 
! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.
 
   CHARACTER (LEN=*) :: str
   CHARACTER (LEN=1) :: ch
   CHARACTER (LEN=len_trim(str)) :: outstr
 
   str = adjustl (str)
   lenstr = len_trim (str)
   outstr = ' '
   isp = 0
   k = 0
 
   DO i = 1, lenstr
    ch = str (i:i)
    ich = iachar (ch)
 
    SELECT CASE (ich)
 
    CASE (9, 32)! space or tab character
     IF (isp == 0) THEN
      k = k + 1
      outstr (k:k) = ' '
     END IF
     isp = 1
 
    CASE (33:)! not a space, quote, or control character
     k = k + 1
     outstr (k:k) = ch
     isp = 0
 
    END SELECT
 
   END DO
 
   str = adjustl (outstr)
 
  END SUBROUTINE compact
 
!**********************************************************************
 
  SUBROUTINE removesp (str)
 
! Removes spaces, tabs, and control characters in string str
 
   CHARACTER (LEN=*) :: str
   CHARACTER (LEN=1) :: ch
   CHARACTER (LEN=len_trim(str)) :: outstr
 
   str = adjustl (str)
   lenstr = len_trim (str)
   outstr = ' '
   k = 0
 
   DO i = 1, lenstr
    ch = str (i:i)
    ich = iachar (ch)
    SELECT CASE (ich)
    CASE (0:32)! space, tab, or control character
     CYCLE
    CASE (33:)
     k = k + 1
     outstr (k:k) = ch
    END SELECT
   END DO
 
   str = adjustl (outstr)
 
  END SUBROUTINE removesp
 
!**********************************************************************
 
  SUBROUTINE value_dr (str, rnum, ios)
 
! Converts number string to a double precision real number
 
   CHARACTER (LEN=*) :: str
   REAL (kr8) :: rnum
   INTEGER :: ios
 
   ilen = len_trim (str)
   ipos = scan (str, 'Ee')
   IF ( .NOT. is_digit(str(ilen:ilen)) .AND. ipos /= 0) THEN
    ios = 3
    RETURN
   END IF
   READ (str,*, IOSTAT=ios) rnum
 
  END SUBROUTINE value_dr
 
!**********************************************************************
 
  SUBROUTINE value_sr (str, rnum, ios)
 
! Converts number string to a single precision real number
 
   CHARACTER (LEN=*) :: str
   REAL (kr4) :: rnum
   REAL (kr8) :: rnumd
 
   CALL value_dr (str, rnumd, ios)
   IF (Abs(rnumd) > huge(rnum)) THEN
    ios = 15
    RETURN
   END IF
   IF (Abs(rnumd) < tiny(rnum)) rnum = 0.0_kr4
   rnum = rnumd
 
  END SUBROUTINE value_sr
 
!**********************************************************************
 
  SUBROUTINE value_di (str, inum, ios)
 
! Converts number string to a double precision integer value
 
   CHARACTER (LEN=*) :: str
   INTEGER (ki8) :: inum
   REAL (kr8) :: rnum
 
   CALL value_dr (str, rnum, ios)
   IF (Abs(rnum) > huge(inum)) THEN
    ios = 15
    RETURN
   END IF
   inum = Nint (rnum, ki8)
 
  END SUBROUTINE value_di
 
!**********************************************************************
 
  SUBROUTINE value_si (str, inum, ios)
 
! Converts number string to a single precision integer value
 
   CHARACTER (LEN=*) :: str
   INTEGER (ki4) :: inum
   REAL (kr8) :: rnum
 
   CALL value_dr (str, rnum, ios)
   IF (Abs(rnum) > huge(inum)) THEN
    ios = 15
    RETURN
   END IF
   inum = Nint (rnum, ki4)
 
  END SUBROUTINE value_si
 
!**********************************************************************
 
  SUBROUTINE shiftstr (str, n)
 
! Shifts characters in in the string 'str' n positions (positive values
! denote a right shift and negative values denote a left shift). Characters
! that are shifted off the end are lost. Positions opened up by the shift
! are replaced by spaces.
 
   CHARACTER (LEN=*) :: str
 
   lenstr = len (str)
   nabs = iabs (n)
   IF (nabs >= lenstr) THEN
    str = repeat (' ', lenstr)
    RETURN
   END IF
   IF (n < 0) str = str (nabs+1:) // repeat (' ', nabs)! shift left
   IF (n > 0) str = repeat (' ', nabs) // str (:lenstr-nabs)! shift right
   RETURN
 
  END SUBROUTINE shiftstr
 
!**********************************************************************
 
  SUBROUTINE insertstr (str, strins, loc)
 
! Inserts the string 'strins' into the string 'str' at position 'loc'.
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. Trailing spaces of 'strins' are
! removed prior to insertion
 
   CHARACTER (LEN=*) :: str, strins
   CHARACTER (LEN=LEN(str)) :: tempstr
 
   lenstrins = len_trim (strins)
   tempstr = str (loc:)
   CALL shiftstr (tempstr, lenstrins)
   tempstr (1:lenstrins) = strins (1:lenstrins)
   str (loc:) = tempstr
   RETURN
 
  END SUBROUTINE insertstr
 
!**********************************************************************
 
  SUBROUTINE delsubstr (str, substr)
 
! Deletes first occurrence of substring 'substr' from string 'str' and
! shifts characters left to fill hole. Trailing spaces or blanks are
! not considered part of 'substr'.
 
   CHARACTER (LEN=*) :: str, substr
 
   lensubstr = len_trim (substr)
   ipos = index (str, substr)
   IF (ipos == 0) RETURN
   IF (ipos == 1) THEN
    str = str (lensubstr+1:)
   ELSE
    str = str (:ipos-1) // str (ipos+lensubstr:)
   END IF
   RETURN
 
  END SUBROUTINE delsubstr
 
!**********************************************************************
 
  SUBROUTINE delall (str, substr)
 
! Deletes all occurrences of substring 'substr' from string 'str' and
! shifts characters left to fill holes.
 
   CHARACTER (LEN=*) :: str, substr
 
   lensubstr = len_trim (substr)
   DO
    ipos = index (str, substr)
    IF (ipos == 0) EXIT
    IF (ipos == 1) THEN
     str = str (lensubstr+1:)
    ELSE
     str = str (:ipos-1) // str (ipos+lensubstr:)
    END IF
   END DO
   RETURN
 
  END SUBROUTINE delall
 
!**********************************************************************
 
  FUNCTION uppercase (str) RESULT (ucstr)
 
! convert string to upper case
 
   CHARACTER (LEN=*) :: str
   CHARACTER (LEN=len_trim(str)) :: ucstr
 
   ilen = len_trim (str)
   ioffset = iachar ('A') - iachar ('a')
   iquote = 0
   ucstr = str
   DO i = 1, ilen
    iav = iachar (str(i:i))
    IF (iquote == 0 .AND. (iav == 34 .OR. iav == 39)) THEN
     iquote = 1
     iqc = iav
     CYCLE
    END IF
    IF (iquote == 1 .AND. iav == iqc) THEN
     iquote = 0
     CYCLE
    END IF
    IF (iquote == 1) CYCLE
    IF (iav >= iachar('a') .AND. iav <= iachar('z')) THEN
     ucstr (i:i) = achar (iav+ioffset)
    ELSE
     ucstr (i:i) = str (i:i)
    END IF
   END DO
   RETURN
 
  END FUNCTION uppercase
 
!**********************************************************************
 
  FUNCTION lowercase (str) RESULT (lcstr)
 
! convert string to lower case
 
   CHARACTER (LEN=*) :: str
   CHARACTER (LEN=len_trim(str)) :: lcstr
 
   ilen = len_trim (str)
   ioffset = iachar ('A') - iachar ('a')
   iquote = 0
   lcstr = str
   DO i = 1, ilen
    iav = iachar (str(i:i))
    IF (iquote == 0 .AND. (iav == 34 .OR. iav == 39)) THEN
     iquote = 1
     iqc = iav
     CYCLE
    END IF
    IF (iquote == 1 .AND. iav == iqc) THEN
     iquote = 0
     CYCLE
    END IF
    IF (iquote == 1) CYCLE
    IF (iav >= iachar('A') .AND. iav <= iachar('Z')) THEN
     lcstr (i:i) = achar (iav-ioffset)
    ELSE
     lcstr (i:i) = str (i:i)
    END IF
   END DO
   RETURN
 
  END FUNCTION lowercase
 
!**********************************************************************
 
  SUBROUTINE readline (nunitr, line, ios)
 
! Reads line from unit=nunitr, ignoring blank lines
! and deleting comments beginning with an exclamation point(!)
 
   CHARACTER (LEN=*) :: line
 
   DO
    READ (nunitr, '(a)', IOSTAT=ios) line ! read input line
    IF (ios /= 0) RETURN
    line = adjustl (line)
    ipos = index (line, '!')
    IF (ipos == 1) CYCLE
    IF (ipos /= 0) line = line (:ipos-1)
    IF (len_trim(line) /= 0) EXIT
   END DO
   RETURN
 
  END SUBROUTINE readline
 
!**********************************************************************
 
  SUBROUTINE match (str, ipos, imatch)
 
! Sets imatch to the position in string of the delimiter matching the delimiter
! in position ipos. Allowable delimiters are (), [], {}, <>.
 
   CHARACTER (LEN=*) :: str
   CHARACTER :: delim1, delim2, ch
 
   lenstr = len_trim (str)
   delim1 = str (ipos:ipos)
   SELECT CASE (delim1)
   CASE ('(')
    idelim2 = iachar (delim1) + 1
    istart = ipos + 1
    iend = lenstr
    inc = 1
   CASE (')')
    idelim2 = iachar (delim1) - 1
    istart = ipos - 1
    iend = 1
    inc = - 1
   CASE ('[', '{', '<')
    idelim2 = iachar (delim1) + 2
    istart = ipos + 1
    iend = lenstr
    inc = 1
   CASE (']', '}', '>')
    idelim2 = iachar (delim1) - 2
    istart = ipos - 1
    iend = 1
    inc = - 1
   CASE DEFAULT
    WRITE (*,*) delim1, ' is not a valid delimiter'
    RETURN
   END SELECT
   IF (istart < 1 .OR. istart > lenstr) THEN
    WRITE (*,*) delim1, ' has no matching delimiter'
    RETURN
   END IF
   delim2 = achar (idelim2)! matching delimiter
 
   isum = 1
   DO i = istart, iend, inc
    ch = str (i:i)
    IF (ch /= delim1 .AND. ch /= delim2) CYCLE
    IF (ch == delim1) isum = isum + 1
    IF (ch == delim2) isum = isum - 1
    IF (isum == 0) EXIT
   END DO
   IF (isum /= 0) THEN
    WRITE (*,*) delim1, ' has no matching delimiter'
    RETURN
   END IF
   imatch = i
 
   RETURN
 
  END SUBROUTINE match
 
!**********************************************************************
 
  SUBROUTINE write_dr (rnum, str, fmt)
 
! Writes double precision real number rnum to string str using format fmt
 
   REAL (kr8) :: rnum
   CHARACTER (LEN=*) :: str, fmt
   CHARACTER (LEN=80) :: formt
 
   formt = '(' // trim (fmt) // ')'
   WRITE (str, formt) rnum
   str = adjustl (str)
 
  END SUBROUTINE write_dr
 
!***********************************************************************
 
  SUBROUTINE write_sr (rnum, str, fmt)
 
! Writes single precision real number rnum to string str using format fmt
 
   REAL (kr4) :: rnum
   CHARACTER (LEN=*) :: str, fmt
   CHARACTER (LEN=80) :: formt
 
   formt = '(' // trim (fmt) // ')'
   WRITE (str, formt) rnum
   str = adjustl (str)
 
  END SUBROUTINE write_sr
 
!***********************************************************************
 
  SUBROUTINE write_di (inum, str, fmt)
 
! Writes double precision integer inum to string str using format fmt
 
   INTEGER (ki8) :: inum
   CHARACTER (LEN=*) :: str, fmt
   CHARACTER (LEN=80) :: formt
 
   formt = '(' // trim (fmt) // ')'
   WRITE (str, formt) inum
   str = adjustl (str)
 
  END SUBROUTINE write_di
 
!***********************************************************************
 
  SUBROUTINE write_si (inum, str, fmt)
 
! Writes single precision integer inum to string str using format fmt
 
   INTEGER (ki4) :: inum
   CHARACTER (LEN=*) :: str, fmt
   CHARACTER (LEN=80) :: formt
 
   formt = '(' // trim (fmt) // ')'
   WRITE (str, formt) inum
   str = adjustl (str)
 
  END SUBROUTINE write_si
 
!***********************************************************************
 
  SUBROUTINE trimzero (str)
 
! Deletes nonsignificant trailing zeroes from number string str. If number
! string ends in a decimal point, one trailing zero is added.
 
   CHARACTER (LEN=*) :: str
   CHARACTER :: ch
   CHARACTER (LEN=10) :: Exp
 
   ipos = scan (str, 'eE')
   IF (ipos > 0) THEN
    Exp = str (ipos:)
    str = str (1:ipos-1)
   END IF
   lstr = len_trim (str)
   DO i = lstr, 1, - 1
    ch = str (i:i)
    IF (ch == '0') CYCLE
    IF (ch == '.') THEN
     str = str (1:i) // '0'
     IF (ipos > 0) str = trim (str) // trim (Exp)
     EXIT
    END IF
    str = str (1:i)
    EXIT
   END DO
   IF (ipos > 0) str = trim (str) // trim (Exp)
 
  END SUBROUTINE trimzero
 
!**********************************************************************
 
  SUBROUTINE writeq_dr (unit, namestr, value, fmt)
 
! Writes a string of the form <name> = value to unit
 
   REAL (kr8) :: value
   INTEGER :: unit
   CHARACTER (LEN=*) :: namestr, fmt
   CHARACTER (LEN=32) :: tempstr
 
   CALL writenum (value, tempstr, fmt)
   CALL trimzero (tempstr)
   WRITE (UNIT,*) trim (namestr) // ' = ' // trim (tempstr)
 
  END SUBROUTINE writeq_dr
 
!**********************************************************************
 
  SUBROUTINE writeq_sr (unit, namestr, value, fmt)
 
! Writes a string of the form <name> = value to unit
 
   REAL (kr4) :: value
   INTEGER :: unit
   CHARACTER (LEN=*) :: namestr, fmt
   CHARACTER (LEN=32) :: tempstr
 
   CALL writenum (value, tempstr, fmt)
   CALL trimzero (tempstr)
   WRITE (UNIT,*) trim (namestr) // ' = ' // trim (tempstr)
 
  END SUBROUTINE writeq_sr
 
!**********************************************************************
 
  SUBROUTINE writeq_di (unit, namestr, ivalue, fmt)
 
! Writes a string of the form <name> = ivalue to unit
 
   INTEGER (ki8) :: ivalue
   INTEGER :: unit
   CHARACTER (LEN=*) :: namestr, fmt
   CHARACTER (LEN=32) :: tempstr
   CALL writenum (ivalue, tempstr, fmt)
   CALL trimzero (tempstr)
   WRITE (UNIT,*) trim (namestr) // ' = ' // trim (tempstr)
 
  END SUBROUTINE writeq_di
 
!**********************************************************************
 
  SUBROUTINE writeq_si (unit, namestr, ivalue, fmt)
 
! Writes a string of the form <name> = ivalue to unit
 
   INTEGER (ki4) :: ivalue
   INTEGER :: unit
   CHARACTER (LEN=*) :: namestr, fmt
   CHARACTER (LEN=32) :: tempstr
   CALL writenum (ivalue, tempstr, fmt)
   CALL trimzero (tempstr)
   WRITE (UNIT,*) trim (namestr) // ' = ' // trim (tempstr)
 
  END SUBROUTINE writeq_si
 
!**********************************************************************
 
  FUNCTION is_letter (ch) RESULT (res)
 
! Returns .true. if ch is a letter and .false. otherwise
 
   CHARACTER :: ch
   LOGICAL :: res
 
   SELECT CASE (ch)
   CASE ('A':'Z', 'a':'z')
    res = .TRUE.
   CASE DEFAULT
    res = .FALSE.
   END SELECT
   RETURN
 
  END FUNCTION is_letter
 
!**********************************************************************
 
  FUNCTION is_digit (ch) RESULT (res)
 
! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise
 
   CHARACTER :: ch
   LOGICAL :: res
 
   SELECT CASE (ch)
   CASE ('0':'9')
    res = .TRUE.
   CASE DEFAULT
    res = .FALSE.
   END SELECT
   RETURN
 
  END FUNCTION is_digit
 
!**********************************************************************
 
  SUBROUTINE split (str, delims, before, sep)
 
! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the
! found delimiter. A delimiter in 'str' is treated like an ordinary
! character if it is preceded by a backslash (\). If the backslash
! character is desired in 'str', then precede it with another backslash.
 
   CHARACTER (LEN=*) :: str, delims, before
   CHARACTER, OPTIONAL :: sep
   LOGICAL :: pres
   CHARACTER :: ch, cha
 
   pres = present (sep)
   str = adjustl (str)
   CALL compact (str)
   lenstr = len_trim (str)
   IF (lenstr == 0) RETURN! string str is empty
   k = 0
   ibsl = 0 ! backslash initially inactive
   before = ' '
   DO i = 1, lenstr
    ch = str (i:i)
    IF (ibsl == 1) THEN ! backslash active
     k = k + 1
     before (k:k) = ch
     ibsl = 0
     CYCLE
    END IF
    IF (ch == '\') THEN ! backslash with backslash inactive
     k = k + 1
     before (k:k) = ch
     ibsl = 1
     CYCLE
    END IF
    ipos = index (delims, ch)
    IF (ipos == 0) THEN ! character is not a delimiter
     k = k + 1
     before (k:k) = ch
     CYCLE
    END IF
    IF (ch /= ' ') THEN ! character is a delimiter that is not a space
     str = str (i+1:)
     IF (pres) sep = ch
     EXIT
    END IF
    cha = str (i+1:i+1)! character is a space delimiter
    iposa = index (delims, cha)
    IF (iposa > 0) THEN ! next character is a delimiter
     str = str (i+2:)
     IF (pres) sep = cha
     EXIT
    ELSE
     str = str (i+1:)
     IF (pres) sep = ch
     EXIT
    END IF
   END DO
   IF (i >= lenstr) str = ''
   str = adjustl (str)! remove initial spaces
   RETURN
 
  END SUBROUTINE split
 
!**********************************************************************
 
  SUBROUTINE removebksl (str)
 
! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.
 
   CHARACTER (LEN=*) :: str
   CHARACTER (LEN=1) :: ch
   CHARACTER (LEN=len_trim(str)) :: outstr
 
   str = adjustl (str)
   lenstr = len_trim (str)
   outstr = ' '
   k = 0
   ibsl = 0 ! backslash initially inactive
 
   DO i = 1, lenstr
    ch = str (i:i)
    IF (ibsl == 1) THEN ! backslash active
     k = k + 1
     outstr (k:k) = ch
     ibsl = 0
     CYCLE
    END IF
    IF (ch == '\') THEN ! backslash with backslash inactive
     ibsl = 1
     CYCLE
    END IF
    k = k + 1
    outstr (k:k) = ch ! non-backslash with backslash inactive
   END DO
 
   str = adjustl (outstr)
 
  END SUBROUTINE removebksl
 
!**********************************************************************
 
END MODULE strings
 
 
