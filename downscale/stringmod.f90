Module strings
 
  Use precision
 
  Private :: value_dr, value_sr, value_di, value_si
  Private :: write_dr, write_sr, write_di, write_si
  Private :: writeq_dr, writeq_sr, writeq_di, writeq_si
 
  Interface value ! Generic operator for converting a number string to a
                 ! number. Calling syntax is 'call value(numstring,number,ios)'
                 ! where 'numstring' is a number string and 'number' is a
                 ! real number or an integer (single or double precision).
    Module Procedure value_dr
    Module Procedure value_sr
    Module Procedure value_di
    Module Procedure value_si
  End Interface
 
  Interface writenum ! Generic  interface for writing a number to a string. The
                    ! number is left justified in the string. The calling syntax
                    ! is 'call writenum(number,string,format)' where 'number' is
                    ! a real number or an integer, 'string' is a character string
                    ! containing the result, and 'format' is the format desired,
                    ! e.g., 'e15.6' or 'i5'.
    Module Procedure write_dr
    Module Procedure write_sr
    Module Procedure write_di
    Module Procedure write_si
  End Interface
 
  Interface writeq ! Generic interface equating a name to a numerical value. The
                  ! calling syntax is 'call writeq(unit,name,value,format)' where
                  ! unit is the integer output unit number, 'name' is the variable
                  ! name, 'value' is the real or integer value of the variable,
                  ! and 'format' is the format of the value. The result written to
                  ! the output unit has the form <name> = <value>.
    Module Procedure writeq_dr
    Module Procedure writeq_sr
    Module Procedure writeq_di
    Module Procedure writeq_si
  End Interface
 
 
!**********************************************************************
 
Contains
 
!**********************************************************************
 
  Subroutine parse (str, delims, args, nargs)
 
! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.
 
    Character (Len=*) :: str, delims
    Character (Len=len_trim(str)) :: strsav
    Character (Len=*), Dimension (:) :: args
 
    strsav = str
    Call compact (str)
    na = size (args)
    Do i = 1, na
      args (i) = ' '
    End Do
    nargs = 0
    lenstr = len_trim (str)
    If (lenstr == 0) Return
    k = 0
 
    Do
      If (len_trim(str) == 0) Exit
      nargs = nargs + 1
      Call split (str, delims, args(nargs))
      Call removebksl (args(nargs))
    End Do
    str = strsav
 
  End Subroutine parse
 
!**********************************************************************
 
  Subroutine compact (str)
 
! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.
 
    Character (Len=*) :: str
    Character (Len=1) :: ch
    Character (Len=len_trim(str)) :: outstr
 
    str = adjustl (str)
    lenstr = len_trim (str)
    outstr = ' '
    isp = 0
    k = 0
 
    Do i = 1, lenstr
      ch = str (i:i)
      ich = iachar (ch)
 
      Select Case (ich)
 
      Case (9, 32)! space or tab character
        If (isp == 0) Then
          k = k + 1
          outstr (k:k) = ' '
        End If
        isp = 1
 
      Case (33:)! not a space, quote, or control character
        k = k + 1
        outstr (k:k) = ch
        isp = 0
 
      End Select
 
    End Do
 
    str = adjustl (outstr)
 
  End Subroutine compact
 
!**********************************************************************
 
  Subroutine removesp (str)
 
! Removes spaces, tabs, and control characters in string str
 
    Character (Len=*) :: str
    Character (Len=1) :: ch
    Character (Len=len_trim(str)) :: outstr
 
    str = adjustl (str)
    lenstr = len_trim (str)
    outstr = ' '
    k = 0
 
    Do i = 1, lenstr
      ch = str (i:i)
      ich = iachar (ch)
      Select Case (ich)
      Case (0:32)! space, tab, or control character
        Cycle
      Case (33:)
        k = k + 1
        outstr (k:k) = ch
      End Select
    End Do
 
    str = adjustl (outstr)
 
  End Subroutine removesp
 
!**********************************************************************
 
  Subroutine value_dr (str, rnum, ios)
 
! Converts number string to a double precision real number
 
    Character (Len=*) :: str
    Real (kr8) :: rnum
    Integer :: ios
 
    ilen = len_trim (str)
    ipos = scan (str, 'Ee')
    If ( .Not. is_digit(str(ilen:ilen)) .And. ipos /= 0) Then
      ios = 3
      Return
    End If
    Read (str,*, IoStat=ios) rnum
 
  End Subroutine value_dr
 
!**********************************************************************
 
  Subroutine value_sr (str, rnum, ios)
 
! Converts number string to a single precision real number
 
    Character (Len=*) :: str
    Real (kr4) :: rnum
    Real (kr8) :: rnumd
 
    Call value_dr (str, rnumd, ios)
    If (Abs(rnumd) > huge(rnum)) Then
      ios = 15
      Return
    End If
    If (Abs(rnumd) < tiny(rnum)) rnum = 0.0_kr4
    rnum = rnumd
 
  End Subroutine value_sr
 
!**********************************************************************
 
  Subroutine value_di (str, inum, ios)
 
! Converts number string to a double precision integer value
 
    Character (Len=*) :: str
    Integer (ki8) :: inum
    Real (kr8) :: rnum
 
    Call value_dr (str, rnum, ios)
    If (Abs(rnum) > huge(inum)) Then
      ios = 15
      Return
    End If
    inum = Nint (rnum, ki8)
 
  End Subroutine value_di
 
!**********************************************************************
 
  Subroutine value_si (str, inum, ios)
 
! Converts number string to a single precision integer value
 
    Character (Len=*) :: str
    Integer (ki4) :: inum
    Real (kr8) :: rnum
 
    Call value_dr (str, rnum, ios)
    If (Abs(rnum) > huge(inum)) Then
      ios = 15
      Return
    End If
    inum = Nint (rnum, ki4)
 
  End Subroutine value_si
 
!**********************************************************************
 
  Subroutine shiftstr (str, n)
 
! Shifts characters in in the string 'str' n positions (positive values
! denote a right shift and negative values denote a left shift). Characters
! that are shifted off the end are lost. Positions opened up by the shift
! are replaced by spaces.
 
    Character (Len=*) :: str
 
    lenstr = len (str)
    nabs = iabs (n)
    If (nabs >= lenstr) Then
      str = repeat (' ', lenstr)
      Return
    End If
    If (n < 0) str = str (nabs+1:) // repeat (' ', nabs)! shift left
    If (n > 0) str = repeat (' ', nabs) // str (:lenstr-nabs)! shift right
    Return
 
  End Subroutine shiftstr
 
!**********************************************************************
 
  Subroutine insertstr (str, strins, loc)
 
! Inserts the string 'strins' into the string 'str' at position 'loc'.
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. Trailing spaces of 'strins' are
! removed prior to insertion
 
    Character (Len=*) :: str, strins
    Character (Len=Len(str)) :: tempstr
 
    lenstrins = len_trim (strins)
    tempstr = str (loc:)
    Call shiftstr (tempstr, lenstrins)
    tempstr (1:lenstrins) = strins (1:lenstrins)
    str (loc:) = tempstr
    Return
 
  End Subroutine insertstr
 
!**********************************************************************
 
  Subroutine delsubstr (str, substr)
 
! Deletes first occurrence of substring 'substr' from string 'str' and
! shifts characters left to fill hole. Trailing spaces or blanks are
! not considered part of 'substr'.
 
    Character (Len=*) :: str, substr
 
    lensubstr = len_trim (substr)
    ipos = index (str, substr)
    If (ipos == 0) Return
    If (ipos == 1) Then
      str = str (lensubstr+1:)
    Else
      str = str (:ipos-1) // str (ipos+lensubstr:)
    End If
    Return
 
  End Subroutine delsubstr
 
!**********************************************************************
 
  Subroutine delall (str, substr)
 
! Deletes all occurrences of substring 'substr' from string 'str' and
! shifts characters left to fill holes.
 
    Character (Len=*) :: str, substr
 
    lensubstr = len_trim (substr)
    Do
      ipos = index (str, substr)
      If (ipos == 0) Exit
      If (ipos == 1) Then
        str = str (lensubstr+1:)
      Else
        str = str (:ipos-1) // str (ipos+lensubstr:)
      End If
    End Do
    Return
 
  End Subroutine delall
 
!**********************************************************************
 
  Function uppercase (str) Result (ucstr)
 
! convert string to upper case
 
    Character (Len=*) :: str
    Character (Len=len_trim(str)) :: ucstr
 
    ilen = len_trim (str)
    ioffset = iachar ('A') - iachar ('a')
    iquote = 0
    ucstr = str
    Do i = 1, ilen
      iav = iachar (str(i:i))
      If (iquote == 0 .And. (iav == 34 .Or. iav == 39)) Then
        iquote = 1
        iqc = iav
        Cycle
      End If
      If (iquote == 1 .And. iav == iqc) Then
        iquote = 0
        Cycle
      End If
      If (iquote == 1) Cycle
      If (iav >= iachar('a') .And. iav <= iachar('z')) Then
        ucstr (i:i) = achar (iav+ioffset)
      Else
        ucstr (i:i) = str (i:i)
      End If
    End Do
    Return
 
  End Function uppercase
 
!**********************************************************************
 
  Function lowercase (str) Result (lcstr)
 
! convert string to lower case
 
    Character (Len=*) :: str
    Character (Len=len_trim(str)) :: lcstr
 
    ilen = len_trim (str)
    ioffset = iachar ('A') - iachar ('a')
    iquote = 0
    lcstr = str
    Do i = 1, ilen
      iav = iachar (str(i:i))
      If (iquote == 0 .And. (iav == 34 .Or. iav == 39)) Then
        iquote = 1
        iqc = iav
        Cycle
      End If
      If (iquote == 1 .And. iav == iqc) Then
        iquote = 0
        Cycle
      End If
      If (iquote == 1) Cycle
      If (iav >= iachar('A') .And. iav <= iachar('Z')) Then
        lcstr (i:i) = achar (iav-ioffset)
      Else
        lcstr (i:i) = str (i:i)
      End If
    End Do
    Return
 
  End Function lowercase
 
!**********************************************************************
 
  Subroutine readline (nunitr, line, ios)
 
! Reads line from unit=nunitr, ignoring blank lines
! and deleting comments beginning with an exclamation point(!)
 
    Character (Len=*) :: line
 
    Do
      Read (nunitr, '(a)', IoStat=ios) line ! read input line
      If (ios /= 0) Return
      line = adjustl (line)
      ipos = index (line, '!')
      If (ipos == 1) Cycle
      If (ipos /= 0) line = line (:ipos-1)
      If (len_trim(line) /= 0) Exit
    End Do
    Return
 
  End Subroutine readline
 
!**********************************************************************
 
  Subroutine match (str, ipos, imatch)
 
! Sets imatch to the position in string of the delimiter matching the delimiter
! in position ipos. Allowable delimiters are (), [], {}, <>.
 
    Character (Len=*) :: str
    Character :: delim1, delim2, ch
 
    lenstr = len_trim (str)
    delim1 = str (ipos:ipos)
    Select Case (delim1)
    Case ('(')
      idelim2 = iachar (delim1) + 1
      istart = ipos + 1
      iend = lenstr
      inc = 1
    Case (')')
      idelim2 = iachar (delim1) - 1
      istart = ipos - 1
      iend = 1
      inc = - 1
    Case ('[', '{', '<')
      idelim2 = iachar (delim1) + 2
      istart = ipos + 1
      iend = lenstr
      inc = 1
    Case (']', '}', '>')
      idelim2 = iachar (delim1) - 2
      istart = ipos - 1
      iend = 1
      inc = - 1
    Case Default
      Write (*,*) delim1, ' is not a valid delimiter'
      Return
    End Select
    If (istart < 1 .Or. istart > lenstr) Then
      Write (*,*) delim1, ' has no matching delimiter'
      Return
    End If
    delim2 = achar (idelim2)! matching delimiter
 
    isum = 1
    Do i = istart, iend, inc
      ch = str (i:i)
      If (ch /= delim1 .And. ch /= delim2) Cycle
      If (ch == delim1) isum = isum + 1
      If (ch == delim2) isum = isum - 1
      If (isum == 0) Exit
    End Do
    If (isum /= 0) Then
      Write (*,*) delim1, ' has no matching delimiter'
      Return
    End If
    imatch = i
 
    Return
 
  End Subroutine match
 
!**********************************************************************
 
  Subroutine write_dr (rnum, str, fmt)
 
! Writes double precision real number rnum to string str using format fmt
 
    Real (kr8) :: rnum
    Character (Len=*) :: str, fmt
    Character (Len=80) :: formt
 
    formt = '(' // trim (fmt) // ')'
    Write (str, formt) rnum
    str = adjustl (str)
 
  End Subroutine write_dr
 
!***********************************************************************
 
  Subroutine write_sr (rnum, str, fmt)
 
! Writes single precision real number rnum to string str using format fmt
 
    Real (kr4) :: rnum
    Character (Len=*) :: str, fmt
    Character (Len=80) :: formt
 
    formt = '(' // trim (fmt) // ')'
    Write (str, formt) rnum
    str = adjustl (str)
 
  End Subroutine write_sr
 
!***********************************************************************
 
  Subroutine write_di (inum, str, fmt)
 
! Writes double precision integer inum to string str using format fmt
 
    Integer (ki8) :: inum
    Character (Len=*) :: str, fmt
    Character (Len=80) :: formt
 
    formt = '(' // trim (fmt) // ')'
    Write (str, formt) inum
    str = adjustl (str)
 
  End Subroutine write_di
 
!***********************************************************************
 
  Subroutine write_si (inum, str, fmt)
 
! Writes single precision integer inum to string str using format fmt
 
    Integer (ki4) :: inum
    Character (Len=*) :: str, fmt
    Character (Len=80) :: formt
 
    formt = '(' // trim (fmt) // ')'
    Write (str, formt) inum
    str = adjustl (str)
 
  End Subroutine write_si
 
!***********************************************************************
 
  Subroutine trimzero (str)
 
! Deletes nonsignificant trailing zeroes from number string str. If number
! string ends in a decimal point, one trailing zero is added.
 
    Character (Len=*) :: str
    Character :: ch
    Character (Len=10) :: Exp
 
    ipos = scan (str, 'eE')
    If (ipos > 0) Then
      Exp = str (ipos:)
      str = str (1:ipos-1)
    End If
    lstr = len_trim (str)
    Do i = lstr, 1, - 1
      ch = str (i:i)
      If (ch == '0') Cycle
      If (ch == '.') Then
        str = str (1:i) // '0'
        If (ipos > 0) str = trim (str) // trim (Exp)
        Exit
      End If
      str = str (1:i)
      Exit
    End Do
    If (ipos > 0) str = trim (str) // trim (Exp)
 
  End Subroutine trimzero
 
!**********************************************************************
 
  Subroutine writeq_dr (unit, namestr, value, fmt)
 
! Writes a string of the form <name> = value to unit
 
    Real (kr8) :: value
    Integer :: unit
    Character (Len=*) :: namestr, fmt
    Character (Len=32) :: tempstr
 
    Call writenum (value, tempstr, fmt)
    Call trimzero (tempstr)
    Write (Unit,*) trim (namestr) // ' = ' // trim (tempstr)
 
  End Subroutine writeq_dr
 
!**********************************************************************
 
  Subroutine writeq_sr (unit, namestr, value, fmt)
 
! Writes a string of the form <name> = value to unit
 
    Real (kr4) :: value
    Integer :: unit
    Character (Len=*) :: namestr, fmt
    Character (Len=32) :: tempstr
 
    Call writenum (value, tempstr, fmt)
    Call trimzero (tempstr)
    Write (Unit,*) trim (namestr) // ' = ' // trim (tempstr)
 
  End Subroutine writeq_sr
 
!**********************************************************************
 
  Subroutine writeq_di (unit, namestr, ivalue, fmt)
 
! Writes a string of the form <name> = ivalue to unit
 
    Integer (ki8) :: ivalue
    Integer :: unit
    Character (Len=*) :: namestr, fmt
    Character (Len=32) :: tempstr
    Call writenum (ivalue, tempstr, fmt)
    Call trimzero (tempstr)
    Write (Unit,*) trim (namestr) // ' = ' // trim (tempstr)
 
  End Subroutine writeq_di
 
!**********************************************************************
 
  Subroutine writeq_si (unit, namestr, ivalue, fmt)
 
! Writes a string of the form <name> = ivalue to unit
 
    Integer (ki4) :: ivalue
    Integer :: unit
    Character (Len=*) :: namestr, fmt
    Character (Len=32) :: tempstr
    Call writenum (ivalue, tempstr, fmt)
    Call trimzero (tempstr)
    Write (Unit,*) trim (namestr) // ' = ' // trim (tempstr)
 
  End Subroutine writeq_si
 
!**********************************************************************
 
  Function is_letter (ch) Result (res)
 
! Returns .true. if ch is a letter and .false. otherwise
 
    Character :: ch
    Logical :: res
 
    Select Case (ch)
    Case ('A':'Z', 'a':'z')
      res = .True.
    Case Default
      res = .False.
    End Select
    Return
 
  End Function is_letter
 
!**********************************************************************
 
  Function is_digit (ch) Result (res)
 
! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise
 
    Character :: ch
    Logical :: res
 
    Select Case (ch)
    Case ('0':'9')
      res = .True.
    Case Default
      res = .False.
    End Select
    Return
 
  End Function is_digit
 
!**********************************************************************
 
  Subroutine split (str, delims, before, sep)
 
! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the
! found delimiter. A delimiter in 'str' is treated like an ordinary
! character if it is preceded by a backslash (\). If the backslash
! character is desired in 'str', then precede it with another backslash.
 
    Character (Len=*) :: str, delims, before
    Character, Optional :: sep
    Logical :: pres
    Character :: ch, cha
 
    pres = present (sep)
    str = adjustl (str)
    Call compact (str)
    lenstr = len_trim (str)
    If (lenstr == 0) Return! string str is empty
    k = 0
    ibsl = 0 ! backslash initially inactive
    before = ' '
    Do i = 1, lenstr
      ch = str (i:i)
      If (ibsl == 1) Then ! backslash active
        k = k + 1
        before (k:k) = ch
        ibsl = 0
        Cycle
      End If
      If (ch == '\') Then ! backslash with backslash inactive
        k = k + 1
        before (k:k) = ch
        ibsl = 1
        Cycle
      End If
      ipos = index (delims, ch)
      If (ipos == 0) Then ! character is not a delimiter
        k = k + 1
        before (k:k) = ch
        Cycle
      End If
      If (ch /= ' ') Then ! character is a delimiter that is not a space
        str = str (i+1:)
        If (pres) sep = ch
        Exit
      End If
      cha = str (i+1:i+1)! character is a space delimiter
      iposa = index (delims, cha)
      If (iposa > 0) Then ! next character is a delimiter
        str = str (i+2:)
        If (pres) sep = cha
        Exit
      Else
        str = str (i+1:)
        If (pres) sep = ch
        Exit
      End If
    End Do
    If (i >= lenstr) str = ''
    str = adjustl (str)! remove initial spaces
    Return
 
  End Subroutine split
 
!**********************************************************************
 
  Subroutine removebksl (str)
 
! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.
 
    Character (Len=*) :: str
    Character (Len=1) :: ch
    Character (Len=len_trim(str)) :: outstr
 
    str = adjustl (str)
    lenstr = len_trim (str)
    outstr = ' '
    k = 0
    ibsl = 0 ! backslash initially inactive
 
    Do i = 1, lenstr
      ch = str (i:i)
      If (ibsl == 1) Then ! backslash active
        k = k + 1
        outstr (k:k) = ch
        ibsl = 0
        Cycle
      End If
      If (ch == '\') Then ! backslash with backslash inactive
        ibsl = 1
        Cycle
      End If
      k = k + 1
      outstr (k:k) = ch ! non-backslash with backslash inactive
    End Do
 
    str = adjustl (outstr)
 
  End Subroutine removebksl
 
!**********************************************************************
 
End Module strings
 
 
