Module utim
 
  Use strings
  Implicit None
 
Contains
 
  Subroutine parse_date (date, year, month, day, hour, Min, sec, error)
    Character (Len=*), Intent (In) :: date
    Integer, Intent (Out) :: sec, Min, hour, day, month, year
    Integer, Intent (Out) :: error
 
    If (len_trim(date) /= 14 .And. len_trim(date) /= 8) Then
      error = 1
      Return
    End If
    Call value (date(7:8), day, error)
    Call value (date(5:6), month, error)
    Call value (date(1:4), year, error)
    If (len_trim(date) == 8) Then
      sec = 0
      Min = 0
      hour = 0
    Else
      Call value (date(13:14), sec, error)
      Call value (date(11:12), Min, error)
      Call value (date(9:10), hour, error)
    End If
  End Subroutine parse_date
 
 
  Subroutine calendar_date (jdate, day, month, year)
 
    Integer, Intent (In) :: jdate
    Integer, Intent (Out) :: day, month, year
!   algorithm from Wikipedia: http://en.wikipedia.org/wiki/Julian_day
!   originally from Richards, E. G. (2013). Calendars. In S. E. Urban & P. K. Seidelmann, eds.
!                   Explanatory Supplement to the Astronomical Almanac, 3rd ed. (pp. 585–624).
!                   Mill Valley, Calif.: University Science Books. ISBN 978-1-89138-985-6
!                   p617-9
    Integer :: y = 4716, j = 1401, m = 2, n = 12, r = 4, p = 1461
    Integer :: v = 3, u = 5, s = 153, w = 2, B = 274277, C = - 38
    Integer :: f, e, g, h
    f = jdate + j + (((4*jdate+B)/146097)*3) / 4 + C
    e = r * f + v
    g = Mod (e, p) / r
    h = u * g + w
    day = Mod (h, s) / u + 1
    month = Mod (h/s+m, n) + 1
    year = e / p - y + (n+m-month) / n
 
!   integer :: a,b,c,d,e,z,alpha
!
!   z = jdate +1
!
!   if (z < 2299161) then
!      a = z
!   else
!      alpha = (z - 1867216.25) / 36524.25
!      a = z + 1 + alpha - (alpha / 4)
!   endif
!
!   b = a + 1524
!   c = (b - 122.1) / 365.25
!   d = 365.25 * c
!   e = (b - d) / 30.6001
!
!   day = b - d - INT(30.6001 * e)
!   if(e < 13.5) then
!      month = e - 1
!   else
!      month = e - 13
!   endif
!   if(month > 2.5) then
!      year = c - 4716
!   else
!      year = c - 4715
!   endif
 
  End Subroutine calendar_date
 
 
  Double Precision Function date_to_unix (date)
 
    Character (Len=*), Intent (In) :: date
    Double Precision :: u_day, i_day, days
    Integer :: sec, Min, hour, day, month, year, error
 
    Call parse_date (date, year, month, day, hour, Min, sec, error)
 
    If (error /= 0) Then
      date_to_unix = - 9999.99
      Return
    End If
 
    u_day = julian_date (1, 1, 1970)
    i_day = julian_date (day, month, year)
    days = i_day - u_day
 
 
 
    date_to_unix = (days*86400) + (hour*3600) + (Min*60) + sec
 
  End Function date_to_unix
 
 
  Subroutine unix_to_date (itime, year, month, day, hour, Min, sec)
 
    Double Precision, Intent (In) :: itime
    Integer, Intent (Out) :: sec, Min, hour, day, month, year
    Integer :: u_day, i_day
 
    u_day = julian_date (1, 1, 1970)
    i_day = Int (itime/86400)
    If (i_day < 0) Then
      i_day = i_day - 1
    End If
    i_day = i_day + 1
 
    Call calendar_date (u_day+i_day, day, month, year)
 
    i_day = Mod (itime, 86400.0)
    If (i_day < 0) Then
      i_day = i_day + 86400
    End If
 
    hour = i_day / 3600
    Min = (i_day/60) - (hour*60)
    sec = Mod (i_day, 60)
 
  End Subroutine unix_to_date
 
 
  Integer Function julian_date (day, month, year)
 
    Integer, Intent (In) :: day, month, year
    Double Precision :: d, m, y
    Integer :: a, B
    Double Precision :: yr_corr = 0.0
 
    d = day
    m = month
    y = year
    a = 0
    B = 0
  ! there is no year 0
    If (year < 0) Then
      y = y + 1
      yr_corr = 0.75
    End If
 
    If (month <= 2) Then
      y = y - 1.0
      m = month + 12.0
    End If
 
    If (y*10000.0+m*100.0+d >= 15821015.0) Then
      a = year / 100
      B = 2 - a + (a/4)
    End If
 
    julian_date = Int (365.25*y-yr_corr) + Int (30.6001*(m+1)) + d + 1720994.0 + B
 
  End Function julian_date
 
End Module utim
 
