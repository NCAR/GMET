module utim
 
  use string_mod
  implicit none
 
contains
 
  subroutine parse_date (date, year, month, day, hour, min, sec, error)
    character (len=*), intent (in) :: date
    integer, intent (out) :: sec, min, hour, day, month, year
    integer, intent (out) :: error
 
    if (len_trim(date) /= 14 .and. len_trim(date) /= 8 .and. len_trim(date) /= 10) then
      error = 1
      return
    end if
    call value (date(7:8), day, error)
    call value (date(5:6), month, error)
    call value (date(1:4), year, error)
    if (len_trim(date) == 8) then
      sec = 0
      min = 0
      hour = 0
    else if(len_trim(date) == 10) then
      sec = 0
      min = 0
      call value (date(9:10), hour, error)
    else
      call value (date(13:14), sec, error)
      call value (date(11:12), min, error)
      call value (date(9:10), hour, error)
    end if
  end subroutine parse_date
 
 
  subroutine calendar_date (jdate, day, month, year)
    !   algorithm from Wikipedia: http://en.wikipedia.org/wiki/Julian_day
    !   originally from Richards, E. G. (2013). Calendars. In S. E. Urban & P. K. Seidelmann, eds.
    !                   Explanatory Supplement to the Astronomical Almanac, 3rd ed. (pp. 585–624).
    !                   Mill Valley, Calif.: University Science Books. ISBN 978-1-89138-985-6
    !                   p617-9

    integer, intent (in) :: jdate
    integer, intent (out) :: day, month, year
    integer :: y = 4716, j = 1401, m = 2, n = 12, r = 4, p = 1461
    integer :: v = 3, u = 5, s = 153, w = 2, b = 274277, c = - 38
    integer :: f, e, g, h
    f = jdate + j + (((4*jdate+b)/146097)*3) / 4 + c
    e = r * f + v
    g = mod (e, p) / r
    h = u * g + w
    day = mod (h, s) / u + 1
    month = mod (h/s+m, n) + 1
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
 
  end subroutine calendar_date
 
 
  double precision function date_to_unix (date)
 
    character (len=*), intent (in) :: date
    double precision :: u_day, i_day, days
    integer :: sec, min, hour, day, month, year, error
 
    call parse_date (date, year, month, day, hour, min, sec, error)
 
    if (error /= 0) then
      date_to_unix = - 9999.99
      return
    end if
 
    u_day = julian_date (1, 1, 1970)
    i_day = julian_date (day, month, year)
    days = i_day - u_day
 
    date_to_unix = (days*86400) + (hour*3600) + (min*60) + sec
 
  end function date_to_unix
 
 
  subroutine unix_to_date (itime, year, month, day, hour, min, sec)
 
    double precision, intent (in) :: itime
    integer, intent (out) :: sec, min, hour, day, month, year
    integer :: u_day, i_day
 
    u_day = julian_date (1, 1, 1970)
    i_day = int (itime/86400)
    if (i_day < 0) then
      i_day = i_day - 1
    end if
    i_day = i_day + 1
 
    call calendar_date (u_day+i_day, day, month, year)

    i_day = mod (itime, 86400.0)
    if (i_day < 0) then
      i_day = i_day + 86400
    end if
 
    hour = i_day / 3600
    min = (i_day/60) - (hour*60)
    sec = mod (i_day, 60)
 
  end subroutine unix_to_date
 
 
  integer function julian_date (day, month, year)
    ! returns Julian Day (JD), where zero is noon on 1 Jan -4712 (eg 4713 BC)
 
    integer, intent (in) :: day, month, year
    double precision :: d, m, y
    integer :: a, b
    double precision :: yr_corr = 0.0
 
    d = day
    m = month
    y = year
    a = 0
    b = 0
    ! there is no year 0
    if (year < 0) then
      y = y + 1
      yr_corr = 0.75
    end if
 
    if (month <= 2) then
      y = y - 1.0
      m = month + 12.0
    end if
 
    if (y*10000.0+m*100.0+d >= 15821015.0) then
      a = year / 100
      b = 2 - a + (a/4)
    end if
 
    julian_date = int (365.25*y-yr_corr) + int (30.6001*(m+1)) + d + 1720994.0 + b
 
  end function julian_date
 
end module utim
 
