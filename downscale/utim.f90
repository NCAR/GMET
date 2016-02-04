MODULE utim
 
  USE strings
  IMPLICIT NONE
 
CONTAINS
 
  SUBROUTINE parse_date (date, year, month, day, hour, Min, sec, error)
   CHARACTER (LEN=*), INTENT (IN) :: date
   INTEGER, INTENT (OUT) :: sec, Min, hour, day, month, year
   INTEGER, INTENT (OUT) :: error
 
   IF (len_trim(date) /= 14 .AND. len_trim(date) /= 8) THEN
    error = 1
    RETURN
   END IF
   CALL value (date(7:8), day, error)
   CALL value (date(5:6), month, error)
   CALL value (date(1:4), year, error)
   IF (len_trim(date) == 8) THEN
    sec = 0
    Min = 0
    hour = 0
   ELSE
    CALL value (date(13:14), sec, error)
    CALL value (date(11:12), Min, error)
    CALL value (date(9:10), hour, error)
   END IF
  END SUBROUTINE parse_date
 
 
  SUBROUTINE calendar_date (jdate, day, month, year)
 
   INTEGER, INTENT (IN) :: jdate
   INTEGER, INTENT (OUT) :: day, month, year
!   algorithm from Wikipedia: http://en.wikipedia.org/wiki/Julian_day
!   originally from Richards, E. G. (2013). Calendars. In S. E. Urban & P. K. Seidelmann, eds.
!                   Explanatory Supplement to the Astronomical Almanac, 3rd ed. (pp. 585–624).
!                   Mill Valley, Calif.: University Science Books. ISBN 978-1-89138-985-6
!                   p617-9
   INTEGER :: y = 4716, j = 1401, m = 2, n = 12, r = 4, p = 1461
   INTEGER :: v = 3, u = 5, s = 153, w = 2, B = 274277, C = - 38
   INTEGER :: f, e, g, h
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
 
  END SUBROUTINE calendar_date
 
 
  DOUBLE PRECISION FUNCTION date_to_unix (date)
 
   CHARACTER (LEN=*), INTENT (IN) :: date
   DOUBLE PRECISION :: u_day, i_day, days
   INTEGER :: sec, Min, hour, day, month, year, error
 
   CALL parse_date (date, year, month, day, hour, Min, sec, error)
 
   IF (error /= 0) THEN
    date_to_unix = - 9999.99
    RETURN
   END IF
 
   u_day = julian_date (1, 1, 1970)
   i_day = julian_date (day, month, year)
   days = i_day - u_day
 
 
 
   date_to_unix = (days*86400) + (hour*3600) + (Min*60) + sec
 
  END FUNCTION date_to_unix
 
 
  SUBROUTINE unix_to_date (itime, year, month, day, hour, Min, sec)
 
   DOUBLE PRECISION, INTENT (IN) :: itime
   INTEGER, INTENT (OUT) :: sec, Min, hour, day, month, year
   INTEGER :: u_day, i_day
 
   u_day = julian_date (1, 1, 1970)
   i_day = Int (itime/86400)
   IF (i_day < 0) THEN
    i_day = i_day - 1
   END IF
   i_day = i_day + 1
 
   CALL calendar_date (u_day+i_day, day, month, year)
 
   i_day = Mod (itime, 86400.0)
   IF (i_day < 0) THEN
    i_day = i_day + 86400
   END IF
 
   hour = i_day / 3600
   Min = (i_day/60) - (hour*60)
   sec = Mod (i_day, 60)
 
  END SUBROUTINE unix_to_date
 
 
  INTEGER FUNCTION julian_date (day, month, year)
 
   INTEGER, INTENT (IN) :: day, month, year
   DOUBLE PRECISION :: d, m, y
   INTEGER :: a, B
   DOUBLE PRECISION :: yr_corr = 0.0
 
   d = day
   m = month
   y = year
   a = 0
   B = 0
  ! there is no year 0
   IF (year < 0) THEN
    y = y + 1
    yr_corr = 0.75
   END IF
 
   IF (month <= 2) THEN
    y = y - 1.0
    m = month + 12.0
   END IF
 
   IF (y*10000.0+m*100.0+d >= 15821015.0) THEN
    a = year / 100
    B = 2 - a + (a/4)
   END IF
 
   julian_date = Int (365.25*y-yr_corr) + Int (30.6001*(m+1)) + d + 1720994.0 + B
 
  END FUNCTION julian_date
 
END MODULE utim
 
