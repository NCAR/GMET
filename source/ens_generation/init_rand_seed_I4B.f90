! ======= generate one random integer seed =======
! Author:  A. Wood, 2016
! Adapted from base fortran resources found on the web

subroutine init_rand_seed_I4B(one_seed)
            use iso_fortran_env, only: int32
            implicit none

            integer(int32), intent(out) :: one_seed
            integer :: i, n, un, istat, dt(8), pid
            integer(int32) :: t
            integer :: getpid  ! AW added

            ! AWW removed
            !            call random_seed(size = n)
            !            allocate(seed(n))

            ! First try if the OS provides a random number generator

            ! AWW -- the file check is too slow ... 
            !   just go straight to time clock/pid approach

!            open(newunit=un, file="/dev/urandom", access="stream", &
!                 form="unformatted", action="read", status="old", iostat=istat)
!            if (istat == 0) then
!               read(un) seed
!               close(un)
!            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  ! changed to int32
                  t = (dt(1) - 1970) * 365_int32 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int32 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int32 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               one_seed = lcg(t)

!            end if
!            call random_seed(put=one_seed)

          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            !  changed to int32 in settings below
            function lcg(s)
              integer :: lcg
              integer(int32) :: s

              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 2147483647_int32)

              end if
              s = mod(s * 139735137_int32, 2147483647_int32)
              lcg = int(mod(s, int(huge(0), int32)), kind(0))
            end function lcg

end subroutine init_rand_seed_I4B
