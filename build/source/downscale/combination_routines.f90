module combination_routines
  use iso_fortran_env
  use type
  implicit none 

contains

  function choose(n,k,err)
    integer*8 :: choose
    integer(I4B), intent(in) :: n,k
    integer(I4B), optional, intent(out) :: err

    !local variables
    integer*8 :: imax, i, imin, ie

    ie = 0
    if((n < 0) .or. (k < 0)) then
      write(ERROR_UNIT,*) "negative in choose"
      choose = 0
      ie = 1
    else
      if(n < k) then
        choose = 1
      else if (n == k) then
        choose = 1
      else
        imax = max(k,n-k)
        imin = min(k,n-k)
        choose = 1
        do i = imax+1,n
          choose = choose * i
        end do
        do i = 2,imin
          choose = choose/i
        end do
      end if
    end if

    if(present(err)) err = ie

  end function choose

  subroutine comb(n,k,nk,co)

    integer(I4B), intent(in)  :: n,k,nk
    integer(I4B), intent(out) :: co(:,:)
!    type(comb_result), dimension(:), pointer, intent(out) :: co

    !local variables
    integer(I4B)   :: hm,s,kx
    integer*8 :: i,ix,t
    integer(I4B) :: err

!    hm = choose(n,k,err)
    err = 0
    hm = nk
    if(err /= 0) then
      return
    end if

    do i = 1,hm
      ix = (i)*(choose(n,k)/hm)-1
!      ix = (choose(n,k))-1-i
      kx = k
      do s = 0, n-1
        if (kx == 0) exit
        t = choose(n-(s+1),kx-1)
        if(ix < t) then
          co(i,kx) = s+1
          kx = kx -1
        else
          ix = ix - t
        end if
      end do
    end do

  end subroutine comb

end module combination_routines
