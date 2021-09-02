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
      !if(n < k) then
      !  choose = 1
      !else if (n == k) then
      !  choose = 1
      if(n <= k) then
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

  ! routine for generating subset combinations of indexes for sampling a larger array
  subroutine comb(n,k,nk,co)

    integer(I4B), intent(in)  :: n,k,nk   ! n_total, n_train, n_trials
    integer(I4B), intent(out) :: co(:,:)  ! output combinations (ntrials, n_train)
    !    type(comb_result), dimension(:), pointer, intent(out) :: co

    ! local variables
    integer(I4B)   :: hm,s,kx
    integer*8      :: i,ix,t
    integer(I4B)   :: err

    !    hm = choose(n,k,err)
    err = 0
    hm = nk
    if(err /= 0) then
      return
    end if

    ! loop through the number of trials (hm=nk)
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
  
  function scrambled_indices(number_of_values) result(array)

    !@(#) M_random::scramble(3f): return integer array of random values 1 to N.
    integer,intent(in)    :: number_of_values
    integer,allocatable   :: array(:)
    integer               :: i, j, k, m, n
    integer               :: temp
    real                  :: u

    array=[(i,i=1,number_of_values)]

    ! The intrinsic RANDOM_NUMBER(3f) returns a real number (or an array
    ! of such) from the uniform distribution over the interval [0,1). (ie.
    ! it includes 0 but not 1.).

    ! To have a discrete uniform distribution on
    ! the integers {n, n+1, ..., m-1, m} carve the continuous distribution
    ! up into m+1-n equal sized chunks, mapping each chunk to an integer.

    ! One way is:
    !   call random_number(u)
    !   j = n + FLOOR((m+1-n)*u)  ! choose one from m-n+1 integers

    n=1
    m=number_of_values
    do k=1,2
      do i=1,m
        call random_number(u)
        j = n + FLOOR((m+1-n)*u)
        ! switch values
        temp=array(j)
        array(j)=array(i)
        array(i)=temp
      enddo
    enddo
  end function scrambled_indices

end module combination_routines
