subroutine compute_station_weights(sta_weight_name,ngrid,nstns,X,Z,search_distance, &
                                   sta_limit,sta_data,tair_data, &
                                   close_meta,close_meta_t,close_loc,close_loc_t, &
                                   close_count,close_count_t,close_weights,close_weights_t,error)

  !subroutine to compute the grid - station weight matrix for every grid point for all stations considered for each grid point

  use type

  implicit none

  !declarations
  !inputs
  character(len=500),intent(in) :: sta_weight_name     !name of station weight binary file
  integer(I4B), intent(in)      :: ngrid               !number of grid points
  integer(I4B), intent(in)      :: nstns               !number of stations
  real(DP), intent(in)          :: Z(:,:)              !grid metadata array
  real(DP), intent(in)          :: X(:,:)              !station metadata array
  real(DP), intent(in)          :: search_distance     !default station search distance
  integer(I4B), intent(in)      :: sta_limit           !maximum number of stations for a grid point
  real(DP), intent(in)          :: sta_data(:,:)       !station data values for precipitation
  real(DP), intent(in)          :: tair_data(:,:,:)    !station air temperature data

  !in/out
  real(DP), intent(inout)     :: close_meta(:,:,:)
  real(DP), intent(inout)     :: close_meta_t(:,:,:)
  integer(I4B), intent(inout) :: close_loc(:,:)
  integer(I4B), intent(inout) :: close_loc_t(:,:)
  integer(I4B), intent(inout) :: close_count(:)
  integer(I4B), intent(inout) :: close_count_t(:)
  real(DP), intent(inout)     :: close_weights(:,:)
  real(DP), intent(inout)     :: close_weights_t(:,:)
  integer(I4B), intent(inout) :: error
  
  !local variables
  real(DP)                    :: min_weight_t
  real(DP)                    :: min_weight
  real(DP), allocatable       :: w_base(:,:)
  integer(I4B)                :: out_loc
  integer(I4B)                :: out_loc_t
  integer(I4B)                :: i,g        ! counter variables

  !code starts below
  
  !allocate base weight array
  allocate(w_base(ngrid,nstns),stat=error)
  if(error/=0)then; print *,'Error allocating w_base',error; return; endif

  !set to zero
  w_base = 0.0d0

  do g = 1, ngrid, 1
    if(z(g,4) > -400.) then
      close_count(g) = 1
      min_weight = 0.0d0
      close_weights(g,:) = 0.0d0

      close_count_t(g) = 1
      min_weight_t = 0.0d0
      close_weights_t(g,:) = 0.0d0

      ! for current grid cell, loop through stations, find distance
      do i = 1, nstns, 1
        !setup distinct weight matrices for precip and temperature
        ! x() are station lonlat; z() are grid lonlat; returns weight (w_base) for grd-to-stn
        call calc_distance_weight(search_distance, X(i,2), X(i,3), Z(g,2), Z(g,3), w_base(g,i))

        !Precipitation
        min_weight = 0.0d0

        if(w_base(g,i) .gt. min_weight .and. sta_data(i,1) .gt. -1.0d0) then
          if(close_count(g) .le. sta_limit) then

            close_weights(g,close_count(g)) = w_base(g,i)
            close_loc(g,close_count(g))  = i

            close_meta(1,g,close_count(g)) = X(i,2)
            close_meta(2,g,close_count(g)) = X(i,3)
            close_meta(3,g,close_count(g)) = Z(g,2)
            close_meta(4,g,close_count(g)) = Z(g,3)
            call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta(5,g,close_count(g)))

            close_count(g) = close_count(g) + 1
          else
            min_weight = minval(close_weights(g,:),1)
            if(w_base(g,i) .gt. min_weight) then
              out_loc = minloc(close_weights(g,:),1)
              close_weights(g,out_loc) = w_base(g,i)
              close_loc(g,out_loc) = i
              close_meta(1,g,out_loc) = X(i,2)
              close_meta(2,g,out_loc) = X(i,3)
              close_meta(3,g,out_loc) = Z(g,2)
              close_meta(4,g,out_loc) = Z(g,3)
              call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta(5,g,out_loc))
            endif
          endif
        endif

        !need to repeat above for temperature since that data is independent of precipitation
        min_weight_t = 0.0d0

        if(w_base(g,i) .gt. min_weight_t .and. tair_data(1,i,1) .gt. -200.0d0) then
          if(close_count_t(g) .le. sta_limit) then
  
            close_weights_t(g,close_count_t(g)) = w_base(g,i)
            close_loc_t(g,close_count_t(g))  = i

            close_meta_t(1,g,close_count_t(g)) = X(i,2)
            close_meta_t(2,g,close_count_t(g)) = X(i,3)
            close_meta_t(3,g,close_count_t(g)) = Z(g,2)
            close_meta_t(4,g,close_count_t(g)) = Z(g,3)
            call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta_t(5,g,close_count_t(g)))

            close_count_t(g) = close_count_t(g) + 1
          else
            min_weight_t = minval(close_weights_t(g,:),1)
            if(w_base(g,i) .gt. min_weight_t) then
              out_loc_t = minloc(close_weights_t(g,:),1)
              close_weights_t(g,out_loc_t) = w_base(g,i)

              close_loc_t(g,out_loc_t) = i

              close_meta_t(1,g,out_loc_t) = X(i,2)
              close_meta_t(2,g,out_loc_t) = X(i,3)
              close_meta_t(3,g,out_loc_t) = Z(g,2)
              close_meta_t(4,g,out_loc_t) = Z(g,3)
              call calc_distance(X(i,2),X(i,3),Z(g,2),Z(g,3),close_meta_t(5,g,out_loc_t))
            endif
          endif
        endif
      enddo  !end station loop
    endif    !end grid point elevation check
  enddo      !end grid point loop

end subroutine compute_station_weights


subroutine write_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !input
                                close_weights_t,close_count,close_count_t,error) !input

  !subroutine to write the grid-station weight matrix to a fortran binary file

  use type

  implicit none

  !declarations
  !inputs
  character(len=500),intent(in) :: sta_weight_name     !name of station weight binary file
  real(DP), intent(in)     :: close_meta(:,:,:)
  real(DP), intent(in)     :: close_meta_t(:,:,:)
  integer(I4B), intent(in) :: close_loc(:,:)
  integer(I4B), intent(in) :: close_loc_t(:,:)
  integer(I4B), intent(in) :: close_count(:)
  integer(I4B), intent(in) :: close_count_t(:)
  real(DP), intent(in)     :: close_weights(:,:)
  real(DP), intent(in)     :: close_weights_t(:,:)
  !in/out
  integer(I4B), intent(inout) :: error

  !output weight variables to file
  open(unit=34,file=trim(sta_weight_name),form='unformatted',iostat=error)

  if(error .ne. 0) then; print *, 'Error opening station weight file ', trim(sta_weight_name), ' ', error; stop; end if

  write(unit=34,iostat=error) close_meta,close_loc,close_weights,close_count, &
                              close_meta_t,close_loc_t,close_weights_t,close_count_t

  if(error .ne. 0) then; print *, 'Error writing station weight file ', trim(sta_weight_name), ' ', error; stop; end if

  close(unit=34)

end subroutine write_station_weights

subroutine read_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !output
                                close_weights_t,close_count,close_count_t,error) !output

  !subroutine to read the grid-station weight matrix from fortran binary file

  use type
  implicit none

  !input
  character(len=500), intent(in)    :: sta_weight_name

  !output
  real(DP), intent(out)      :: close_meta(:,:,:)
  real(DP), intent(out)      :: close_meta_t(:,:,:)
  integer(I4B), intent(out)  :: close_loc(:,:)
  integer(I4B), intent(out)  :: close_loc_t(:,:)
  real(DP), intent(out)      :: close_weights(:,:)
  real(DP), intent(out)      :: close_weights_t(:,:)
  integer(I4B), intent(out)  :: close_count(:)
  integer(I4B), intent(out)  :: close_count_t(:)

  !in/out
  integer(I4B), intent(inout):: error

  !code starts below

  open(unit=34,file=trim(sta_weight_name),form='unformatted',iostat=error)

  if(error .ne. 0) then; print *, 'Error opening station weight file ', trim(sta_weight_name), ' ', error; stop; end if

  read(unit=34,iostat=error) close_meta,close_loc,close_weights,close_count, &
                             close_meta_t,close_loc_t,close_weights_t,close_count_t


  if(error .ne. 0) then; print *, 'Error reading station weight file ', trim(sta_weight_name), ' ', error; stop; end if

  close(unit=34)

end subroutine read_station_weights
