subroutine read_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !output
                                close_weights_t,close_count,close_count_t,error) !output

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

  open(unit=34,file=trim(sta_weight_name),form='unformatted',status='old',iostat=error)

  if(error .ne. 0) then; print *, 'Error opening station weight file ', trim(sta_weight_name), ' ', error; return; end if

  read(unit=34,iostat=error) close_meta,close_loc,close_weights,close_count, &
                             close_meta_t,close_loc_t,close_weights_t,close_count_t


  if(error .ne. 0) then; print *, 'Error reading station weight file ', trim(sta_weight_name), ' ', error; return; end if

  close(unit=34)

end subroutine read_station_weights
