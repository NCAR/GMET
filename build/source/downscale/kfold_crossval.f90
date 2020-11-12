subroutine kfold_crossval(X,Y,W,kfold_trials,kfold_nsamp,kfold_hold,xval_combinations,varUncert)
  use type
  implicit none

! ===== start interface =======
  interface
    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares
  end interface


  !inputs/outputs 
  real(dp), intent(in)      :: X(:,:)                    ! full station attribute array
  real(dp), intent(in)      :: Y(:)                    ! full station variable value vector
  real(dp), intent(in)      :: W(:,:)                  ! full station diagonal weight matrix
  integer(I4B), intent(in)  :: kfold_trials            ! number of kfold xval trials
  integer(I4B), intent(in)  :: kfold_nsamp             ! number of samples to draw from
  integer(I4B), intent(in)  :: kfold_hold              ! number of stations to withhold from regression
  integer(I4B), intent(in)  :: xval_combinations(:,:)  ! array of sampling combinations (integer array indicies)
  real(sp), intent(out)     :: varUncert               ! output: uncertainty estimate from kfold trials
  
  !local variables
  real(dp), allocatable     :: w_1d(:)               ! 1d vector of the weights from diagonal weight matrix W
  real(dp)                  :: weightSum             ! summation of all station weights in error estimate
  real(dp)                  :: errorSum              ! summation of errors from kfold trials
  real(dp)                  :: varTmp                ! temporary variable value
  real(dp)                  :: meanVal               ! mean value of withheld obs
  real(dp), allocatable     :: x_xval(:,:)           ! station predictor array for kfold xval
  real(dp), allocatable     :: y_xval(:)             ! station variable values for kfold xval
  real(dp), allocatable     :: w_xval(:,:)           ! diagonal weight matrix for trial stations
  real(dp), allocatable     :: tx_xval(:,:)          ! transposed station predictor array for kfold xval
  real(dp), allocatable     :: twx_xval(:,:)         ! transposed matmul(tx,w) array for kfold xval
  real(dp), allocatable     :: vv(:)                 ! temporary matrix to test for singularities
  real(dp), allocatable     :: B(:)                  ! regression coefficients for predictors

  integer(I4B)              :: i, j                  ! counter variables
  integer(I4B)              :: nPredict              ! number of predictors in attribute arrays
  integer(I4B), allocatable :: xval_sta_list(:)      ! vector of one sampling combination (integer array indicies)
  integer(I4B), allocatable :: xval_sta_list_fill(:) ! vector of sampling combination zero padded to length kfold_nsamp
  integer(I4B), allocatable :: withheld_sta_list(:)  ! indices of withheld stations
  integer(I4B), allocatable :: tmp_inds(:)           ! temporary vector of indices
  integer(I4B), allocatable :: all_inds(:)           ! vector of possible station indices (1-sta_limit)


  !initialize variables
  weightSum = 0.0
  errorSum  = 0.0

  nPredict = size(X,2)

  !allocations
  allocate(xval_sta_list(kfold_nsamp-kfold_hold))
  allocate(xval_sta_list_fill(kfold_nsamp))
  allocate(withheld_sta_list(kfold_hold))
  allocate(tmp_inds(kfold_nsamp))
  allocate(all_inds(kfold_nsamp))
  allocate(x_xval(kfold_nsamp-kfold_hold,nPredict))
  allocate(y_xval(kfold_nsamp-kfold_hold))
  allocate(w_xval(kfold_nsamp-kfold_hold,kfold_nsamp-kfold_hold))
  allocate(w_1d(kfold_nsamp))
  allocate(tx_xval(nPredict,kfold_nsamp-kfold_hold))
  allocate(twx_xval(nPredict,kfold_nsamp-kfold_hold))
  allocate(B(nPredict))
  allocate(vv(nPredict))

  !initalize other variables
  do i = 1,kfold_nsamp,1
    all_inds(i) = i
    w_1d(i) = W(i,i)
  end do

  !loop through the kfold trials
  do i = 1,kfold_trials
!  print *,'kfold:', i
    !current kfold trial sampling combination
    xval_sta_list = xval_combinations(i,:)
    !create a zero padded vector of above
    xval_sta_list_fill(1:kfold_nsamp-kfold_hold) = xval_sta_list

    !find withheld stations
    tmp_inds = pack(all_inds,all_inds .ne. xval_sta_list_fill)
    withheld_sta_list = tmp_inds(1:kfold_hold)


    !station predictor array
    x_xval = X(xval_sta_list,:)
    !station values
    y_xval = Y(xval_sta_list)
    !station diagonal weight matrix
    w_xval = W(xval_sta_list,xval_sta_list)
    !create transposed matrices for least squares
    tx_xval  = transpose(x_xval)
    twx_xval = matmul(tx_xval,w_xval)
    
    vv = maxval(abs(twx_xval),dim=2)

    if(any(vv==0.0)) then 
       meanVal = sum(Y(withheld_sta_list))/kfold_hold
       do j = 1,kfold_hold
         errorSum = errorSum + w_1d(withheld_sta_list(j))*(meanVal - Y(withheld_sta_list(j)))**2
       end do
       weightSum = weightSum + sum(W_1d(withheld_sta_list))
    else
      do j = 1,kfold_hold
        call least_squares(x_xval,y_xval,twx_xval,B)
        !regression precip
        varTmp = real( dot_product( X(withheld_sta_list(j),:),B ), kind(sp) )
        !error total
        errorSum = errorSum + w_1d(withheld_sta_list(j))*(varTmp - Y(withheld_sta_list(j)))**2
      end do

      !create weight sum of in regression station weights
      weightSum = weightSum + sum(W_1d(withheld_sta_list))
    end if
  end do
  !mean uncertainty estimate from all kfold trials
  varUncert = real((errorSum/weightSum)**(1.0/2.0),kind(sp))

end subroutine kfold_crossval