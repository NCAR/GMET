! performs least squares regression in an iterative process over k folds of the data, returns predictive error

subroutine kfold_crossval(X, Y, W, kfold_trials, n_total, n_train, max_n_test, xval_indices, varUncert)
  use type
  implicit none

  ! ===== start interface =======
  interface
    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in)               :: x (:, :)
      real (dp), intent (in)               :: y (:)
      real (dp), intent (in)               :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares
  end interface

  ! inputs/outputs 
  real(dp), intent(in)      :: X(:,:)                ! full station attribute array
  real(dp), intent(in)      :: Y(:)                  ! full station variable value vector
  real(dp), intent(in)      :: W(:,:)                ! full station diagonal weight matrix
  integer(I4B), intent(in)  :: kfold_trials          ! number of kfold xval trials
  integer(I4B), intent(in)  :: n_total               ! size of station sample to draw from (sta_limit)
  integer(I4B), intent(in)  :: n_train               ! number in training sample
  integer(i4b), intent(in)  :: max_n_test            ! maximum test sample size over all trials
  integer(I4B), intent(in)  :: xval_indices(:,:)     ! array of sampling combinations (integer array indices)
  real(sp), intent(out)     :: varUncert             ! output: uncertainty estimate from kfold trials
  
  ! local variables
  real(dp), allocatable     :: w_1d(:)               ! 1d vector of the station weights from diagonal weight matrix W
  real(dp)                  :: weightSum             ! summation of all station weights in error estimate
  real(dp)                  :: errorSum              ! summation of errors from kfold trials
  real(dp)                  :: varTmp                ! temporary variable value
  real(dp)                  :: meanVal               ! mean value of withheld obs
  real(dp), allocatable     :: x_xval_train(:,:)     ! station predictor array for kfold xval
  real(dp), allocatable     :: y_xval_train(:)       ! station variable values for kfold xval
  real(dp), allocatable     :: w_xval_train(:,:)     ! diagonal weight matrix for trial stations
  real(dp), allocatable     :: tx_xval_train(:,:)    ! transposed station predictor array for xval
  real(dp), allocatable     :: twx_xval_train(:,:)   ! transposed matmul(tx,w) array for kfold xval
  real(dp), allocatable     :: vv(:)                 ! temporary matrix to test for singularities
  real(dp), allocatable     :: B(:)                  ! regression coefficients for predictors
  
  integer(I4B), allocatable :: train_indices(:)      ! vector of training indices
  integer(I4B), allocatable :: test_indices(:)       ! vector of testing indices
  
  integer(I4B)              :: i, j                   ! counter variables
  integer(I4B)              :: nPredictors            ! number of predictors in attribute arrays
  integer(i4b)              :: trial_n_test, trial_n_train ! test & train sample sizes in each trial (varies)

  ! initialize variables
  weightSum = 0.0; errorSum  = 0.0

  ! assign some parameters
  nPredictors = size(X,2)                             ! number of predictors (static + dynamic)
  
  ! allocations
  allocate(x_xval_train(n_train, nPredictors))        ! this is the max trial_n_train could be
  allocate(y_xval_train(n_train))
  allocate(w_xval_train(n_train, n_train))
  allocate(w_1d(n_total))
  allocate(tx_xval_train(nPredictors, n_train))
  allocate(twx_xval_train(nPredictors, n_train))
  allocate(B(nPredictors))
  allocate(vv(nPredictors))
  allocate(train_indices(n_train))
  allocate(test_indices(max_n_test))

  ! initalize weight vector variable
  do i = 1,n_total,1
    w_1d(i) = W(i,i)
  end do

  ! loop through the kfold trials
  do i = 1, kfold_trials
    trial_n_test  = xval_indices(i, 1)       ! number of test records; the last sample may be larger than nominal
    trial_n_train = n_total - trial_n_test   ! number of training records:  last sample may be smaller
    
    ! extract vectors of training and test indices for each trial (for code readability)
    train_indices = xval_indices(i, 2:(trial_n_train+1))
    test_indices  = xval_indices(i, (trial_n_train+2):(trial_n_train+1+trial_n_test))
    
    ! station predictor array, values and weights -- training sample
    x_xval_train = X(train_indices, :)
    y_xval_train = Y(train_indices)
    w_xval_train = W(train_indices, train_indices)
    
    ! create transposed matrices for least squares and non-singular matrix check
    tx_xval_train  = transpose(x_xval_train)
    twx_xval_train = matmul(tx_xval_train, w_xval_train)
    vv = maxval(abs(twx_xval_train), dim=2)

    ! if matrix is singular, default to calculating weighted error around the mean value of the obs test sample (eg variance)
    if(any(vv==0.0)) then 
      meanVal = sum(Y(test_indices(1:trial_n_test)))/trial_n_test
      do j = 1, trial_n_test
        errorSum = errorSum + w_1d(test_indices(j))*(meanVal - Y(test_indices(j)))**2  ! sum up trial *test* sample errors
      end do

    ! otherwise regression can be performed for the current training sample
    else
      call least_squares(x_xval_train, y_xval_train, twx_xval_train, B)         ! regress on training sample
      ! loop over test sample and calculate error
      do j = 1, trial_n_test
        !call least_squares(x_xval, y_xval, twx_xval, B)
        ! regression prediction for test record
        varTmp = real( dot_product( X(test_indices(j), :), B ), kind(sp) )
        ! sum error for tests
        errorSum = errorSum + (w_1d(test_indices(j)) * (varTmp - Y(test_indices(j)))**2)
      end do

    end if
    
    ! sum up normalizing weight for aggregated test samples
    weightSum = weightSum + sum(w_1d(test_indices))        ! sum up trial *test* sample weights
    
  end do
  
  ! normalize and convert error variance into error standard deviation 
  !varUncert = real(errorSum**0.5, kind(sp))/weightSum
  varUncert = real((errorSum/weightSum)**0.5, kind(sp))      ! output cross-validated prediction error
    
end subroutine kfold_crossval
