MODULE dat_2dgrid
 USE nrtype                                               ! variable types, etc.
 USE linkstruct                                           ! linkage structures -- GRD2BAS
 IMPLICIT NONE
 SAVE
 ! --------------------------------------------------------------------------------------
 ! structure to hold data for the model simulation
 TYPE GENDAT
  TYPE(DES_IX)                             :: IDX         ! indices of data in input files
  LOGICAL(LGT)                             :: INIT        ! .TRUE. if first call (must read space-time data)
  LOGICAL(LGT)                             :: LFLG        ! .TRUE. if have read linkage information already
  LOGICAL(LGT)                             :: BEX         ! .TRUE. if time bounds exist
  INTEGER(I4B)                             :: IDAT0       ! index in TIM associated with first element in data array
  REAL(DP), DIMENSION(:), POINTER          :: TB0         ! start of time interval (in data file)
  REAL(DP), DIMENSION(:), POINTER          :: TB1         ! end of time interval (in data file)
  REAL(DP), DIMENSION(:), POINTER          :: TIM         ! time stamp (in data file)
  REAL(DP), DIMENSION(:,:), POINTER        :: LAT         ! latitude
  REAL(DP), DIMENSION(:,:), POINTER        :: LON         ! longitude
  REAL(DP), DIMENSION(:,:), POINTER        :: ELV         ! elevation
  TYPE(INTERP), DIMENSION(:), POINTER      :: GRD2BAS     ! grid-basin relationship
  REAL(DP), DIMENSION(:,:), POINTER        :: CEA_DAT     ! calib/eval/assim data (time res of observation)
  REAL(DP), DIMENSION(:,:), POINTER        :: ANC_DAT     ! ancillary data
  REAL(DP), DIMENSION(:), POINTER          :: ERROR_CEA   ! Assumed uncertainty in CEA_DAT
  LOGICAL(LGT), DIMENSION(:), POINTER      :: MASK_CEA    ! True if CEA_DAT is used in sub-basin (otherwise FALSE)
  REAL(DP), DIMENSION(:,:,:), POINTER      :: RAW_DAT     ! raw data (time res of the data)
  REAL(DP), DIMENSION(:,:,:,:), POINTER    :: ENS_DAT     ! ensemble data (time res of the data)
  REAL(DP), DIMENSION(:,:,:,:), POINTER    :: ENS_ONE     ! ensemble data (time res of the data)
  REAL(DP), DIMENSION(:,:,:,:), POINTER    :: ENS_TWO     ! ensemble data (time res of the data)
  REAL(DP), DIMENSION(:,:,:,:), POINTER    :: SIM_DAT     ! forcing data (time res of the model sim)
  CHARACTER(LEN=120)                       :: TSERIES_FILE ! name of time series file
  CHARACTER(LEN=120)                       :: STANDARD_NAME ! standard name of variable
  CHARACTER(LEN=120)                       :: CELL_METHOD   ! cell method for variable
 END TYPE GENDAT
 ! --------------------------------------------------------------------------------------
END MODULE dat_2dgrid
