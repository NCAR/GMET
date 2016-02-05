Module dat_2dgrid
  Use nrtype ! variable types, etc.
  Use linkstruct ! linkage structures -- GRD2BAS
  Implicit None
  Save
 ! --------------------------------------------------------------------------------------
 ! structure to hold data for the model simulation
  Type GENDAT
    Type (DES_IX) :: IDX ! indices of data in input files
    Logical (LGT) :: INIT ! .TRUE. if first call (must read space-time data)
    Logical (LGT) :: LFLG ! .TRUE. if have read linkage information already
    Logical (LGT) :: BEX ! .TRUE. if time bounds exist
    Integer (I4B) :: IDAT0 ! index in TIM associated with first element in data array
    Real (DP), Dimension (:), Pointer :: TB0 ! start of time interval (in data file)
    Real (DP), Dimension (:), Pointer :: TB1 ! end of time interval (in data file)
    Real (DP), Dimension (:), Pointer :: TIM ! time stamp (in data file)
    Real (DP), Dimension (:, :), Pointer :: LAT ! latitude
    Real (DP), Dimension (:, :), Pointer :: LON ! longitude
    Real (DP), Dimension (:, :), Pointer :: ELV ! elevation
    Type (INTERP), Dimension (:), Pointer :: GRD2BAS ! grid-basin relationship
    Real (DP), Dimension (:, :), Pointer :: CEA_DAT ! calib/eval/assim data (time res of observation)
    Real (DP), Dimension (:, :), Pointer :: ANC_DAT ! ancillary data
    Real (DP), Dimension (:), Pointer :: ERROR_CEA ! Assumed uncertainty in CEA_DAT
    Logical (LGT), Dimension (:), Pointer :: MASK_CEA ! True if CEA_DAT is used in sub-basin (otherwise FALSE)
    Real (DP), Dimension (:, :, :), Pointer :: RAW_DAT ! raw data (time res of the data)
    Real (DP), Dimension (:, :, :, :), Pointer :: ENS_DAT ! ensemble data (time res of the data)
    Real (DP), Dimension (:, :, :, :), Pointer :: ENS_ONE ! ensemble data (time res of the data)
    Real (DP), Dimension (:, :, :, :), Pointer :: ENS_TWO ! ensemble data (time res of the data)
    Real (DP), Dimension (:, :, :, :), Pointer :: SIM_DAT ! forcing data (time res of the model sim)
    Character (Len=120) :: TSERIES_FILE ! name of time series file
    Character (Len=120) :: STANDARD_NAME ! standard name of variable
    Character (Len=120) :: CELL_METHOD ! cell method for variable
  End Type GENDAT
 ! --------------------------------------------------------------------------------------
End Module dat_2dgrid
