module dat_2dgrid
  use nrtype ! variable types, etc.
  use linkstruct ! linkage structures -- GRD2BAS
  implicit none
  save
 ! --------------------------------------------------------------------------------------
 ! structure to hold data for the model simulation
  type gendat
    type (des_ix) :: idx ! indices of data in input files
    logical (lgt) :: init ! .TRUE. if first call (must read space-time data)
    logical (lgt) :: lflg ! .TRUE. if have read linkage information already
    logical (lgt) :: bex ! .TRUE. if time bounds exist
    integer (i4b) :: idat0 ! index in TIM associated with first element in data array
    real (dp), dimension (:), pointer :: tb0 ! start of time interval (in data file)
    real (dp), dimension (:), pointer :: tb1 ! end of time interval (in data file)
    real (dp), dimension (:), pointer :: tim ! time stamp (in data file)
    real (dp), dimension (:, :), pointer :: lat ! latitude
    real (dp), dimension (:, :), pointer :: lon ! longitude
    real (dp), dimension (:, :), pointer :: elv ! elevation
    type (interp), dimension (:), pointer :: grd2bas ! grid-basin relationship
    real (dp), dimension (:, :), pointer :: cea_dat ! calib/eval/assim data (time res of observation)
    real (dp), dimension (:, :), pointer :: anc_dat ! ancillary data
    real (dp), dimension (:), pointer :: error_cea ! Assumed uncertainty in CEA_DAT
    logical (lgt), dimension (:), pointer :: mask_cea ! True if CEA_DAT is used in sub-basin (otherwise FALSE)
    real (dp), dimension (:, :, :), pointer :: raw_dat ! raw data (time res of the data)
    real (dp), dimension (:, :, :, :), pointer :: ens_dat ! ensemble data (time res of the data)
    real (dp), dimension (:, :, :, :), pointer :: ens_one ! ensemble data (time res of the data)
    real (dp), dimension (:, :, :, :), pointer :: ens_two ! ensemble data (time res of the data)
    real (dp), dimension (:, :, :, :), pointer :: sim_dat ! forcing data (time res of the model sim)
    character (len=120) :: tseries_file ! name of time series file
    character (len=120) :: standard_name ! standard name of variable
    character (len=120) :: cell_method ! cell method for variable
  end type gendat
 ! --------------------------------------------------------------------------------------
end module dat_2dgrid
