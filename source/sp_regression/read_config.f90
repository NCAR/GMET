! modified AWW Dec 2015
subroutine read_config (fname, n, config_names, values)
 
  use string_mod
  implicit none
 
  character (len=*)   :: fname
  integer             :: n
  character (len=500) :: config_names (n)
  character (len=500) :: values (n)
  character (len=500) :: settings (2)
  character (1000) line
  integer ipos, stat, nsettings, i
  logical fexist
 
  inquire (file=fname, exist=fexist)
  if ( .not. fexist) then
    print *, "Cannot find config file: ", trim (fname)
    return
  end if

  ! define config file names
  config_names(1)  = "MODE"                ! only mode 2 is currently supported
  config_names(2)  = "START_DATE"          ! start/end dates of regression
  config_names(3)  = "END_DATE"
  config_names(4)  = "SITE_LIST"           ! list of station sites to process
  config_names(5)  = "SITE_VAR"
  config_names(6)  = "STATION_VAR"
  config_names(7)  = "PERTURBATION"
  config_names(8)  = "FORECAST"
  config_names(9)  = "NUMBER_VARS"
  config_names(10) = "FILE_VARIABLE"
  config_names(11) = "VARIABLE_NAME"        ! 
  config_names(12) = "OUTPUT_FILE"          ! output regression file
  config_names(13) = "GRID_LIST"            ! list of grid files needed for ens. generation
  config_names(14) = "MAX_DISTANCE"         ! entered in km
  config_names(15) = "SITE_VAR_T"           ! whether to use cross-correlation with temperature
  config_names(16) = "DATA_DIRECTORY"       ! free data dir from site list path
  config_names(17) = "STN_START_DATE"       ! add station data period limits (start/end)
  config_names(18) = "STN_END_DATE"         !  
  config_names(19) = "GEN_STA_WEIGHTS"      ! T/F whether station weights need to be generated versus read from file
  config_names(20) = "STA_WEIGHT_NAME"      ! filename for binary station weights file
  config_names(21) = "USE_STN_WEIGHTS"      ! use station weights in forming regression (TRUE/FALSE)
  config_names(22) = "N_DYN_PREDICTORS"     ! total number of dynamic predictors
  config_names(23) = "DYN_VAR_NAMES"        ! variable list for dynamic predictors
  config_names(24) = "DYN_PRCP_VAR_NAME"    ! variable name of dynamic precipitation predictor
  config_names(25) = "DYN_INPUT_FILE_LIST"  ! list of dynamic predictor input files
  config_names(26) = "NUM_STATIONS"         ! number of stations to include in regression for each point
  config_names(27) = "KFOLD_TRIALS"         ! number of kfold xval trials to run [2-50]; 0 means cross-val is not run

  ! initialize entries
  do i = 1, n, 1
    values (i) = ""
  end do
 
  open (11, file=fname, status='old', access='sequential')
 
  do
    read (11, "(A)", iostat=stat) line
    if (stat < 0) exit
    line = adjustl (line)
    if (line(1:1) .ne. "!" .and. line(1:1) .ne. "#" .and. len(trim(line)) /= 0) then
      ipos = index (line, '=')
      if (ipos > 0) then
        call parse (line, "=", settings, nsettings)
        do i = 1, n, 1
          if (index(settings(1), config_names(i)) == 1) values (i) = settings (2)
        end do
      end if
    end if
  end do
 
end subroutine read_config
