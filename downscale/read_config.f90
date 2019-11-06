! modified AWW Dec 2015
subroutine read_config (fname, n, config_names, values)
 
  use string_mod
  implicit none
 
  character (len=*) :: fname
  integer :: n
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


  config_names(1) = "MODE"
  config_names(2) = "START_DATE"
  config_names(3) = "END_DATE"
  config_names(4) = "SITE_LIST"
  config_names(5) = "SITE_VAR"
  config_names(6) = "STATION_VAR"
  config_names(7) = "PERTURBATION"
  config_names(8) = "FORECAST"
  config_names(9) = "NUMBER_VARS"
  config_names(10) = "FILE_VARIABLE"
  config_names(11) = "VARIABLE_NAME"
  config_names(12) = "OUTPUT_FILE"
  config_names(13) = "GRID_LIST"
  config_names(14) = "MAX_DISTANCE"
  config_names(15) = "SITE_VAR_T" !modified AJN Sept 2013
  config_names(16) = "DATA_DIRECTORY" ! AWW-feb2016, free data dir from site list path
  config_names(17) = "STN_START_DATE" ! AWW-apr2016, add station period limits
  config_names(18) = "STN_END_DATE"   !  
  config_names(19) = "GEN_STA_WEIGHTS"   !  
  config_names(20) = "STA_WEIGHT_NAME"   !  
  config_names(21) = "NPREDICT"          !total number of predictors (6+number of NWP predictors)
  config_names(22) = "NWP_VAR_NAMES"      !Variable list for NWP predictors
  config_names(23) = "NWP_INPUT_FILE_LIST"   !list of NWP input files

 
  do i = 1, n, 1
    values (i) = ""
  end do
 
  open (11, file=fname, status='old', access='sequential')
 
  do
    read (11, "(A)", iostat=stat) line
    if (stat < 0) exit
    line = adjustl (line)
    if (line(1:1) .ne. "!" .and. len(trim(line)) /= 0) then
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
