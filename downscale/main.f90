
program precip
  use type
  use strings
  implicit none

  interface
     subroutine get_time_list(startdate, enddate, T)
       use utim
       use type
       character (len = 100), intent(in) :: startdate, enddate
       real(DP), allocatable, intent(out) :: T(:)
     end subroutine get_time_list

     subroutine read_refcst(startdate, enddate, file_var, perturbation, var_name, forecast, V, X, Y, T, error)
      use type
       character (len = 100), intent(in) :: startdate, enddate, file_var, perturbation
       character (len = *), intent(in) :: var_name
       integer(I4B), intent(in) :: forecast
       real(DP), allocatable, intent(out) :: V(:,:), X(:), Y(:)
       real(DP), allocatable, intent(out) :: T(:)
       integer, intent(out) :: error
    end subroutine read_refcst

    subroutine read_station_list(file_name, id, name, lat, lon, alt, sslp_n, sslp_e, n_stations, error)
      use type
      character(len=500), intent(in) :: file_name
      character(len=100), allocatable, intent(out) :: id(:), name(:)
      real(DP), allocatable, intent(out) :: lat(:), lon(:), alt(:), sslp_n(:), sslp_e(:)
      integer(I4B), intent(out) :: n_stations
      integer, intent(out) :: error
    end subroutine read_station_list

    subroutine read_grid_list(file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
      use type
      character(len=500), intent(in) :: file_name
      real(DP), allocatable, intent(out) :: lats(:), lons(:), alts(:), slp_n(:), slp_e(:)
      integer(I4B), intent(out) :: nx, ny
      integer, intent(out) :: error
    end subroutine read_grid_list

    subroutine read_nc_grid(file_name,lat,lon,elev,grad_n,grad_e,mask,nx,ny,error)
      use netcdf
      use type

      character(len=500), intent(in)		:: file_name
      real(DP), allocatable, intent(out)	:: lat(:,:),lon(:,:),elev(:,:),grad_n(:,:),grad_e(:,:),mask(:,:)
      integer(I4B), intent(out)			:: nx,ny
      integer, intent(out)			:: error
    end subroutine read_nc_grid


    subroutine estimate_coefficients(D, nvars, Lats, Lons, Times, stnid, stnlat, stnlon, &
         stnalt, stnvar, site_var, site_list, C, POC, error)
      use type
      real(DP), intent(in) :: D(:,:,:), Lats(:), Lons(:)
      real(DP), intent(in) :: Times(:)
      integer(I4B), intent(in) :: nvars
      character (len = 100), intent(in) :: stnid(:)
      real(DP), intent(in) :: stnlat(:), stnlon(:), stnalt(:)
      character(len=100), intent(in) :: stnvar
      character(len=100), intent(in) :: site_var
      character(len=500), intent(in) :: site_list
      real(DP), allocatable, intent(out) :: C(:,:,:), POC(:,:,:)
      integer, intent(out) :: error
    end subroutine estimate_coefficients

    subroutine save_coefficients(n_vars, var_names, coefs, startdate, enddate, times, &
         site_list, station_var, stnid, stnlat, stnlon, stnalt, forecast, file, error)
      use netcdf
      use type
      character (len = 100), intent(in) :: var_names(:), stnid(:)
      integer, intent(in) :: n_vars, forecast
      character (len = 100), intent(in) :: startdate, enddate, station_var
      character (len = 500), intent(in) :: file, site_list
      real(DP), intent(in) :: stnlat(:), stnlon(:), stnalt(:)
      real(DP), intent(in) :: coefs(:,:,:)
      real(DP), intent(in) :: times(:)
      integer, intent(out) :: error
    end subroutine save_coefficients

    subroutine estimate_precip(X, Z, nsta, ngrid, maxDistance, Times,  &
     stnid, stnvar, site_var, site_var_t, site_list, PCP, POP, PCPERR, &
     tmean,tmean_err,trange,trange_err,mean_autocorr, mean_tp_corr,    &
     y_mean,y_std,y_std_all,y_min,y_max,error,pcp_2,pop_2,pcperr_2, &
      tmean_2,tmean_err_2,trange_2,trange_err_2)
      use type
      real(DP), intent(in) :: X(:,:), Z(:,:)
      real(DP), intent(in) :: maxDistance
      integer(I4B), intent(in) :: nsta, ngrid
      real(DP), intent(in) :: Times(:)
      character (len = 100), intent(in) :: stnid(:)
      character(len=100), intent(in) :: stnvar, site_var, site_var_t
      character(len=500), intent(in) :: site_list
      real(SP), allocatable, intent(out) :: PCP(:,:), POP(:,:), PCPERR(:,:)
      real(SP), allocatable, intent(out) :: tmean(:,:),tmean_err(:,:)   !OLS tmean estimate and error
      real(SP), allocatable, intent(out) :: trange(:,:),trange_err(:,:) !OLS trange estimate and error

      real(SP), allocatable, intent(out) :: PCP_2(:,:), POP_2(:,:), PCPERR_2(:,:)
      real(SP), allocatable, intent(out) :: tmean_2(:,:),tmean_err_2(:,:)   !OLS tmean estimate and error
      real(SP), allocatable, intent(out) :: trange_2(:,:),trange_err_2(:,:) !OLS trange estimate and error


      integer, intent(out) :: error
      real(DP),intent(out) :: mean_autocorr(:)  !mean autocorrelation from all stations over entire time period
      real(DP),intent(out) :: mean_tp_corr(:)  !mean correlation for mean temp and precip

      real(DP), intent(out) :: y_mean(:,:), y_std(:,:), y_std_all(:,:)  !std and mean of transformed time step precipitation
      real(DP), intent(out) :: y_min(:,:), y_max(:,:)  !min,max of normalized time step precip
    end subroutine estimate_precip

    subroutine save_precip(pcp, pop, pcperror, tmean, tmean_err, trange, trange_err, &
   			   nx, ny, grdlat, grdlon, grdalt, Times, mean_autocorr, mean_tp_corr, &
                           y_mean, y_std,y_std_all,y_min,y_max,file, error, &
			 pcp_2, pop_2, pcperror_2, tmean_2, tmean_err_2, trange_2, trange_err_2)
      use netcdf
      use type

      real(SP), intent(in) :: pcp(:,:), pop(:,:), pcperror(:,:)
      real(SP), intent(in) :: tmean(:,:),tmean_err(:,:),trange(:,:),trange_err(:,:)

      real(SP), intent(in) :: pcp_2(:,:), pop_2(:,:), pcperror_2(:,:)
      real(SP), intent(in) :: tmean_2(:,:),tmean_err_2(:,:),trange_2(:,:),trange_err_2(:,:)

      integer(I4B), intent(in) :: nx, ny
      real(DP), intent(in) :: grdlat(:), grdlon(:), grdalt(:)
      real(DP), intent(in) :: Times(:)
      real(DP), intent(in) :: mean_autocorr(:),mean_tp_corr(:)

      real(DP), intent(in) :: y_mean(:,:), y_std(:,:) , y_std_all(:,:)
      real(DP), intent(in) :: y_min(:,:), y_max(:,:)
      character (len = 500), intent(in) :: file
      integer, intent(out) :: error
    end subroutine save_precip

  end interface

  character (len = 100) :: config_file
  integer, parameter :: nconfigs = 15       
  character (len = 500) :: config_names(nconfigs)
  character (len = 500) :: config_values(nconfigs)
  character (len = 500) :: site_list, output_file, output_file2, grid_list
  character (len = 100) :: startdate, enddate, perturbation, station_var, site_var, site_var_t
  character (len = 100), allocatable :: file_var(:), var_name(:)
  character (len = 100), allocatable :: stnid(:), stnname(:)

  character(len=2000)		:: arg  !command line arg for configuration file

  integer :: i, error, n_vars, nfile_var, nvar_name, forecast, mode
  integer :: nstations, lenfile

  real(DP), allocatable :: Y(:,:,:), Vals(:,:), Lats(:), Lons(:)
  real(DP), allocatable :: coefs(:,:,:), prob_coefs(:,:,:)
  real(DP), allocatable :: stnlat(:), stnlon(:), stnalt(:), stn_slp_n(:), stn_slp_e(:)
  real(DP), allocatable :: Times(:)

  real(DP), allocatable :: mean_autocorr(:)  !mean auto correlation for all stations over entire time period
  real(DP), allocatable :: mean_tp_corr(:)   !mean correlation between precip and trange (31-day moving avg anomaly)
  real(DP), allocatable :: y_mean(:,:)
  real(DP), allocatable :: y_std(:,:)
  real(DP), allocatable :: y_std_all(:,:)
  real(DP), allocatable :: y_max(:,:)
  real(DP), allocatable :: y_min(:,:)

  real(DP), allocatable :: lat(:,:)
  real(DP), allocatable :: lon(:,:)
  real(DP), allocatable :: elev(:,:)
  real(DP), allocatable :: grad_n(:,:)
  real(DP), allocatable :: grad_e(:,:)
  real(DP), allocatable :: mask(:,:)


  real(DP), allocatable :: X(:,:), Z(:,:)
  integer :: ngrid

  integer(I4B) :: nx, ny, ntimes
  real(DP) :: maxDistance
  real(DP), allocatable :: grdlat(:), grdlon(:), grdalt(:), grd_slp_n(:), grd_slp_e(:), mask_1d(:)
  real(SP), allocatable :: pcp(:,:), pop(:,:), pcperror(:,:)
  real(SP), allocatable :: tmean(:,:),tmean_err(:,:)
  real(SP), allocatable :: trange(:,:), trange_err(:,:)

  real(SP), allocatable :: pcp_2(:,:), pop_2(:,:), pcperror_2(:,:)
  real(SP), allocatable :: tmean_2(:,:),tmean_err_2(:,:)
  real(SP), allocatable :: trange_2(:,:), trange_err_2(:,:)

!code starts below

!mode 2
   config_file = "config_pnw.txt"
!mode 1
!  config_file = "config_prcp.txt"


!get config_file filename from command line
  i = 0
  do
    call get_command_argument(i,arg)
    if(i .eq. 1) config_file=arg
    if(LEN_TRIM(arg) == 0) EXIT
    i = i + 1
  end do

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
  config_names(15) = "SITE_VAR_T" 

  error = 0
  n_vars = 0
  forecast = -1

  call read_config(config_file, nconfigs, config_names, config_values)
  call value(config_values(1), mode, error)
  if(error /= 0) then
     print *, "ERROR: Failed to read mode from config file."
     return
  endif

  startdate = config_values(2)
  enddate = config_values(3)
  site_list = config_values(4)
  site_var = config_values(5)
  site_var_t = config_values(15)
  station_var = config_values(6)
  output_file = config_values(12)


  if(len(trim(startdate)) == 0 .OR. len(trim(enddate)) == 0 .OR. len(trim(site_list)) == 0 .OR. &
       len(trim(output_file)) == 0) then
     print *, "ERROR: Failed to read in one more more required config parameters. (START_DATE, END_DATE, SITE_LIST or OUTPUT_FILE)"
     return
  endif

  if(mode == 1) then
     ! Model Variables
     perturbation = config_values(7)
     call value(config_values(8), forecast, error)
     if(error /= 0) forecast = -1
     call value(config_values(9), n_vars, error)
     if(error /= 0) n_vars = 0
     
     if(len(trim(perturbation)) == 0 .OR. n_vars == 0 .OR. forecast == -1) then
        print *, "ERROR: Failed to read in one or more required model config parameters. (PERTURBATION, NUMBER_VARS or FORECAST)"
        return
     endif
     
     if(n_vars > 1) then
        allocate(file_var(n_vars))
        allocate(var_name(n_vars))
        call parse(config_values(10), ",", file_var, nfile_var)
        call parse(config_values(11), ",", var_name, nvar_name)
        if(nfile_var /= n_vars .OR. nvar_name /= n_vars) then
           print *, "ERROR: Number of variables in config file does not match."
           return
        endif
     else
        allocate(file_var(1))
        allocate(var_name(1))
        file_var(1) = config_values(10)
        var_name(1) = config_values(11)
        if(len(trim(file_var(1))) == 0 .OR. len(trim(var_name(1))) == 0) then
           print *, "ERROR: Failed to read in one or more required model config parameters. (FILE_VARIABLE or VARIABLE_NAME)"
           return
        endif
        
     endif

     do i = 1, n_vars, 1
        call read_refcst(startdate, enddate, file_var(i), perturbation, var_name(i), &
             forecast, Vals, Lats, Lons, Times, error)
        if(error /= 0) return
        if(i == 1) then
           allocate(Y(n_vars, size(Lons) * size(Lats), size(Times)))
        endif
        Y(i,:,:) = Vals(:,:)
        deallocate(Vals)        
     enddo 
     
     call read_station_list(site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, stn_slp_e, nstations, error)
     if(error /= 0) return

     call estimate_coefficients(Y, n_vars, Lats, Lons, Times, stnid, stnlat, stnlon, &
          stnalt, station_var, site_var, site_list, coefs, prob_coefs, error)
     if(error /= 0) return


     if(trim(station_var) .EQ. "PRCP") then
        lenfile=len_trim(output_file)
        output_file2(:) = " "
        output_file2(1:5) = "prob_"
        output_file2(6:lenfile+6) = output_file(1:lenfile)

        call save_coefficients(n_vars, var_name, prob_coefs, startdate, enddate, Times, &
             site_list, station_var, stnid, stnlat, stnlon, stnalt, forecast, output_file2, error)
        if(error /= 0) return
     end if
     
     call save_coefficients(n_vars, var_name, coefs, startdate, enddate, Times, &
          site_list, station_var, stnid, stnlat, stnlon, stnalt, forecast, output_file, error)
     if(error /= 0) return

  else
     if(mode == 2) then

        grid_list = config_values(13)
        call value(config_values(14), maxDistance, error)
        maxDistance = maxDistance * 0.539957
        if(error /= 0) maxDistance = -1

        if(len(trim(grid_list)) == 0) then
           print *, "ERROR: Failed to read in one more more required model config parameters. (GRID_LIST)"
           return
        endif

	
        call read_station_list(site_list, stnid, stnname, stnlat, stnlon, stnalt, stn_slp_n, stn_slp_e, nstations, error)
        if(error /= 0) return

	call read_nc_grid(grid_list,lat,lon,elev,grad_n,grad_e,mask,nx,ny,error)

      !allocate 1-d grid variables
      allocate(grdlat(nx*ny))
      allocate(grdlon(nx*ny))
      allocate(grdalt(nx*ny))
      allocate(grd_slp_n(nx*ny))
      allocate(grd_slp_e(nx*ny))
      allocate(mask_1d(nx*ny))


	grdlat = reshape(lat,(/nx*ny/))
	grdlon = reshape(lon,(/nx*ny/))
	grdalt = reshape(elev,(/nx*ny/))
	grd_slp_n = reshape(grad_n,(/nx*ny/))
	grd_slp_e = reshape(grad_e,(/nx*ny/))
	mask_1d = reshape(mask,(/nx*ny/))

        ngrid = nx*ny
        allocate(X(nstations,6))
        allocate(Z(ngrid,6))
        X(:,1) = 1.0
        X(:,2) = stnlat(:)
        X(:,3) = stnlon(:)
        X(:,4) = stnalt(:)
	X(:,5) = stn_slp_n(:)
	X(:,6) = stn_slp_e(:)

        Z(:,1) = 1.0
        Z(:,2) = grdlat(:)
        Z(:,3) = grdlon(:)
        Z(:,4) = grdalt(:)
	Z(:,5) = grd_slp_n(:)
	Z(:,6) = grd_slp_e(:)


        call get_time_list(startdate, enddate, Times)
        
	ntimes = size(Times)

	allocate(mean_autocorr(ntimes))
	allocate(mean_tp_corr(ntimes))
        allocate(y_mean(ngrid,ntimes))
        allocate(y_std(ngrid,ntimes))
	allocate(y_std_all(ngrid,ntimes))
        allocate(y_max(ngrid,ntimes))
        allocate(y_min(ngrid,ntimes))

        call estimate_precip(X, Z, nstations, ngrid, maxDistance, Times,  &
			    stnid, station_var, site_var, site_var_t, site_list, pcp, pop, pcperror, &
			    tmean, tmean_err, trange,trange_err, mean_autocorr, mean_tp_corr, &
                            y_mean,y_std,y_std_all,y_min,y_max,error, &
			    pcp_2, pop_2, pcperror_2,tmean_2, tmean_err_2, trange_2,trange_err_2)
        if(error /= 0) return


        print *,'Creating output file'

        call save_precip(pcp, pop, pcperror, tmean, tmean_err, trange, trange_err, &
                         nx, ny, grdlat, grdlon, grdalt, Times, mean_autocorr, mean_tp_corr, &
                         y_mean,y_std,y_std_all,y_min,y_max,output_file, error, &
			 pcp_2, pop_2, pcperror_2,tmean_2, tmean_err_2, trange_2,trange_err_2)
        if(error /= 0) return


     endif
  endif


end program precip



subroutine get_time_list(startdate, enddate, Times)
  use utim
  use type
  implicit none

  character (len = 100), intent(in) :: startdate, enddate
  real(DP), allocatable, intent(out) :: Times(:)

  integer(I4B) :: t, ntimes, sday, eday, error
  integer(I4B) :: sec, min, hour, day, month, year
  real(DP) :: utime  

  call parse_date(startdate, year, month, day, hour, min, sec, error)
  sday = julian_date(day, month, year)
  call parse_date(enddate, year, month, day, hour, min, sec, error)
  eday = julian_date(day, month, year)
  ntimes = eday - sday +1
  allocate(Times(ntimes))


  utime= date_to_unix(startdate)
  do t = 1, ntimes, 1
     if(utime > date_to_unix(enddate)) exit
     Times(t) = utime
     utime = utime+86400
  enddo

end subroutine get_time_list
