module inputdat2d
  use nrtype ! variable types (i4b, dp, etc.)
  use dat_2dgrid ! generic data structure for 2-d grid
  implicit none
  save
 ! --------------------------------------------------------------------------------------
 ! data structures for forcing and (calibration, evaluation, assimilation) data
  type (gendat), pointer :: aprecip ! precipitation (kg m-2 dt-1)
  type (gendat), pointer :: staprec ! station precipitation data (kg m-2 dt-1)
  type (gendat), pointer :: avgtemp ! average temperature (K)
  type (gendat), pointer :: rel_hum ! relative humidity (%)
  type (gendat), pointer :: swrndwn ! downwelling shortwave radiation (W m-2)
  type (gendat), pointer :: lwrndwn ! downwelling longwave radiation (W m-2)
  type (gendat), pointer :: windspd ! windspeed (m/s)
  type (gendat), pointer :: airpres ! air pressure (hPa)
  type (gendat), pointer :: streamq ! streamflow (m3 s-1)
  type (gendat), pointer :: snwstor ! snow water equivalent (mm)
  type (gendat), pointer :: snwarea ! fractional snow covered area (-)
  type (gendat), pointer :: cldarea ! percent of area classified as clouds
  type (gendat), pointer :: soilh2o ! soil moisture (mm)
  type (gendat), pointer :: lakstag ! lake stage (lake level) (m)
  type (gendat), pointer :: lakstrq ! lake outflow (for managed lakes) (m3 s-1)
  type (gendat), pointer :: srfrain ! mean precipitation surface (arbitrary units)
 ! --------------------------------------------------------------------------------------
 ! additional metadata for streamflow
  type sfinfo
  ! station attributes
    character (len=120) :: stnname ! station name
    integer (i4b) :: tide_id ! station IDs for flow stns
    integer (i4b) :: rch_rec ! REC IDs for flow stns
    integer (i4b) :: rch_idx ! index for flow stns (in network topology)
    integer (i4b) :: rch_idv ! index for valid flow stns (in network topology)
  ! flow statistics
    real (dp) :: flowann ! mean annual flood
    real (dp) :: flow100 ! 100-year flood
    real (dp) :: flowmax ! maximum flow ever recorded
  ! data assimilation control parameters
    logical (lgt) :: stn_use ! assimilate data from station
    integer (i4b) :: num_avg ! number of points to average
    real (dp) :: stn_err ! error in station obervations
  ! list of indices of basins upstream of each gauge
    integer (i4b), dimension (:), pointer :: rch_ups ! reaches upstream of each gauge
  end type sfinfo
 ! allow space for multiple gauges in the basin/region
  type (sfinfo), dimension (:), pointer :: sfmeta
 ! --------------------------------------------------------------------------------------
! additional metadata for lake level
  type lakeinfo
  ! station attributes
    character (len=120) :: stnname ! station name
    integer (i4b) :: tide_id ! station IDs for lake observation locations
    integer (i4b) :: rch_rec ! REC reach IDs for lake observation locations
    integer (i4b) :: rch_idx ! index for lake observation locations (in network topology)
    integer (i4b) :: lake_id ! lake IDs for lake observation locations (in network topology)
    integer (i4b) :: lake_ix ! lake index for lake observation locations (in network topology)
! data assimilation control parameters
    logical (lgt) :: stn_use ! assimilate data from station
    integer (i4b) :: num_avg ! number of points to average
    real (dp) :: stn_err ! error in station obervations
  end type lakeinfo
 ! allow space for multiple gauges in the basin/region
  type (lakeinfo), dimension (:), pointer :: lakemeta
 ! --------------------------------------------------------------------------------------
end module inputdat2d
