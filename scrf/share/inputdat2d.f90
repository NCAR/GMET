MODULE inputdat2d
  USE nrtype ! variable types (i4b, dp, etc.)
  USE dat_2dgrid ! generic data structure for 2-d grid
  IMPLICIT NONE
  SAVE
 ! --------------------------------------------------------------------------------------
 ! data structures for forcing and (calibration, evaluation, assimilation) data
  TYPE (GENDAT), POINTER :: APRECIP ! precipitation (kg m-2 dt-1)
  TYPE (GENDAT), POINTER :: STAPREC ! station precipitation data (kg m-2 dt-1)
  TYPE (GENDAT), POINTER :: AVGTEMP ! average temperature (K)
  TYPE (GENDAT), POINTER :: REL_HUM ! relative humidity (%)
  TYPE (GENDAT), POINTER :: SWRNDWN ! downwelling shortwave radiation (W m-2)
  TYPE (GENDAT), POINTER :: LWRNDWN ! downwelling longwave radiation (W m-2)
  TYPE (GENDAT), POINTER :: WINDSPD ! windspeed (m/s)
  TYPE (GENDAT), POINTER :: AIRPRES ! air pressure (hPa)
  TYPE (GENDAT), POINTER :: STREAMQ ! streamflow (m3 s-1)
  TYPE (GENDAT), POINTER :: SNWSTOR ! snow water equivalent (mm)
  TYPE (GENDAT), POINTER :: SNWAREA ! fractional snow covered area (-)
  TYPE (GENDAT), POINTER :: CLDAREA ! percent of area classified as clouds
  TYPE (GENDAT), POINTER :: SOILH2O ! soil moisture (mm)
  TYPE (GENDAT), POINTER :: LAKSTAG ! lake stage (lake level) (m)
  TYPE (GENDAT), POINTER :: LAKSTRQ ! lake outflow (for managed lakes) (m3 s-1)
  TYPE (GENDAT), POINTER :: SRFRAIN ! mean precipitation surface (arbitrary units)
 ! --------------------------------------------------------------------------------------
 ! additional metadata for streamflow
  TYPE SFINFO
  ! station attributes
   CHARACTER (LEN=120) :: STNNAME ! station name
   INTEGER (I4B) :: TIDE_ID ! station IDs for flow stns
   INTEGER (I4B) :: RCH_REC ! REC IDs for flow stns
   INTEGER (I4B) :: RCH_IDX ! index for flow stns (in network topology)
   INTEGER (I4B) :: RCH_IDV ! index for valid flow stns (in network topology)
  ! flow statistics
   REAL (DP) :: FLOWANN ! mean annual flood
   REAL (DP) :: FLOW100 ! 100-year flood
   REAL (DP) :: FLOWMAX ! maximum flow ever recorded
  ! data assimilation control parameters
   LOGICAL (LGT) :: STN_USE ! assimilate data from station
   INTEGER (I4B) :: NUM_AVG ! number of points to average
   REAL (DP) :: STN_ERR ! error in station obervations
  ! list of indices of basins upstream of each gauge
   INTEGER (I4B), DIMENSION (:), POINTER :: RCH_UPS ! reaches upstream of each gauge
  END TYPE SFINFO
 ! allow space for multiple gauges in the basin/region
  TYPE (SFINFO), DIMENSION (:), POINTER :: SFMETA
 ! --------------------------------------------------------------------------------------
! additional metadata for lake level
  TYPE LAKEINFO
  ! station attributes
   CHARACTER (LEN=120) :: STNNAME ! station name
   INTEGER (I4B) :: TIDE_ID ! station IDs for lake observation locations
   INTEGER (I4B) :: RCH_REC ! REC reach IDs for lake observation locations
   INTEGER (I4B) :: RCH_IDX ! index for lake observation locations (in network topology)
   INTEGER (I4B) :: LAKE_ID ! lake IDs for lake observation locations (in network topology)
   INTEGER (I4B) :: LAKE_IX ! lake index for lake observation locations (in network topology)
! data assimilation control parameters
   LOGICAL (LGT) :: STN_USE ! assimilate data from station
   INTEGER (I4B) :: NUM_AVG ! number of points to average
   REAL (DP) :: STN_ERR ! error in station obervations
  END TYPE LAKEINFO
 ! allow space for multiple gauges in the basin/region
  TYPE (LAKEINFO), DIMENSION (:), POINTER :: LAKEMETA
 ! --------------------------------------------------------------------------------------
END MODULE inputdat2d
