Module inputdat2d
  Use nrtype ! variable types (i4b, dp, etc.)
  Use dat_2dgrid ! generic data structure for 2-d grid
  Implicit None
  Save
 ! --------------------------------------------------------------------------------------
 ! data structures for forcing and (calibration, evaluation, assimilation) data
  Type (GENDAT), Pointer :: APRECIP ! precipitation (kg m-2 dt-1)
  Type (GENDAT), Pointer :: STAPREC ! station precipitation data (kg m-2 dt-1)
  Type (GENDAT), Pointer :: AVGTEMP ! average temperature (K)
  Type (GENDAT), Pointer :: REL_HUM ! relative humidity (%)
  Type (GENDAT), Pointer :: SWRNDWN ! downwelling shortwave radiation (W m-2)
  Type (GENDAT), Pointer :: LWRNDWN ! downwelling longwave radiation (W m-2)
  Type (GENDAT), Pointer :: WINDSPD ! windspeed (m/s)
  Type (GENDAT), Pointer :: AIRPRES ! air pressure (hPa)
  Type (GENDAT), Pointer :: STREAMQ ! streamflow (m3 s-1)
  Type (GENDAT), Pointer :: SNWSTOR ! snow water equivalent (mm)
  Type (GENDAT), Pointer :: SNWAREA ! fractional snow covered area (-)
  Type (GENDAT), Pointer :: CLDAREA ! percent of area classified as clouds
  Type (GENDAT), Pointer :: SOILH2O ! soil moisture (mm)
  Type (GENDAT), Pointer :: LAKSTAG ! lake stage (lake level) (m)
  Type (GENDAT), Pointer :: LAKSTRQ ! lake outflow (for managed lakes) (m3 s-1)
  Type (GENDAT), Pointer :: SRFRAIN ! mean precipitation surface (arbitrary units)
 ! --------------------------------------------------------------------------------------
 ! additional metadata for streamflow
  Type SFINFO
  ! station attributes
    Character (Len=120) :: STNNAME ! station name
    Integer (I4B) :: TIDE_ID ! station IDs for flow stns
    Integer (I4B) :: RCH_REC ! REC IDs for flow stns
    Integer (I4B) :: RCH_IDX ! index for flow stns (in network topology)
    Integer (I4B) :: RCH_IDV ! index for valid flow stns (in network topology)
  ! flow statistics
    Real (DP) :: FLOWANN ! mean annual flood
    Real (DP) :: FLOW100 ! 100-year flood
    Real (DP) :: FLOWMAX ! maximum flow ever recorded
  ! data assimilation control parameters
    Logical (LGT) :: STN_USE ! assimilate data from station
    Integer (I4B) :: NUM_AVG ! number of points to average
    Real (DP) :: STN_ERR ! error in station obervations
  ! list of indices of basins upstream of each gauge
    Integer (I4B), Dimension (:), Pointer :: RCH_UPS ! reaches upstream of each gauge
  End Type SFINFO
 ! allow space for multiple gauges in the basin/region
  Type (SFINFO), Dimension (:), Pointer :: SFMETA
 ! --------------------------------------------------------------------------------------
! additional metadata for lake level
  Type LAKEINFO
  ! station attributes
    Character (Len=120) :: STNNAME ! station name
    Integer (I4B) :: TIDE_ID ! station IDs for lake observation locations
    Integer (I4B) :: RCH_REC ! REC reach IDs for lake observation locations
    Integer (I4B) :: RCH_IDX ! index for lake observation locations (in network topology)
    Integer (I4B) :: LAKE_ID ! lake IDs for lake observation locations (in network topology)
    Integer (I4B) :: LAKE_IX ! lake index for lake observation locations (in network topology)
! data assimilation control parameters
    Logical (LGT) :: STN_USE ! assimilate data from station
    Integer (I4B) :: NUM_AVG ! number of points to average
    Real (DP) :: STN_ERR ! error in station obervations
  End Type LAKEINFO
 ! allow space for multiple gauges in the basin/region
  Type (LAKEINFO), Dimension (:), Pointer :: LAKEMETA
 ! --------------------------------------------------------------------------------------
End Module inputdat2d
