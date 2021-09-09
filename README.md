# Gridded Meteorological Ensemble Tool (GMET)

# Overview

GMET is a software for created gridded meteorological datasets for precipitation and temperature.  The current applications have been at a daily timestep, yielding daily total precipitation, mean temperature and temperature range. The algorithm is based on locally-weighted spatial regression, applied independently for each day and each output grid cell to a sample of nearby stations.  The approach yields ensemble meteorological fields for which the mean and spread vary in time and space.  

GMET is organized into two programs that are run in sequence: (1) sp_regression, and (2) ens_generation.  The first applies the spatial regression, generating an output file for the domain that contains all the regression coefficients.  The second uses this file to create ensemble members by using spatially correlated random fields to sample the regression uncertainty. 

Several test cases are included (in tar bundles) to give examples of the application of GMET. These are compatible with GMET v2. 
 * Example 1 is for a 2 degree by 2 degree domain and does not use cross-validation or dynamic predictors
 * Example 2 is for a domain including much of northern California, uses cross validation and dynamic predictors (from HRRR).  This test case is described in the paper:        Bunn et al. (2021)
 
GMET has been used successfully in a number of applications, including supporting real-time and retrospective hydrologic modeling and forecasting and research in the western US. 

# Existing GMET Datasets
 * 1/8th degree 100-member NLDAS-domain implementation (Newman et al, 2015), version 2, https://doi.org/10.5065/D6TH8JR2
 * 1 km 100-member Hawaii dataset from 1990-2014, https://doi.org/10.5065/D6SB44JV
 * 2 km 100-member Alaska dataset from 1980-2012, https://doi.org/10.5065/D61Z42T0 
 * 1/8 degree 100-member NLDAS-domain ensemble dressing for NLDASv2 (EDN2), ds613.0 | https://doi.org/10.5065/KARJ-0E19 
 * 1/8th and 1/16th degree 36-member ensembles for a western US domain (125-92W, 25-53N), 1970-2020 (with some real-time components).  Contact: Andy Wood, NCAR

# Notes
 * This code is a work in progress and is provided without guarantee of fitness for any particular application.  
 * The master branch contains the most vetted release
 * The SHARP branch is used operationally to generate real-time ensemble model forcings at NCAR. The only direct pushes to the SHARP branch, however, are made from the operational clone after testing in the operational context.
 * The develop branch will contain up to date general improvements between new release updates in the master branch

# References

## GMET development and major datasets
 * Bunn, PTW, AW Wood, AJ Newman, H Chang, CL Castro, MP Clark and JR Arnold (2021) Improving in situ observation-based ensemble surface meteorological analyses using numerical weather prediction: A case study of the Oroville Dam crisis precipitation event. AMS J. Hydromet. (in review)
 * Newman, A. J. et al. (2020) ‘Probabilistic Spatial Meteorological Estimates for Alaska and the Yukon’, Journal of Geophysical Research: Atmospheres, 125(22), pp. 1–21. doi: 10.1029/2020JD032696.
 * Newman, A. J. et al. (2019) ‘Use of daily station observations to produce high-resolution gridded probabilistic precipitation and temperature time series for the Hawaiian Islands’, Journal of Hydrometeorology, 20(3), pp. 509–529. doi: 10.1175/JHM-D-18-0113.1.
 * Newman, AJ, MP Clark, J Craig, B Nijssen, AW Wood, E Gutmann, N Mizukami, L Brekke, and JR Arnold, 2015, Gridded Ensemble Precipitation and Temperature Estimates for the Contiguous United States, J. Hydromet., doi: http://dx.doi.org/10.1175/JHM-D-15-0026.1
 * Clark, M. P. and Slater, A. G. (2006) ‘Probabilistic Quantitative Precipitation Estimation in Complex Terrain’, Hydrometeorology, Journal O F, (2000), pp. 3–22.

## Applications or Research using GMET
 * Huang, C, AJ Newman, MP Clark, AW Wood and X Zheng, 2016, Evaluation of snow data assimilation using the ensemble Kalman Filter for seasonal streamflow prediction in the Western United States, Hydrol. Earth Syst. Sci. 21, 635-650, 2017, http://www.hydrol-earth-syst-sci.net/21/635/2017/, doi:10.5194/hess-21-635-2017
 * Mendoza, PA, AW Wood, EA Clark, E Rothwell, MP Clark, B Nijssen, LD Brekke, and JR Arnold, 2017, An intercomparison of approaches for improving predictability in operational seasonal streamflow forecasting, Hydrol. Earth Syst. Sci., 21, 3915–3935, 2017
 * Liu, Hongli, AW Wood, AJ Newman and MP Clark, 2021, Ensemble dressing of meteorological fields: using spatial regression to estimate uncertainty in deterministic gridded meteorological datasets, AMS J. Hydromet. (submitted)