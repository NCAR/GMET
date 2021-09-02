# GMET
Gridded Meteorological Ensemble Tool

# Overview

GMET is a software for created gridded meteorological datasets for precipitation and temperature.  The current applications have been at a daily timestep, yielding daily total precipitation, mean temperature and temperature range. The algorithm is based on locally-weighted spatial regression, applied independently for each day and each output grid cell to a sample of nearby stations.  The approach yields ensemble meteorological fields for which the mean and spread vary in time and space.  

GMET is organized into two programs:  (1) gmet_regression, and (2) gmet_ensemble

# Notes

* This code is a work in progress and is provided without guarantee of fitness for any particular application.  
* The SHARP branch is used operationally to generate real-time ensemble model forcings at NCAR.  It is likely to be the most 
up-to-date and complete.  The only direct pushes to the SHARP branch, however, are made from the operational clone after 
testing in the operational context.
