GMET example cases

Example 1 Details:
  - A. Wood, 20180701
  - updated for GMET v2
  - 2 degree domain subset domain:  39 to 41N, 107.5 to 109.5 W
  - 1 month of daily station input data
  To compile:  enter build/ directory and make each program (sp_regression and ens_generation)
  - tested with gnu fortran on the Cheyenne HPC system at NCAR
  To run: link executables from the /bin/ directory into run/ directory
  - run ensemble regression with
    ~> sp_regression.exe config.ens_regr.txt
  - run ensemble generation with
    ~> generate_ensemble.exe namelist.ens_forc.txt

Example 2 Details;
  - A. Wood, 20210903
  - GMET v2
  - capabilties:  kfold_xval, station weighting options, and dynamic predictors (HRRR)
  - California domain subset
  - 1 month of daily station input data and HRRR grids
  To compile:  enter build/ directory and make each program (sp_regression and ens_generation)
  - tested with gnu fortran on the Cheyenne HPC system at NCAR
  To run: link executables from the /bin/ directory into run/ directory
  - run ensemble regression with
    ~> sp_regression.exe config.ens_regr.txt
  - run ensemble generation with
    ~> generate_ensemble.exe namelist.ens_forc.txt
  - data and capabilities described in 
      Bunn, PTW, AW Wood, AJ Newman, H Chang, CL Castro, MP Clark, and JR Arnold, 2021, 
      Improving in situ observation-based ensemble surface meteorological analyses using numerical weather prediction:  
      A case study of the Oroville Dam crisis precipitation event. AMS J. Hydromet. (in review)
