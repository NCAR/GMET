GMET example case
A. Wood, 20180701

Example 1 Details:
  - 2 degree domain subset domain:  39 to 41N, 107.5 to 109.5 W
  - 1 month of daily station input data
  - link executables for 'downscale/' and 'scrf/' programs into run/ directory
  - run ensemble regression with
    ~> downscale.exe config.ens_regr.txt
  - run ensemble generation with
    ~> generate_ensemble.exe namelist.ens_forc.txt

Example2 Details;
  - kfold_xval and NWP (HRRR) version merged into SHARP
  - California domain subset
  - 1 month of daily station input data
  - HRRR data for same month (first 15 days)
  - link executables for 'downscale/' and 'scrf/' programs into run/ directory if desired
  - run ensemble regression with
    ~> downscale.exe config.ens_regr.txt
  - run ensemble generation with
    ~> generate_ensemble.exe namelist.ens_forc.txt

