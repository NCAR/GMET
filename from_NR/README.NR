The current version of GMET still relies on some code from Numerical Recipes

The 'downscale' program requires: 
   * numrec.f90

The 'scrf' program requires:
   * nrutil.f90  
   * nrtype.f90  
   * nr.f90
   * a modified version of 'ran_state.f90', called 'ran_state_modified_NR.f90' (A.Wood, Jul-2016)
     This modified version generates a fully random seed with each start -- ran_init() is altered
     and depends on a routine in the main scrf directory

   * also       ran1.f90 
                gser.f90 
                gcf.f90 
                gamdev.f90 
                gammln.f90 
                gammp.f90 
                gammq.f90 
                lubksb.f90 
                ludcmp.f90 
                erf.f90 
                erfc.f90 
                erfcc.f90
                gasdev.f90  

