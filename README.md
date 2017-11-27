Project owner: Mykhaylo Shumko
Institution: Montana State University and The Aerospace Corporation.

This project code and a small amount of data is for analyzing microbursts observed with the Van Allen Probes.

Programs:
find_psd.py - Use the MagEphem data to look for times when RBSP was within
    a user-defined L and MLT range. The top panel shows B_l/B_sc as an 
    indicator of how far RBSP is off of the equator.

rabbit_holes/
    dc_efield.py - Plot the EFW instrument's y and z compoents of the E field.
        This script will also plot |E| and do a running average of |E|.        
        NOTE: RBSP-A's EFW instrument is noisy and should not be trusted!
    calc_sc_sep.py - Calculates the eclidian seperation between the RBSP
        spacecraft.
        

The code contains:
- diffusion and resonance curve analysis, 
- PAD fitting scripts
- scripts to generate figures for the manuscript.

This code will contain:
- Scripts to calculate the electron bounce period
