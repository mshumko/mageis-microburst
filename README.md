Project owner: Mykhaylo Shumko
Institutions: Montana State University and The Aerospace Corporation.

This project code and a small amount of data is for analyzing microbursts observed with the Van Allen Probes.

Programs:
find_psd.py - Use the MagEphem data to look for times when RBSP was within
    a user-defined L and MLT range. The top panel shows B_l/B_sc as an 
    indicator of how far RBSP is off of the equator.

figure_scripts/
- fig1_plot.py - This script plots figure 1 for the microburst paper.
        It is pretty dissorganized since I was adding in plotting 
        functionality here and there.
- fig2_plot.py - This scipt plots the MagEIS and RBSPICE data for figure
        2 of the microburst paper.
- fig3_plot.py - This scipt plots the AC6 data during the RBSP-AC6 
        conjunction. The time and distance lags are hard coded since the 
        data values are errorous.

rabbit_holes/
    dc_efield.py - Plot the EFW instrument's y and z compoents of the E field.
        This script will also plot |E| and do a running average of |E|.        
        NOTE: RBSP-A's EFW instrument is noisy and should not be trusted!
    calc_sc_sep.py - Calculates the eclidian seperation between the RBSP
        spacecraft.
    mirror_alt.py - This script calculates the morror point of particles that 
        have pitch angle alpha at position X during Kp given in maginut 
        (T89 model).
    numerical_bounce_period.py - This script calculates the electron 
        bounce periods in a 2D grid defined by local pitch angles and
        kinetic energies.
    analytic_bounce_period.py - Similar style script as numerical_bounce_period
        but uses the Schulz & Lanzerotti derived analytic approach.
    energy_spectra.py - Calculates the energy spectra of the microbursts 
        assuming they fall off exponentially. The MagEIS flux can be 
        background subtracted first. NOTE: This code has not been 
        refactored to the new MagEIS plotting code!
    url_wave_wavelet.py - This script will load the EMFISIS hires 
        magnetometer data to calculate a wavelet power spectrum using a DOG
        m = 2 wavelet to look for ULF waves. NOTE: The file and indicies are 
        hard coded since the filtering takes so long with these files.
        
The code contains:
- diffusion and resonance curve analysis, 
- PAD fitting scripts
- scripts to generate figures for the manuscript.

This code will contain:
- Scripts to calculate the electron bounce period
