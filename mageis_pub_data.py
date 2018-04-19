# This script will load in the MagEIS internal data, "highrate" from RBSP-A
# and the normal "int" MagEIS data product from RBSP-B, and convert it to 
# a cdf file with time, PA, and channel fluxes. In the flux data, include the
# channel energy ranges.

import sys
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

from spacepy import pycdf

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp/')
import plot_mageis

START_TIME = datetime(2017, 3, 31, 11, 15, 0)
END_TIME = datetime(2017, 3, 31, 11, 20, 0)
SC_ID = 'A'

if SC_ID == 'A':
    fluxConvert = 'True'
else:
    fluxConvert = 'False'

f = plot_mageis.PlotMageis(SC_ID, START_TIME, 'highrate', 
                        tRange=[START_TIME, END_TIME], 
                        instrument='low')
f.getFluxTimeseries(fluxConvert=fluxConvert) # Get flux and times
alphas = f.magEISdata[f.alphaKey][:, :f.n_sectors].flatten() # get alphas

def get_mageis_params(ids):
    """
    This function returns a dictionary of energy and geometric conversion arrays
    that work with the level 3 magEIS data.
    """
    assert ids.upper() in ['A', 'B'], 'ERROR: Incorrect RBSP id!'
        
    # Define magEIS detector constants
    if ids.upper() == 'A':
        Emid = [34, 54, 78, 108, 143, 182, 223] # keV
        Elow = [29, 46, 68, 95, 126, 164, 206] # keV
        Ehigh = [41, 66, 92, 126, 164, 204, 247] # keV
        # Units of (keV cm^2 sr)
        G0dE = [4.13E-2, 5.73E-2, 6.056E-2, 6.88E-2, 7.35E-2, 
            6.90E-2, 5.98E-2]
        
    if ids.upper() == 'B':
        Emid = [32, 51, 74, 101, 132, 168, 208] # keV
        Elow = [27, 43, 63, 88, 117, 152, 193] # keV
        Ehigh = [39, 63, 88, 117, 150, 188] # keV
        # Units of (keV cm^2 sr)
        G0dE = [4.33E-2, 5.41E-2, 5.926E-2, 6.605E-2, 6.460E-2,
            6.23E-2, 5.96E-2]
    return {'Emid':Emid, 'Elow':Elow, 'Ehigh':Ehigh, 'G0dE':G0dE}

### WRITE TO FILE ###
fName = '{}_rbsp{}_mageis_low_highrate.cdf'.format(
                        START_TIME.date(), SC_ID.lower())
with pycdf.CDF(fName, '') as outF:
    outF['Epoch'] = f.times
    outF['FEDU'] = f.flux
    outF['FEDU_Alpha'] = alphas

    # Copy attributes
    outF.attrs = f.magEISdata.attrs
    outF['Epoch'].attrs = f.magEISdata['Epoch'].attrs
    outF['FEDU'].attrs = f.magEISdata['FEDU'].attrs
    outF['FEDU_Alpha'].attrs = f.magEISdata['FEDU_Alpha'].attrs

