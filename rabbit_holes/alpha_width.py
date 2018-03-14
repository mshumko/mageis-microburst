# This script will make scatter plots of pitch angle vs time for the RBSPICE
# EBR data, and MagEIS data.

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
#import matplotlib.ticker as ticker
#import matplotlib.colors
from datetime import datetime, timedelta

import spacepy.pycdf

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp/')
import plot_mageis
#import plot_rbspice
#import fig1_plot

sys.path.insert(0, '/home/mike/research/mageis-microburst/rabbit_holes')
import fit_rbspice_pa
import rbspice_alpha_fit_popt


burstTimes = np.array([datetime(2017, 3, 31, 11, 17, 9, 777000),
                    datetime(2017, 3, 31, 11, 17, 10, 250000),
                    datetime(2017, 3, 31, 11, 17, 10, 860000),
                    datetime(2017, 3, 31, 11, 17, 11, 500000)])

# Load MagEIS-LOW data for pitch angle information
mageisObj = plot_mageis.PlotMageis('A', datetime(2017, 3, 31), 'highrate', 
    tRange=None, instrument='low')
mageisT = mageisObj._resolveSpinTimes(True) # Get falttened timestamps
alphas = mageisObj.magEISdata['HighRate_Alpha'].flatten()
print(len(mageisT), len(alphas))

# Load RBSPICE EBR data for count rate
d = spacepy.pycdf.CDF('/home/mike/research/rbsp/data/rbspice/rbspa/'
                        'rbsp-a-rbspice_lev-1_EBR_20170331_v1.1.2-01')
rbspiceT = np.asarray(d['Epoch'])

# Load parameters to convert time to PA for each telescope
popt = rbspice_alpha_fit_popt.popt
t0 = rbspice_alpha_fit_popt.t0

for t in burstTimes: 
    # Find RBSPICE index for each time in burstTimes.
    idT = np.where(rbspiceT > t)[0][0]
    EBR = np.array(d['EBR'][:][idT])
    a = [None]*6
    
    for tel in range(5):
        # Calculate pitch angle for each telescope at time t.
        a[tel] = fit_rbspice_pa.sin_fit(
            mdates.date2num(rbspiceT[idT])-t0, *popt[tel])
    # Plot results
    plt.scatter(a, d['EBR'][:][idT])
    plt.xlabel(r'$\alpha_L$ [deg]')
    plt.ylabel('RBSPICE EBR [counts/s]')
    plt.title('RBSPICE and MagEIS PA at {}'.format(t))
    plt.savefig('{}_rbspice_alpha.png'.format(t))
    plt.close()












