# This script will make scatter plots of pitch angle vs time for the RBSPICE
# EBR data during and before the microburst times.

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

# Load RBSPICE EBR data for count rate
d = spacepy.pycdf.CDF('/home/mike/research/rbsp/data/rbspice/rbspa/'
                        'rbsp-a-rbspice_lev-1_EBR_20170331_v1.1.2-01')
rbspiceT = np.asarray(d['Epoch'])

# Load parameters to convert time to PA for each telescope
popt = rbspice_alpha_fit_popt.popt
t0 = rbspice_alpha_fit_popt.t0

c = ['r', 'b', 'g', 'c', 'k']

for (i, t) in enumerate(burstTimes): 
    # Find RBSPICE index for each time in burstTimes.
    idT = np.where(rbspiceT > t)[0][0]
    EBR = np.array(d['EBR'][:][idT])
    a = [None]*5
    
    for tel in range(5):
        # Calculate pitch angle for each telescope at time t.
        a[tel] = fit_rbspice_pa.sin_fit(
            mdates.date2num(rbspiceT[idT])-t0, *popt[tel])
            
    # Sort the pitch angles
    a, counts = zip(*sorted(zip(a, d['EBR'][idT][:-1])))
    # Calculate pitch angle and count rate one spin prior to t.
    qIdT = np.where(rbspiceT > t - timedelta(seconds=11.16))[0][0]
    qEBR = np.array(d['EBR'][:][qIdT])
    qa = [None]*6

    for tel in range(5):
        # Calculate pitch angle for each telescope one spin prior to t.
        qa[tel] = fit_rbspice_pa.sin_fit(
            mdates.date2num(rbspiceT[qIdT])-t0, *popt[tel])
    qa, qCounts = zip(*sorted(zip(qa, d['EBR'][qIdT][:-1])))
    
    # Plot results
    ratio = np.divide(list(counts), list(qCounts))
    plt.errorbar(qa, ratio, yerr=None, fmt='o--',
                color=c[i], label='Microburst #{}'.format(i+1))

plt.xlim(90, 180)
plt.xlabel(r'$\alpha_L$ [deg]')
plt.ylabel('RBSPICE EBR Ratio')
plt.title('RBSPICE microburst/quiet ratio')# {}\n{}'.format(i+1, t))


plt.legend()
plt.savefig('/home/mike/Dropbox/0_grad_work/'
                'mageis_microburst/plots/RBSPICE/'
                'rbspice_pad_ratio.png')
plt.close()