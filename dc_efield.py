# This script will plot the DC electric field from the EFW insturment at the 
# time of the microburst.

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

import spacepy.pycdf

tRange = [datetime(2017, 3, 31, 11, 10), datetime(2017, 3, 31, 11, 30)]
sc_id = 'b'
d = spacepy.pycdf.CDF('/home/mike/research/rbsp/data/efw/rbsp{}/'
    'rbsp{}_efw-l2_e-spinfit-mgse_20170331_v02.cdf'.format(sc_id, sc_id))

# Now filter the data by times and valid values
validIdt = (d['epoch'][:] > tRange[0]) & (d['epoch'][:] < tRange[1])

# Calc |E| for the y and z components (E_x is all error values)
Etot = np.sqrt(np.array(d['efield_spinfit_mgse'])[validIdt, 1]**2 + 
    np.array(d['efield_spinfit_mgse'])[validIdt, 2]**2)

# Now do a running average of the electric field.
dt = 11 # Kernel size
avgE = np.convolve(np.ones(dt)/dt, Etot, mode = 'same')

#### E_x is all error values
plt.plot(np.array(d['epoch'])[validIdt], 
    np.array(d['efield_spinfit_mgse'])[validIdt, 1], 
    label=d['efield_mgse_LABL_1'][1])
plt.plot(np.array(d['epoch'])[validIdt], 
    np.array(d['efield_spinfit_mgse'])[validIdt, 2], 
    label=d['efield_mgse_LABL_1'][2])
plt.plot(np.array(d['epoch'])[validIdt], Etot, label='|EyEz|')
plt.plot(np.array(d['epoch'])[validIdt], avgE, 
    label='E_avg over {} s'.format(dt*11))
plt.ylabel('efield_spinfit_mgse (mV/m)')
plt.xlabel('UTC')
plt.title('VAP-{} EFW Electric field'.format(sc_id.upper()))
plt.legend(fontsize = 10)
plt.ylim(-100, 500)
plt.show()
