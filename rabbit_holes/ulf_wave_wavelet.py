# This script looks at the wavelet power spectrum for the EMFISIS magnetometer
# data.
import copy
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

from spacepy import pycdf

sys.path.insert(0, '/home/mike/research/microburst_detector/wavelets')
import waveletAnalysis

# Read in the data (it's big so filter it as well)
d = pycdf.CDF(os.path.join('/home/mike/research/rbsp/data/emfisis/rbspa', 
    'rbsp-a_magnetometer_hires-geo_emfisis-L3_20170331_v1.6.1.cdf'))
# Usefull keys: 'Epoch', 'Magnitude'
t = np.array(d['Epoch'][2592500:2611699]) # Shortcut to filter times.
B = np.array(d['Magnitude'][2592500:2611699])
#tRange = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 20)]
cadence = 1/64
#validIdt = np.where((t > tRange[0]) & (t < tRange[1]))[0]

### WAVELET TRANSFORM ###
w = waveletAnalysis.WaveletDetector(B, t, cadence, mother='DOG', j1=60)
w.waveletTransform()

fig, ax = plt.subplots(2, sharex=True, figsize=(13, 10))
ax[0].plot(t, B)
ax[0].set(ylabel='|B| (nT)', title='VAP-A EMFISIS Magnetometer')

w.plotPower(ax=ax[1])
ax[1].set(xlabel='UTC')
fig.autofmt_xdate()

plt.tight_layout()
plt.savefig('vapa_emfisis_dog_wavelet_spectrum.png')
plt.show()
