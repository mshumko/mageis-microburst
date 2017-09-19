import numpy as np
import os
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime
import dateutil.parser

import spacepy.datamodel

import mageis_diffusion_curves
# This script is used to find the combination of parameters that can help extrapolate the microburst PSD.

rb_id = 'a'
magEphemDir = '/home/mike/research/rbsp/magephem/rbsp{}'.format(rb_id)
magEphemName = 'rbsp{}_def_MagEphem_TS04D_20170331_v1.0.0.txt'.format(rb_id)

magData = spacepy.datamodel.readJSONheadedASCII(os.path.join(
    magEphemDir, magEphemName))

# Convert times
magData['DateTime'] = np.array(list(map(
        lambda x: dateutil.parser.parse(x).replace(tzinfo=None),
        magData['DateTime'])))


validL = np.where(magData['Lstar'][:, 0] != -1E31)[0]

MLTbound = [18, 20]
Lbound = [5.3, 6]

validLowMLT = (magData['EDMAG_MLT'] > MLTbound[0])
validHighMLT = (magData['EDMAG_MLT'] < MLTbound[1])
validLowL = (magData['Lstar'][:, 0] > Lbound[0])
validHighL = (magData['Lstar'][:, 0] < Lbound[1])

validInd = np.where(validLowMLT & validHighMLT & validLowL & validHighL)[0]

fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(2,1)
bPlt = fig.add_subplot(gs[0, 0], facecolor='w')
posPlt = fig.add_subplot(gs[1, 0], facecolor='w', sharex = bPlt)

bPlt.plot(magData['DateTime'], magData['BoverBeq'], 'r', label = 'B/Beq')
bPlt.scatter(magData['DateTime'][validInd], magData['BoverBeq'][validInd], c = 'b', label = 'B/Beq')
posPlt.plot(magData['DateTime'], magData['EDMAG_MLT'], label = 'MLT')
posPlt.plot(magData['DateTime'][validL], magData['Lstar'][validL, 0], label = 'L')
plt.legend()
plt.show()

### PLOT STUFF ###

###def sinAlpha(alpha, A, n):
###    """
###    This function will return a value from the function A*sin(alpha)^n. 
###    This is used for fitting the equatorial pitch angle distribution.
###    """
###    return A*np.sin(np.deg2rad(alpha))**n


###tBoundsDict = {'q':[datetime(2017, 3, 31, 11, 15, 0), datetime(2017, 3, 31, 11, 17, 0)], 'm1':[datetime(2017, 3, 31, 11, 17, 0), datetime(2017, 3, 31, 11, 17, 20)], 'm2':[datetime(2017, 3, 31, 11, 17, 10), datetime(2017, 3, 31, 11, 17, 20)],
###    'bigOne':[datetime(2017, 3, 31, 11, 17, 13), datetime(2017, 3, 31, 11, 17, 18)],
###    'smallOne':[datetime(2017, 3, 31, 11, 17, 9, 500000), datetime(2017, 3, 31, 11, 17, 10, 500000)], 'None':None}
###    tBounds = tBoundsDict['None']
