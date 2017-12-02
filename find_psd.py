# This script is used to find other times to help extrapolate the microburst 
# PSD. The MLTbound and Lbound variables define the combination of parameters
# used to filter out times to investigate.

import numpy as np
import os
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime
import dateutil.parser

import spacepy.datamodel

rb_id = 'a'
MLTbound = [18, 20]
Lbound = [5.3, 6]
date = 20170331

# Load magnetic ephemeris.
magEphemDir = '/home/mike/research/rbsp/magephem/rbsp{}'.format(rb_id)
magEphemName = 'rbsp{}_def_MagEphem_TS04D_{}_v1.0.0.txt'.format(rb_id, date)
magData = spacepy.datamodel.readJSONheadedASCII(os.path.join(
    magEphemDir, magEphemName))
# Convert times
magData['DateTime'] = np.array(list(map(
        lambda x: dateutil.parser.parse(x).replace(tzinfo=None),
        magData['DateTime'])))
validL = np.where(magData['Lstar'][:, 0] != -1E31)[0]

# FIlter the magnetic ephemeris data by L and MLT.
validLowMLT = (magData['EDMAG_MLT'] > MLTbound[0])
validHighMLT = (magData['EDMAG_MLT'] < MLTbound[1])
validLowL = (magData['Lstar'][:, 0] > Lbound[0])
validHighL = (magData['Lstar'][:, 0] < Lbound[1])
validInd = np.where(validLowMLT & validHighMLT & validLowL & validHighL)[0]

# Visualize L, MLT, and B_l/B_sc to determine how close RBSP is to the magnetic
# equator.
fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(2,1)
bPlt = fig.add_subplot(gs[0, 0], facecolor='w')
posPlt = fig.add_subplot(gs[1, 0], facecolor='w', sharex = bPlt)

bPlt.plot(magData['DateTime'], magData['BoverBeq'], 'r', label='B/Beq')
bPlt.scatter(magData['DateTime'][validInd], magData['BoverBeq'][validInd], 
    c='b', label='B/Beq (L & MLT satisfied)')
posPlt.plot(magData['DateTime'], magData['EDMAG_MLT'], label='MLT')
posPlt.plot(magData['DateTime'][validL], magData['Lstar'][validL, 0], 
    label='L')
bPlt.legend()
posPlt.legend()

bPlt.set(title='Investigating RBSP-{} MagEphem position on {}\n{} < L < {}  |  {} < MLT < {}'.format(
    rb_id.upper(), date, *Lbound, *MLTbound), ylabel=r'$B_L/B_eq$')
posPlt.set(xlabel='UTC')

fig.autofmt_xdate()
plt.show()