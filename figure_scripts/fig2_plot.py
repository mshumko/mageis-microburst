"""
This script and supporting functions are made to produce figure 2 in the
MagEIS microburst paper

Mykhaylo Shumko
Last modified: 2017-10-18
"""

import numpy as np
import sys
import os
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')))
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from datetime import datetime, timedelta

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp/')
import plot_mageis
import plot_rbspice
import fig1_plot

# "Interactive time range selection."
tKey = 'muBurst'
rb_id = 'B'
times = {'muBurst':[datetime(2017, 3, 31, 11, 17, 0), 
                    datetime(2017, 3, 31, 11, 17, 20)],
            'later':[datetime(2017, 3, 31, 11, 35, 0), 
                    datetime(2017, 3, 31, 11, 38)],
            'all':[datetime(2017, 3, 31, 11, 15), 
                    datetime(2017, 3, 31, 11, 20)]}
tBounds = times[tKey]
cmin = 2E4
cmax = 2E5

# Load the RBSPICE data
rbspiceObj = plot_rbspice.plot_rbspice(rb_id, tBounds[0], tBounds=tBounds)
rbspiceObj.loadData()
#### Load MagEIS data
###mageisObj = plot_mageis.magEISspectra(rb_id, tBounds[0], dataLevel = 3)
###mageisObj.tBounds = tBounds
###mageisObj.loadMagEIS(instrument = 'LOW', highrate = False)

# Set up three panels to plot MagEIS timeseries around the microburst,
# Highlight the times used for the PSD analysis, show MagEIS and RBSPICE 
# Pitch angle evolution, and EMFISIS magnetometer data.
npanels = 2
fig = plt.figure(figsize = (11, 11), dpi = 80)
gs = gridspec.GridSpec(npanels, 10)
ax = npanels*[None]
ax[0] = fig.add_subplot(gs[0,:-1])
for j in range(1, npanels):
    ax[j] = fig.add_subplot(gs[j, :-1], sharex = ax[0])
alphaCbar = fig.add_subplot(gs[1, -1])

# Plot MagEIS
###mageisObj.plotHighRateTimeSeries(ax=ax[0], smooth=10,
###    chLegend=False) # flux
###mageisObj.plotHighRateSpectra(E_ch=1, scatterS=50, ax=ax[1], 
###    plotCb=False, pltTitle=False, pltXlabel=False, cmin=cmin, cmax=cmax) # Alpha
###ax[0].legend(bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0.)
###ax[0].set(ylabel=('MagEIS-LOW electron flux \n ' + 
###    r'$(cm^2 \ sr \ s \ keV)^{-1}$'))
fig1_plot.plot_mageis(rb_id.upper(), 't', tBounds, False, ax = ax[0])
fig1_plot.plot_mageis(rb_id.upper(), 'a', tBounds, False, ax = ax[1])


# plot RBSPICE
ax[1], p = rbspiceObj.plotTelecopeAlphaScatter(range(14, 20), ax=ax[1], 
    cmin=cmin, cmax=cmax)
plt.colorbar(p, ax=ax[1], cax=alphaCbar, label=r'Flux $(keV \ cm^2 \ s \ sr)^{-1}$')
ax[1].set(facecolor='k', title='', ylabel=r'$\alpha_{sc}$')
ax[-1].set_xlabel('UTC')
plt.suptitle("RBSP-{} from {}".format(rb_id.upper(), tBounds[0].date()))

# Ticks at every second
second5s = mdates.SecondLocator(bysecond = range(0, 60, 5), 
                interval = 1)
for a in ax:
    #a.xaxis.set_major_locator(second5s)
    a.xaxis.set_minor_locator(mdates.SecondLocator())
    a.xaxis.set_tick_params(which='major', width=2, length=7)
    a.xaxis.set_tick_params(which='minor', width=2, length=4)

fig.autofmt_xdate()
gs.tight_layout(fig)
plt.show()
    


