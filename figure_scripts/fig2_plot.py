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
import spacepy.pycdf

# "Interactive time range selection."
tKey = 'muBurst2'
rb_id = 'A'
if rb_id == 'A':
    highrate = True
else:
    highrate = False
times = {'muBurst':[datetime(2017, 3, 31, 11, 16, 40), 
                    datetime(2017, 3, 31, 11, 17, 20)],
            'muBurst2':[datetime(2017, 3, 31, 11, 17), 
                    datetime(2017, 3, 31, 11, 18, 20)],
            'rbspb_b_peaks':[datetime(2017, 3, 31, 11, 17), 
                    datetime(2017, 3, 31, 11, 18, 30)],
            'later':[datetime(2017, 3, 31, 11, 35, 0), 
                    datetime(2017, 3, 31, 11, 38)],
            'all':[datetime(2017, 3, 31, 11, 15), 
                    datetime(2017, 3, 31, 11, 20)]}
tBounds = times[tKey]
cmin = 2E4
cmax = 5E5

# Load the RBSPICE data
rbspiceObj = plot_rbspice.plot_rbspice(rb_id, tBounds[0], tBounds=tBounds)
rbspiceObj.loadData()
#### Load MagEIS data
mageisObj = plot_mageis.PlotMageis(rb_id, tBounds[0], 'highrate', 
    tRange=tBounds, instrument='low')

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
t, j = mageisObj.getFluxTimeseries(smooth=10) # Get flux
print(j.shape)
for ee in range(j.shape[1]-1):
    validj = np.where(j[:, ee] > 0)[0]
    ax[0].plot(t[validj], j[validj, ee], label='{}-{} keV'.format(mageisObj.Elow[ee], mageisObj.Ehigh[ee]))
mageisObj.plotAlpha(E_ch=0, scatterS=50, ax=ax[1], 
    plotCb=False, pltLabels=False, cmin=cmin, cmax=cmax, 
    downSampleAlpha=5) # Alpha
ax[0].legend(bbox_to_anchor=(0.95, 1), loc=2, borderaxespad=0.)
ax[0].set(ylabel=('MagEIS-A J \n ' + 
    r'$(cm^2 \ sr \ s \ keV)^{-1}$'), yscale='log')

# Plot arrow to microbursts
# Arrow params
t = np.array([datetime(2017, 3, 31, 11, 17, 9, 777000),
            datetime(2017, 3, 31, 11, 17, 10, 250000),
            datetime(2017, 3, 31, 11, 17, 10, 860000),
            datetime(2017, 3, 31, 11, 17, 11, 500000)])
j = np.array([
    [2E5, 1E6], [2E5, 1E6], [2E5, 1E6], [3E5, 1E6]
])

for (tt, jj) in zip(t, j):
    ax[0].annotate('',
            xy=(mdates.date2num(tt), jj[0]), xycoords='data',
            xytext=(mdates.date2num(tt), jj[1]), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )

# plot RBSPICE
ax[1], p = rbspiceObj.plotTelecopeAlphaScatter(range(6, 12), ax=ax[1], 
    cmin=cmin, cmax=cmax, telescopes=range(4), Elabel=False, plotColorbar=False)
plt.colorbar(p, ax=ax[1], cax=alphaCbar, label=r'Flux $(keV \ cm^2 \ s \ sr)^{-1}$')
ax[1].set(facecolor='k', title='', ylabel='MagEIS 29-41 keV | RBSPICE 28-40 keV\n' + r'$\alpha_{L}$ (degrees)')

# plot EMFISIS magnetometer data
###b = spacepy.pycdf.CDF('/home/mike/research/rbsp/data/emfisis/rbsp{}/'
###    'rbsp-{}_magnetometer_4sec-geo_emfisis-L3_20170331_v1.6.1.cdf'.format(
###    rb_id.lower(), rb_id.lower()))
###idB = np.where((b['Epoch'][:] > tBounds[0]) & (b['Epoch'][:] < tBounds[1]+timedelta(minutes=1)))[0] 
###ax[-1].plot(np.array(b['Epoch'])[idB], np.array(b['Magnitude'])[idB])
###ax[-1].set_ylabel('|B| (nT)')

ax[-1].set_xlabel('UTC')
ax[0].set_title("RBSP-{} from {}".format(rb_id.upper(), tBounds[0].date()))

# Ticks at every second
second5s = mdates.SecondLocator(bysecond = range(0, 60, 5), 
                interval = 1)
                
abcLabels = ['(a)', '(b)']
abcColors = ['k', 'w']

for i, a in enumerate(ax):
    #a.xaxis.set_major_locator(second5s)
    a.xaxis.set_minor_locator(mdates.SecondLocator())
    a.xaxis.set_tick_params(which='major', width=2, length=7)
    a.xaxis.set_tick_params(which='minor', width=2, length=4)
    
    # Add subplot labels
    a.text(.01, 0.95, abcLabels[i], transform=a.transAxes, va='top', 
            color=abcColors[i])     
    
# Draw vertical lines to show the resonant diffusion analysis time range.
ax[0].axvline(datetime(2017, 3, 31, 11, 17, 2), c='k',
    ls='--')
ax[0].axvline(datetime(2017, 3, 31, 11, 17, 13), c='k',
    ls='--')
ax[-1].axvline(datetime(2017, 3, 31, 11, 17, 2), c='w',
    ls='--')
ax[-1].axvline(datetime(2017, 3, 31, 11, 17, 13), c='w',
    ls='--')
if tKey == 'muBurst':
    print(tBounds[0] - timedelta(seconds=6))
    ax[0].set_xlim(tBounds[0] + timedelta(seconds=17), tBounds[1])
fig.autofmt_xdate()
gs.tight_layout(fig)
plt.show()
    


