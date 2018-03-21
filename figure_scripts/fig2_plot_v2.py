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
import matplotlib.ticker as ticker
import matplotlib.colors
from datetime import datetime, timedelta

import spacepy.pycdf

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp/')
import plot_mageis
import plot_rbspice
import fig1_plot

sys.path.insert(0, '/home/mike/research/mageis-microburst/rabbit_holes')
import fit_rbspice_pa
import rbspice_alpha_fit_popt

# "Interactive time range selection."
tKey = 'muBurst'
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
MAGEIS_MIN = 2E4
MAGEIS_MAX = 5E5
RBSPICE_MIN = 5E3
RBSPICE_MAX = 1E4

# Load the RBSPICE data
# rbspiceObj = plot_rbspice.plot_rbspice(rb_id, tBounds[0], tBounds=tBounds)
# rbspiceObj.loadData()
#### Load MagEIS data
mageisObj = plot_mageis.PlotMageis(rb_id, tBounds[0], 'highrate', 
    tRange=tBounds, instrument='low')

# Set up three panels to plot MagEIS timeseries around the microburst,
# Highlight the times used for the PSD analysis, show MagEIS and RBSPICE 
# Pitch angle evolution, and EMFISIS magnetometer data.
npanels = 3
fig = plt.figure(figsize = (11, 11), dpi = 80)
gs = gridspec.GridSpec(npanels, 15)
ax = npanels*[None]
ax[0] = fig.add_subplot(gs[0,:-2])
#axx = ax[0].twinx()
for j in range(1, npanels):
    ax[j] = fig.add_subplot(gs[j, :-2], sharex = ax[0])
alphaCbar = fig.add_subplot(gs[-1, -2:])

# Plot MagEIS
t, j = mageisObj.getFluxTimeseries(smooth=10) # Get flux
for ee in range(j.shape[1]-1):
    validj = np.where(j[:, ee] > 0)[0]
    ax[0].plot(t[validj], j[validj, ee], label='{}-{} keV'.format(mageisObj.Elow[ee], mageisObj.Ehigh[ee]))
mageisObj.plotAlpha(E_ch=0, scatterS=50, ax=ax[-1], 
    plotCb=False, pltLabels=False, cmin=MAGEIS_MIN, cmax=MAGEIS_MAX, 
    downSampleAlpha=5) # Alpha
ax[0].legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0., fontsize=10)
ax[0].set(ylabel=('MagEIS-A J \n ' + 
    r'$(cm^2 \ sr \ s \ keV)^{-1}$'), yscale='log')

# Plot arrow to microbursts
# Arrow params
arrowT = np.array([datetime(2017, 3, 31, 11, 17, 9, 777000),
            datetime(2017, 3, 31, 11, 17, 10, 250000),
            datetime(2017, 3, 31, 11, 17, 10, 860000),
            datetime(2017, 3, 31, 11, 17, 11, 500000)])
arrowJ = np.array([[2E5, 1E6], [2E5, 1E6], [2E5, 1E6], [3E5, 1E6]])
arrowCounts = np.array([[4E3, 7E3], [4E3, 7E3], [4E3, 7E3], [4E3, 7E3]])

for (tt, jj, cc) in zip(arrowT, arrowJ, arrowCounts):
    # Draw arrows to show the microburst times in the MagEIS data.
    ax[0].annotate('',
            xy=(mdates.date2num(tt), jj[0]), xycoords='data',
            xytext=(mdates.date2num(tt), jj[1]), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
    # Draw arrows to show the microburst times in the RBSPICE data.
    ax[1].annotate('',
            xy=(mdates.date2num(tt), cc[0]), xycoords='data',
            xytext=(mdates.date2num(tt), cc[1]), textcoords='data',
            arrowprops=dict(arrowstyle="<-", facecolor='black',
                            connectionstyle="arc3"),
            )

# plot RBSPICE
# ax[1], p = rbspiceObj.plotTelecopeAlphaScatter(range(6, 12), ax=ax[1], 
#     cmin=cmin, cmax=cmax, telescopes=range(4), Elabel=False, plotColorbar=False)

# PLOT EBR time series
d = spacepy.pycdf.CDF('/home/mike/research/rbsp/data/rbspice/rbspa/'
                        'rbsp-a-rbspice_lev-1_EBR_20170331_v1.1.2-01')
idT = np.where((tBounds[0] < np.array(d['Epoch'][:])) & 
                    (tBounds[1] > np.array(d['Epoch'][:])))[0]

t = np.array(d['Epoch'][idT[0]:idT[-1]])
EBR = np.array(d['EBR'][:][idT[0]:idT[-1]])
for (i) in range(5): #enumerate(np.array(d['EBR'])[:, idT].T):
    ax[1].plot(t, EBR[:,i], label='Telescope {}'.format(i+1))        

ax[1].legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0., fontsize=10)

popt = rbspice_alpha_fit_popt.popt
t0 = rbspice_alpha_fit_popt.t0
# Plot EBR pitch angles.
for i in range(5):
    sc = ax[-1].scatter(t, fit_rbspice_pa.sin_fit(mdates.date2num(t)-t0, *popt[i]), 
                        c=EBR[:, i], vmin=RBSPICE_MIN, vmax=RBSPICE_MAX, cmap=plt.get_cmap('rainbow'),
                        norm=matplotlib.colors.LogNorm())

ax[1].set(ylabel='RBSPICE-A EBR\n(counts/s)', yscale='log', ylim=(4000, 20000))
ax[-1].set(facecolor='k', title='', ylabel=r'$\alpha_{L}$ (degrees)')
ax[-1].legend()

ax[-1].set_xlabel('UTC')
ax[0].set_title("RBSP-{} from {}".format(rb_id.upper(), tBounds[0].date()))


def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

### COLORBAR CODE ###
# Plot RBSPICE colorbar
plt.colorbar(sc, ax=ax[-1], cax=alphaCbar,
    format=ticker.FuncFormatter(fmt)) #r'Flux $(keV \ cm^2 \ s \ sr)^{-1}$')
alphaCbar.set_ylabel('RBSPICE EBR (counts/s)', fontsize=10)
alphaCbar.tick_params(labelsize=10) 
alphaCbar.yaxis.set_label_position("left")


# Plot the MagEIS colorbar
alphaCbar2 = alphaCbar.twinx()
alphaCbar2.set_yscale('log')
alphaCbar2.set_ylim(MAGEIS_MIN, MAGEIS_MAX)
alphaCbar2.tick_params(labelsize=10) 
alphaCbar2.yaxis.set_label_position("right")
alphaCbar2.set_ylabel('MagEIS J '+r'$(keV \ cm^2 \ s \ sr)^{-1}$', fontsize=10)

### ADJUST ticks to be at every second ###
second5s = mdates.SecondLocator(bysecond = range(0, 60, 5), 
                interval = 1)
                
abcLabels = ['(a)', '(b)', '(c)']
abcColors = ['k', 'k', 'w']

for i, a in enumerate(ax):
    #a.xaxis.set_major_locator(second5s)
    a.xaxis.set_minor_locator(mdates.SecondLocator())
    a.xaxis.set_tick_params(which='major', width=2, length=7)
    a.xaxis.set_tick_params(which='minor', width=2, length=4)
    
    # Add subplot labels
    a.text(.01, 0.95, abcLabels[i], transform=a.transAxes, va='top', 
            color=abcColors[i])     
    # Draw vertical lines to show the resonant diffusion analysis time range.
    if i in [0, 1]:
        c = 'k'
    else:
        c = 'w'
    a.axvline(datetime(2017, 3, 31, 11, 17, 2), c=c,
        ls='--')
    a.axvline(datetime(2017, 3, 31, 11, 17, 13), c=c,
        ls='--')
    
# adjust plot time range
if tKey == 'muBurst':
    ax[0].set_xlim(tBounds[0] + timedelta(seconds=17), tBounds[1])
fig.autofmt_xdate()
gs.tight_layout(fig)
plt.show()