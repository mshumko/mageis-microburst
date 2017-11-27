# Plot magEIS counts timeseries

import sys
sys.path.insert(0, '/home/ms30715/research_code/auxiliary')
import plot_emfisis_spectrogram
import plot_mageis_spectra
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import matplotlib.gridspec as gridspec
import numpy as np

rb_id = 'A'
date = datetime(2017, 3, 31)
#tBounds = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 25)]
tBounds = [datetime(2017, 3, 31, 11, 17), datetime(2017, 3, 31, 11, 17, 15)]

# magEIS energy boundaries
Elow = np.array([29,46,68,95,126,164,206])
Ehigh = np.array([41,66,92,126,164,204,247])
fluxE = np.arange(4)

fluxObj = plot_mageis_spectra.magEISspectra(rb_id, date, dataLevel = 3)
fluxObj.tBounds = tBounds
fluxObj.loadMagEIS(instrument = 'LOW', highrate = True)

### PLOT STUFFS ###
fig = plt.figure(figsize=(10, 5), dpi = 80, facecolor = 'white')
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0, 0], facecolor='k')

#fig2 = plt.figure(figsize=(5, 5), dpi = 200, facecolor = 'white')
#gsx = gridspec.GridSpec(1, 1)
#bx = fig.add_subplot(gsx[0, 0], facecolor='k')

fluxObj.resolveSpinTimes(flattenTime = True)
for E in range(len(fluxE)):
    flatCountsSec = fluxObj.magEISdata['HighRate'][:, :, fluxE[E]].flatten()
    flatCountsBin = (10.9/1000)*flatCountsSec # 10.9 s spin period, 1000 sectors per spin
    validIdf = np.where(flatCountsSec != -1E31)[0]

    ax.errorbar(fluxObj.times[validIdf], flatCountsBin[validIdf], 
        yerr = np.sqrt(flatCountsBin[validIdf]),
        label = '{} < E < {} keV'.format(Elow[fluxE[E]], Ehigh[fluxE[E]]), lw = 1)
        
        
ax.set(title = 'magEIS electron counts during microburst', xlabel = 'UTC', 
    ylabel = 'counts/sector')
ax.set_ylim(bottom = 10)
ax.set_yscale('log')
ax[-1].xaxis.set_minor_locator(mdates.SecondLocator())
ax.legend()
plt.show()


