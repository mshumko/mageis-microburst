# This script will plot the magEIS time series during the interesting event on
# March 31st, 2017 where magEIS saw the loss cone filling with an accopaning 
# wave.

import sys
import matplotlib.pylab as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import dates
from datetime import datetime, timedelta
import numpy as np

sys.path.insert(0, '/home/mike/Dropbox/0_grad_work/mission_tools')
import plot_mageis_spectra

def loadEmfisisMag(path, tBounds = None):
    import spacepy.pycdf
    magData = spacepy.pycdf.CDF(path)
    
    if tBounds is not None: # Filter by times if provided.
        idt = np.where((magData['Epoch'][:] >= tBounds[0]) & 
            (magData['Epoch'][:] <= tBounds[1]))[0]
        assert len(idt) > 0, ('ERROR: no filtered spectrometer data found in '
                'the time range specified! Check tBounds keyword')

        outData = magData.copy()
        # Now filter the data
        outData['Epoch'] = outData['Epoch'][idt]
        outData['Magnitude'] = outData['Magnitude'][idt]
        outData['Mag'] = outData['Mag'][idt, :]         
    return outData

rb_id = 'A'
da = 10 # Smoothing parameter
Elow = np.array([29,46,68,95,126,164,206])
Ehigh = np.array([41,66,92,126,164,204,247])
vmin = 10**2
vmax = 10**4 #10**5/2 # During microburst time.

energies = np.arange(2)

tKey = 'quiet2'
date = datetime(2017, 3, 31)
times = {'muBurst':[datetime(2017, 3, 31, 11, 15), 
                    datetime(2017, 3, 31, 11, 25)],
           'quiet0':[datetime(2017, 3, 31, 2), 
                    datetime(2017, 3, 31, 2, 56)],
           'quiet1':[datetime(2017, 3, 31, 11),
                    datetime(2017, 3, 31, 11, 17)],
           'quiet2':[datetime(2017, 3, 31, 19, 45), 
                    datetime(2017, 3, 31, 20, 54)]            
        }
#tBounds = [datetime(2017, 3, 31, 11, 17, 1), datetime(2017, 3, 31, 11, 17, 30)]
fluxObj = plot_mageis_spectra.magEISspectra(rb_id, date, dataLevel = 3)
fluxObj.tBounds = times[tKey]
fluxObj.loadMagEIS(instrument = 'LOW', highrate = True)
fluxObj.loadMagEphem()

# Load EMFISIS magnetometer data
###path = ('/home/ms30715/ssd_data/rbsp/data/emfisis/magnetometer/L3/rbspa/' + 
###        'rbsp-a_magnetometer_1sec-gsm_emfisis-L3_20170331_v1.6.1.cdf')
###bField = loadEmfisisMag(path, tBounds = tBounds)

# Get Plot stuffs
fig = plt.figure(figsize=(10, 10), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(len(energies)+1, 15)

ax = (len(energies)+1)*[None]
ax[0] = fig.add_subplot(gs[0, :-1], facecolor='black')
ax[1:] = [fig.add_subplot(gs[i, :-1], facecolor='black', sharex = ax[0]) for i in range(1, len(ax))]
cbar_ax = fig.add_subplot(gs[:-2, -1])

maxFlux = 0

for E in energies:
    fluxObj.plotHighRateSpectra(ax = ax[E], E_ch = E, scatterS = 20, 
        pltTitle = False, plotCb = False, vmin = vmin, vmax = vmax)
    
    flatTime = fluxObj.times.flatten()
    flatFlux = fluxObj.magEISdata['HighRate'][:, :, E].flatten()
    validIdf = np.where(flatFlux != -1E31)[0]

    smoothedFlux = np.convolve(np.ones(da)/da,
                    flatFlux, mode = 'same')
    
    # Plot the local loss cone.
    ax[E].axhline(y = 180 - np.mean(fluxObj.magEphem['Loss_Cone_Alpha_s']), 
        c = 'w', ls = '--', lw = 1)
    ax[E].axhline(y = 180 + np.mean(fluxObj.magEphem['Loss_Cone_Alpha_s']), 
        c = 'w', ls = '--', lw = 1)
    ax[E].axhline(y = 360 - np.mean(fluxObj.magEphem['Loss_Cone_Alpha_n']), 
        c = 'w', ls = '--', lw = 1)
    ax[E].axhline(y = 0 + np.mean(fluxObj.magEphem['Loss_Cone_Alpha_n']), 
        c = 'w', ls = '--', lw = 1)
    
    ax[-1].plot(flatTime[validIdf], smoothedFlux[validIdf], 
        label = '{}-{} keV'.format(Elow[E], Ehigh[E]), lw = 1)
        
#ax[-2].plot(bField['Epoch'], bField['Magnitude'], label = '|B|')
#bLabel = ('Bx GSM', 'By GSM', 'Bz GSM')
#for i in range(3):
#    ax[-2].plot(bField['Epoch'], bField['Mag'][:, i], label = bLabel[i])
#ax[-2].legend(loc = 1)
#ax[-2].set_ylabel('EMFISIS B field (nT)')

ax[-1].set_yscale('log')
ax[-2].legend(loc = 1, fontsize = 8, bbox_to_anchor=(1.1, 1.05), fancybox=True, shadow=True)  
ax[-1].legend(loc = 1, fontsize = 8, bbox_to_anchor=(1.18, 1.05), fancybox=True, shadow=True)  
ax[0].set_title('RBSP-A magEIS HighRate data from {}'.format(
date.date().isoformat()))
for i in ax[:-1]:
    plt.setp(i.get_xticklabels(), visible=False)
    i.set_xlabel('')
ax[-1].set(xlabel = 'UTC', ylabel = 'counts/s')
    
for i in energies:
    ax[i].set_ylabel('pitch angle \n ({} < E < {} keV)'.format(Elow[i], Ehigh[i]))

#ax[-1].set_xlim(datetime(2017, 3, 31, 11, 17, 1), datetime(2017, 3, 31, 11, 17, 30))
# Draw the colorbar
plt.colorbar(fluxObj.sc, cax=cbar_ax, label = 'counts/s')

#secondsLoc = dates.SecondLocator()
#ax[-1].xaxis.set_minor_locator(secondsLoc)
ax[-1].set_ylim(bottom = 10**3)

gs.tight_layout(fig)
plt.show()

#print(fluxObj.times.flatten())
#print(fluxObj.magEISdata['HighRate'][:, :, E_ch].flatten())
