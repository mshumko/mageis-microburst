import sys
sys.path.insert(0, '/home/mike/Dropbox/0_grad_work/mission_tools')
import plot_emfisis_spectrogram
import plot_mageis_spectra
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import matplotlib.gridspec as gridspec
import numpy as np

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

### PLOTTING PARAMETERS ###
rb_id = 'A'
date = datetime(2017, 3, 31)
tBounds = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 25)]
#tBounds = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 20, 10)]

vmin = 100
vmax = 10**5

# magEIS energy boundaries
Elow = np.array([29,46,68,95,126,164,206])
Ehigh = np.array([41,66,92,126,164,204,247])
fluxE = np.arange(2)
#fluxE = [0, 2]

### Load EMFISIS magnetometer data ###
path = ('/home/ms30715/ssd_data/rbsp/data/emfisis/magnetometer/L3/rbspa/' + 
        'rbsp-a_magnetometer_1sec-gsm_emfisis-L3_20170331_v1.6.1.cdf')
bField = loadEmfisisMag(path, tBounds = tBounds)

fig = plt.figure(figsize=(10, 10), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec((len(fluxE)+2), 15)
ax = (len(fluxE)+2)*[None]
ax[0] = fig.add_subplot(gs[0, :-1], facecolor='k')
for i in range(1, len(fluxE)+2):
    ax[i] = fig.add_subplot(gs[i, :-1], facecolor='k', sharex = ax[0])
    
cbar_ax = (len(fluxE)+1)*[None]
for i in range(len(cbar_ax)):
    cbar_ax[i] = fig.add_subplot(gs[i, -1])

# PLOT magEIS HighRate data
fluxObj = plot_mageis_spectra.magEISspectra(rb_id, date, dataLevel = 3)
fluxObj.tBounds = tBounds
fluxObj.loadMagEIS(instrument = 'LOW', highrate = True)
fluxObj.loadMagEphem()
#for i in fluxE:
#    fluxObj.plotHighRateSpectra(E_ch = i, ax = ax[i+1], downsampleAlpha = 10, vmin = vmin, vmax = None, scatterS = 40)
#    
    
### Plot EMFISIS data spectral data ###
pObj = plot_emfisis_spectrogram.EMFISISspectra(rb_id, date, tBounds = tBounds)
pObj.loadWFRSpectra()
pObj.loadMagEphem()
pObj.plotWFRSpectra(ax = ax[0], spectraMax = 10**-2, lowF = 10, grid = False,
    plotCb = 'vertical', cAspect = 10, cax = cbar_ax[0])
    
### PLOT TIME SERIES ###
for E in range(len(fluxE)):
    fluxObj.plotHighRateSpectra(E_ch = fluxE[E], ax = ax[E+1], downsampleAlpha = 10, 
        vmin = vmin, vmax = None, scatterS = 40, cax = cbar_ax[E+1])
    
    flatTime = fluxObj.times.flatten()
    flatCountsSec = fluxObj.magEISdata['HighRate'][:, :, fluxE[E]].flatten()
    flatCountsBin = (10.9/1000)*flatCountsSec # 10.9 s spin period, 1000 sectors per spin
    validIdf = np.where(flatCountsSec != -1E31)[0]
    
##    ax[E].axhline(y = 180, c = 'w', ls = '--') # Southern center of loss cone.
##     Plot the local loss cone.
    ax[E+1].axhline(y = 180 - np.mean(fluxObj.magEphem['Loss_Cone_Alpha_s']), 
        c = 'w', ls = '--', lw = 1)
    ax[E+1].axhline(y = 180 + np.mean(fluxObj.magEphem['Loss_Cone_Alpha_s']), 
        c = 'w', ls = '--', lw = 1)
    ax[E+1].axhline(y = 360 - np.mean(fluxObj.magEphem['Loss_Cone_Alpha_n']), 
        c = 'w', ls = '--', lw = 1)
    ax[E+1].axhline(y = 0 + np.mean(fluxObj.magEphem['Loss_Cone_Alpha_n']), 
        c = 'w', ls = '--', lw = 1)
    
#    ax[-1].plot(flatTime[validIdf], flatCountsSec[validIdf], 
#        label = '{} < E < {} keV'.format(Elow[fluxE[E]], Ehigh[fluxE[E]]), lw = 1)
    ax[-1].errorbar(flatTime[validIdf], flatCountsBin[validIdf], 
        yerr = np.sqrt(flatCountsBin[validIdf]),
    label = '{} < E < {} keV'.format(Elow[fluxE[E]], Ehigh[fluxE[E]]), lw = 1)

ax[-1].set(yscale = 'log')
#ax[-1].set_ylim(bottom = 10**3)

### DO ALL OF THE TICK MANIPULATIONS HERE ###
# fig.suptitle(
ax[0].set_title('RBSP-A EMFISIS and magEIS data from {}'.format(
date.date().isoformat()))
#ax[0].set_title('RBSP-A EMFISIS and magEIS data from {}'.format(
#date.date().isoformat()))
for i in ax[:-1]:
    plt.setp(i.get_xticklabels(), visible=False)
    i.set_xlabel('')
    
for i in ax[1:]: 
    i.set_title('')
    
ax[0].set_ylabel('EMFISIS WFR \n frequency (Hz)')
ax[-1].set_ylabel('counts/sector')
#ax[-1].legend(loc = 1)
ax[-1].legend(loc = 1, fontsize = 8, bbox_to_anchor=(1.18, 1.05), fancybox=True, shadow=True)  
ax[-1].set_xlabel('UTC')

for i in range(len(fluxE)):
    ax[i+1].set_ylabel('magEIS Pitch angle \n ({} < E < {} keV)'.format(Elow[fluxE[i]], Ehigh[fluxE[i]]))

ax[0].set_xlim([datetime(2017, 3, 31, 11, 17), datetime(2017, 3, 31, 11, 17, 12)])

# Adjust the ticks to be seconds
#ax.xaxis.set_major_locator(mdates.SecondsLocator())
#ax.xaxis.set_major_formatter(mdates.DateFormatter('%a %d'))
ax[-1].xaxis.set_minor_locator(mdates.SecondLocator())

gs.tight_layout(fig) # h_pad=0
#plt.subplots_adjust(right=0.99)
#plt.savefig('test.png', dpi = 600)
plt.show()
