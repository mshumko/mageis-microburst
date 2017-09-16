import sys
sys.path.insert(0, '/home/ms30715/research_code/auxiliary')
#import plot_emfisis_spectrogram
import plot_mageis_spectra
import matplotlib.pylab as plt
from datetime import datetime
import matplotlib.gridspec as gridspec
import numpy as np

rb_id = 'A'
E_ch = 0
Elow = np.array([29,46,68,95,126,164,206])
Ehigh = np.array([41,66,92,126,164,204,247])
date = datetime(2017, 3, 31)
tBounds = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 25)]

fig = plt.figure(figsize=(10, 10), dpi=80, facecolor = 'white')
gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0, 0], facecolor='white')

# PLOT magEIS HighRate data
fluxObj = plot_mageis_spectra.magEISspectra(rb_id, date, dataLevel = 3)
fluxObj.tBounds = tBounds
fluxObj.loadMagEIS(instrument = 'LOW', highrate = True)
#print(fluxObj.magEISdata['HighRate_Alpha360'].shape)

for dt in range(len(fluxObj.magEISdata['Epoch'])):
    validInd = np.where(fluxObj.magEISdata['HighRate'][dt, :, E_ch] != -1E31)[0]
    ax.errorbar(fluxObj.magEISdata['HighRate_Alpha360'][dt, validInd], 
        fluxObj.magEISdata['HighRate'][dt, validInd, E_ch],
        yerr = np.sqrt(fluxObj.magEISdata['HighRate'][dt, validInd, E_ch]),
        c = 'k', lw = 0, marker = 'o')  
    ax.text(0.01, .95, fluxObj.magEISdata['Epoch'][dt].isoformat(), 
                    transform = ax.transAxes, fontsize = 15)
    #ax.set_ylim(0, 50000)              
    ax.set_yscale('log')
    #ax.set_ylim((1000, 50000))
    ax.set_ylim(100, np.max(fluxObj.magEISdata['HighRate'][:, :, E_ch]))                
    ax.set(title = 'RBSP-A magEIS LOW electron flux \n {} < E < {} keV'.format(Elow[E_ch], Ehigh[E_ch]), xlabel = 'Pitch angle', ylabel = 'Counts/s')
    gs.tight_layout(fig)
    #plt.show()
    plt.savefig('../plots/movies/raw/ch{}/magEIS_eflux_t_{}_ch_{}.png'.format(
        E_ch, fluxObj.magEISdata['Epoch'][dt].strftime('%Y%m%d_%H%M%S'), E_ch))   
        
    plt.cla()

#fluxObj.plotHighRateSpectra(E_ch = E_ch, ax = ax, downsampleAlpha = 10, vmin = vmin, vmax = vmax, scatterS = 40)

### PLOT magEIS cutpouts ###

# Need to run this command to make it into a gif "convert -delay 25 *.png test.gif"

