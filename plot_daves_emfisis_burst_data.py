import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.dates
from datetime import datetime

def plotDavesBurstSpectra(tRange = None, vmin = None, vmax = None, ax = None, lowF = 10):

    # Load Dave's data
    times = np.load('../data/emfisis_burst_times.npy')
    spectra = np.load('../data/emfisis_burst_spectra.npy')
    frequencies = 1000*np.load('../data/emfisis_burst_frequencies.npy')

    if ax is None:
        fig = plt.figure(figsize=(15, 10), dpi=80, facecolor = 'w')
        gs = gridspec.GridSpec(1,1)
        ax = fig.add_subplot(gs[0, 0], facecolor='k')
    else:
        ax = ax

    # Filter data by the time range
    if tRange is not None:
        idt = np.where((times[:-1] > tRange[0]) & (times[:-1] < tRange[1]))[0]
        times = times[idt]
        spectra = spectra[idt, :]

    # Break up the data into chunks with continuous burst data.
    dt = np.nan*np.ones(len(times) - 2, dtype = float)
    for i in range(len(dt)):
        dt[i] = (times[i+1] - times[i]).total_seconds()

    tBreaks = np.where(dt > 1)[0]
    tBreaks = np.insert(tBreaks, 0, 0)
    tBreaks = np.append(tBreaks, len(times) - 1)

    for ts in range(len(tBreaks)-1):
        wTT, wFF = np.meshgrid(times[tBreaks[ts]+1:tBreaks[ts+1]], frequencies)
        #wTT = np.broadcast_to(times, spectra.shape)

        cs = ax.pcolormesh(wTT, wFF, np.transpose(spectra[tBreaks[ts]+1:tBreaks[ts+1], :]),
                     cmap = plt.get_cmap('gnuplot2'), norm=colors.LogNorm(),
                     vmin = vmin, vmax = vmax)
        ax.set_ylim(bottom = lowF)
    return cs, ax

if __name__ == '__main__':
    tRange = [datetime(2017, 3, 31, 11, 6), datetime(2017, 3, 31, 12, 30)]
    cs, ax = plotDavesBurstSpectra(tRange)   
    plt.show()        
