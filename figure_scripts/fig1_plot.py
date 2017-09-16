"""
This script and supporting functions are made to produce figure 1 in the
magEIS microburst paper

Mykhaylo Shumko
Last modified: 2017-08-10
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import plot_daves_emfisis_burst_data
import matplotlib.dates as mdates
#import matplotlib.colors
from datetime import datetime, timedelta
import matplotlib.gridspec as gridspec
#import operator

sys.path.insert(0, '/home/mike/Dropbox/0_grad_work/mission_tools')
import plot_emfisis_spectrogram
import plot_mageis_spectra
   
def get_mageis_params(ids):
    """
    This function returns a dictionary of energy and geometric conversion arrays
    that work with the level 3 magEIS data.
    """
    assert ids.upper() in ['A', 'B'], 'ERROR: Incorrect RBSP id!'
        
    # Define magEIS detector constants
    if ids.upper() == 'A':
        Emid = [34, 54, 78, 108, 143, 182, 223] # keV
        Elow = [29, 46, 68, 95, 126, 164, 206] # keV
        Ehigh = [41, 66, 92, 126, 164, 204, 247] # keV
        # Units of (keV cm^2 sr)
        G0dE = [4.13E-2, 5.73E-2, 6.056E-2, 6.88E-2, 7.35E-2, 
            6.90E-2, 5.98E-2]
        
    if ids.upper() == 'B':
        Emid = [32, 51, 74, 101, 132, 168, 208] # keV
        Elow = [27, 43, 63, 88, 117, 152, 193] # keV
        Ehigh = [39, 63, 88, 117, 150, 188] # keV
        # Units of (keV cm^2 sr)
        G0dE = [4.33E-2, 5.41E-2, 5.926E-2, 6.605E-2, 6.460E-2,
            6.23E-2, 5.96E-2]
    return {'Emid':Emid, 'Elow':Elow, 'Ehigh':Ehigh, 'G0dE':G0dE}
                
def plot_mageis(rb_id, plotType, tRange, highrate, ax = None, Alpha = None, 
        channels = None, Nsmooth = None, cax = None, downSampleAlpha = 1):
    fluxObj = plot_mageis_spectra.magEISspectra(rb_id, tRange[0], dataLevel = 3)
    fluxObj.tBounds = tRange
    fluxObj.loadMagEIS(instrument = 'LOW', highrate = highrate)
    
    if plotType == 't':
        plot_mageis_timseries(rb_id, fluxObj, ax = ax, Nsmooth = Nsmooth)
    elif plotType == 'a':
        plot_mageis_alpha(rb_id, fluxObj, channels, ax, downSampleAlpha = downSampleAlpha, cax = cax)
        
    return
    
def plot_mageis_timseries(rb_id, fluxObj, channels = None, ax = None, Nsmooth = None):
    
    ### GET magEIS PARAMS ###
    fluxObj.resolveSpinTimes(flattenTime = True)
    mageis_params = get_mageis_params(rb_id)
    
    if channels is None:
        channels = range(fluxObj.magEISdata[fluxObj.fluxKey].shape[2]-3)
        
    if ax is not None:
        for ch in channels:
            flatCountsSec = fluxObj.magEISdata[fluxObj.fluxKey][:, :, ch].flatten()
            # Smooth data with n points
            if Nsmooth is not None:
                kernel = np.ones(Nsmooth)/Nsmooth
                flatCountsSec = np.convolve(kernel, flatCountsSec, mode = 'same')
                
            validIdf = np.where(flatCountsSec != -1E31)[0]
            # If the entire interval is errors.
            if len(validIdf) != len(flatCountsSec): 
                continue
                
            ax.plot(fluxObj.times[validIdf], 
                flatCountsSec[validIdf]/mageis_params['G0dE'][ch], 
                label = '{}-{} keV'.format(mageis_params['Elow'][ch], 
                mageis_params['Ehigh'][ch]))

    # Label y-axis
    ax.set_ylabel('VAP-{} Electron flux \n'.format(rb_id.upper()) + \
        r'$(keV \ cm^2 \ sr \ s)^{-1}$')
    ax.set_yscale('log')
    return ax
    
def plot_mageis_alpha(rb_id, fluxObj, channel, ax, downSampleAlpha = 1, cax = None):
    fluxObj.plotHighRateSpectra(E_ch = channel, ax = ax, downsampleAlpha = downSampleAlpha, cax = cax, pltTitle = False,
        vmin = None, vmax = None, scatterS = 40) #cax = cbar_ax[E+1]
    mageis_params = get_mageis_params(rb_id)
    ax.set_ylabel('VAP-A MagEIS {}-{} keV\nlocal pitch angle (deg)'.format(mageis_params['Elow'][channel], 
                mageis_params['Ehigh'][channel]))

    return ax
    
def plot_emfisis(rb_id, date, tBounds, ax, cax, vmin = 10**-10, vmax = 10**-2, 
        burst_plot = False, lowF = 50):
    # Increase time range to get data for the entire desired time.
    tBounds = [tBounds[0] - timedelta(seconds = 5), tBounds[1] + timedelta(seconds = 5)]
    pObj = plot_emfisis_spectrogram.EMFISISspectra(rb_id, date, tBounds = tBounds)
    pObj.loadWFRSpectra()
    pObj.loadMagEphem()
    pObj.plotSpectra(ax = ax, spectraMax = vmax, spectraMin = vmin, lowF = lowF, 
        grid = False, plotCb = 'vertical', cAspect = 100, cax = cax, 
        printTitle = False)

    if burst_plot:
         plot_daves_emfisis_burst_data.plotDavesBurstSpectra(tBounds, vmin = vmin, 
            vmax = vmax, ax = ax, lowF = lowF)
    return ax, pObj.magEphem['EDMAG_MLT'], pObj.magEphem['Lstar'], pObj.magEphem['EDMAG_MLAT']
    
    
if __name__ == '__main__':
    ### TOP LEVEL USER INPUT ###
    zoomedT = True
    plotBurst = True
    savePlot = True
    tKey = 'muBurst'
    times = {'muBurst':[datetime(2017, 3, 31, 11, 15, 0), datetime(2017, 3, 31, 11, 18, 10)],
            'later':[datetime(2017, 3, 31, 11, 35, 0), datetime(2017, 3, 31, 11, 38)]}

    # Figure 1 time range
    tRange = times[tKey]
    plt.rcParams.update({'font.size': 10})

    # Panels you would like to plot. The value in the key:value pairs is the subplot 
    # position to plot. If you dont want to plot a panel, use a negative number.
    panelDict = {'vapa_mageis_timeseries':0, 'vapa_mageis_alpha':1,
        'vapa_emfisis_wfr':2, 'vapa_emfisis_wfr_burst':-1,
        'vapb_mageis_timeseries':3, 'vapb_mageis_alpha':4,
        'vapb_emfisis_wfr':5, 'vapa_emfisis_wfr_burst':-1}   
    # Figure out how many panels to plot.
    vals = np.array(list(panelDict.values()))
    # See how many positive numbers the panel Dict has
    Npanels = len(np.where(np.abs(vals) == vals)[0]) 

    ### Create plotting environment ###
    fig = plt.figure(figsize = (8.5, 11), dpi = 80)
    gs = gridspec.GridSpec(Npanels, 20)
    ax = Npanels*[None]
    ax[0] = fig.add_subplot(gs[0,:-1])
    for j in range(1, Npanels):
        ax[j] = fig.add_subplot(gs[j, :-1], sharex = ax[0])
    
   # Create colorbar subplots (so the time axis will line up).
    if bool(panelDict['vapa_mageis_alpha']):
        cax_mageis_vapa = fig.add_subplot(gs[panelDict['vapa_mageis_alpha'], -1])
    if bool(panelDict['vapb_mageis_alpha']):
        cax_mageis_vapb = fig.add_subplot(gs[panelDict['vapb_mageis_alpha'], -1])
    if bool(panelDict['vapa_emfisis_wfr']+1) == True:
        cax_emfisis_vapa = fig.add_subplot(gs[panelDict['vapa_emfisis_wfr'], -1])
    if bool(panelDict['vapb_emfisis_wfr']+1) == True:
        cax_emfisis_vapb = fig.add_subplot(gs[panelDict['vapb_emfisis_wfr'], -1])
        
    ### PLOT MAGEIS DATA ###
    if bool(panelDict['vapa_mageis_timeseries']+1) == True:
        axi = ax[panelDict['vapa_mageis_timeseries']]
        plot_mageis('A', 't', tRange, True, ax = axi, Nsmooth = 1)
        axi.legend(loc = 1, bbox_to_anchor=(1.2, 1))
   
    if bool(panelDict['vapb_mageis_timeseries']+1) == True:  
        axi = ax[panelDict['vapb_mageis_timeseries']]
        plot_mageis('B', 't', tRange, False, ax = axi)
        axi.legend(loc = 1, bbox_to_anchor=(1.2, 1))
        
    if bool(panelDict['vapa_mageis_alpha']+1) == True:
        axi = ax[panelDict['vapa_mageis_alpha']]
        plot_mageis('A', 'a', tRange, True, ax = axi, channels = 1, 
            downSampleAlpha = 10, cax = cax_mageis_vapa)
#        axi.set_ylabel('VAP-A MagEIS {} keV\nlocal pitch angle (deg)'.format(mageis_params['Elow'][ch], 
#                mageis_params['Ehigh'][ch]))
        axi.legend()
    if bool(panelDict['vapb_mageis_alpha']+1) == True:
        axi = ax[panelDict['vapb_mageis_alpha']]
        plot_mageis('B', 'a', tRange, False, ax = axi, channels = 1, 
            cax = cax_mageis_vapb)
#        axi.set_ylabel('VAP-B MagEIS {}\nlocal pitch angle (deg)')
        axi.legend()
    if bool(panelDict['vapa_emfisis_wfr']+1) == True:
        axi = ax[panelDict['vapa_emfisis_wfr']]
        zz, MLT_A, Lstar_A, MLAT_A = plot_emfisis('A', tRange[0], tRange, axi,
            cax_emfisis_vapa, burst_plot = plotBurst, vmax = 10**-2, vmin = 10**-10)
        # Plot the EMFISIS burst data
        axi.set_ylabel('VAP-A EMFISIS WFR \n frequency (Hz)')
    if bool(panelDict['vapb_emfisis_wfr']+1) == True:
        axi = ax[panelDict['vapb_emfisis_wfr']]
        zz, MLT_B, Lstar_B, MLAT_B = plot_emfisis('B', tRange[0], tRange, axi,
            cax_emfisis_vapb)
        axi.set_ylabel('VAP-B EMFISIS WFR \n frequency (Hz)')
        
    ax[panelDict['vapa_mageis_timeseries']].set_ylim(bottom = 10**4)
    ax[panelDict['vapb_mageis_timeseries']].set_ylim(bottom = 10**5)
    
    # Annotate the panels with position information
    rbspaText = 'L* = {}, MLT = {}, MLAT = {}'.format(round(np.mean(Lstar_A), 1), 
        round(np.mean(MLT_A), 1), round(np.mean(MLAT_A), 1))
    rbspbText = 'L* = {}, MLT = {}, MLAT = {}'.format(round(np.mean(Lstar_B), 1), 
        round(np.mean(MLT_B), 1), round(np.mean(MLAT_B), 1))
    ax[0].text(.01, .9, rbspaText, transform=ax[0].transAxes)
    ax[3].text(.01, .9, rbspbText, transform=ax[3].transAxes)
        
    # Turn off x-axis labels for all but last subplot
    for a in ax[:-1]:
        a.set_xlabel('')
    ax[-1].set_xlabel('UTC')
        
    # Draw vertical line indicating the northward and 
    # southward traveling microbursts.
    for a in ax:
        if not zoomedT:
            a.axvline(datetime(2017, 3, 31, 11, 17, 16, 650000), c = 'k')
            a.axvline(datetime(2017, 3, 31, 11, 17, 9, 800000), c = 'k')
        
        # Format the time axis
        if zoomedT:
            second5s = mdates.SecondLocator(bysecond = range(0, 60, 5), 
                interval = 1)
            a.xaxis.set_major_locator(second5s)
            a.xaxis.set_minor_locator(mdates.SecondLocator())
            a.xaxis.set_tick_params(which = 'both', width = 2, length = 5)
                
    # Label the panels
    ax[0].set_title('MagEIS and EMFISIS data from March 31st, '
        '2017 microburst event')

    if tKey == 'muBurst':
        if zoomedT:
            ax[0].set_xlim(datetime(2017, 3, 31, 11, 17, 8), 
                datetime(2017, 3, 31, 11, 17, 20))
        else:
            ax[0].set_xlim(datetime(2017, 3, 31, 11, 16, 0), datetime(2017, 3, 31, 11, 18, 10))
    
    #gs.tight_layout(fig)
    plt.subplots_adjust(left=None, bottom=0.04, right=None, top=0.97,
                wspace=0.1, hspace=0.06)
                
    if savePlot:
        if zoomedT:
            plt.savefig('fig2.pdf')
            plt.savefig('fig2.png', dpi = 200)
        else:
            plt.savefig('fig1.pdf')
            plt.savefig('fig1.png', dpi = 200)
    plt.show()

