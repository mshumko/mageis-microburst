"""
This script and supporting functions are made to produce figure 1 in the
magEIS microburst paper

Mykhaylo Shumko
Last modified: 2017-08-10
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')))
import matplotlib.pyplot as plt
import plot_daves_emfisis_burst_data
import matplotlib.dates as mdates
#import matplotlib.colors
from datetime import datetime, timedelta
import matplotlib.gridspec as gridspec
#import operator

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp/')
import plot_mageis as plot_mageis_lib
import plot_emfisis_spectra
   
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
        channels = None, Nsmooth = None, cax = None, downSampleAlpha = 1, vmin=None, vmax=None):
    fluxObj = plot_mageis_lib.PlotMageis(rb_id, tRange[0], 'highrate', tRange=tRange, instrument='low')
    
    if plotType == 't':
        ax = fluxObj.plotTimeseries(ax=ax, pltLabels=False)
        ax.set_ylabel('MagEIS LOW electron flux \n'.format(rb_id.upper()) + \
        r'$(keV \ cm^2 \ sr \ s)^{-1}$')
    elif plotType == 'a':
        ax = fluxObj.plotAlpha(ax=ax, cax=cax, cmin=vmin, cmax=vmax, E_ch=channels, pltLabels=False, scatterS=40)
        ax.set_ylabel('MagEIS {}-{} keV\n'.format(
            fluxObj.Elow[channels], fluxObj.Ehigh[channels]) +
            r'$\alpha_{L}$ (deg)')  
    return ax
    
def plot_emfisis(rb_id, date, tBounds, ax, cax, vmin = 10**-10, vmax = 10**-2, 
        burst_plot = False, lowF = 50):
    # Increase time range to get data for the entire desired time.
    tBounds = [tBounds[0] - timedelta(seconds = 5), tBounds[1] + timedelta(seconds = 5)]
    pObj = plot_emfisis_spectra.EMFISISspectra(rb_id, date, tBounds = tBounds)
    pObj.loadWFRSpectra()
    pObj.loadMagEphem(Bmodel='TS04D')
    pObj.plotSpectra(ax = ax, spectraMax = vmax, spectraMin = vmin, lowF = lowF, 
        grid = False, plotCb = 'vertical', cAspect = 100, cax = cax, 
        printTitle = False)

    if burst_plot:
         plot_daves_emfisis_burst_data.plotDavesBurstSpectra(tBounds, vmin = vmin, 
            vmax = vmax, ax = ax, lowF = lowF)
    return ax, pObj.magEphem['EDMAG_MLT'], pObj.magEphem['Lstar'], pObj.magEphem['EDMAG_MLAT']
    
    
if __name__ == '__main__':
    ### TOP LEVEL USER INPUT ###
    zoomedT = False
    plotBurst = True
    savePlot = False
    tKey = 'muBurst'
    times = {'muBurst':[datetime(2017, 3, 31, 11, 15, 0), datetime(
                        2017, 3, 31, 11, 18, 10)],
            'later':[datetime(2017, 3, 31, 11, 35, 0), datetime(
                        2017, 3, 31, 11, 38)],
            'all':[datetime(2017, 3, 31, 11, 15), datetime(
                        2017, 3, 31, 11, 20)]}

    # Figure 1 time range
    tRange = times[tKey]
    plt.rcParams.update({'font.size': 10})

    # Panels you would like to plot. The value in the key:value pairs is the subplot 
    # position to plot. If you dont want to plot a panel, use a negative number.
    panelDict = {'rbspa_mageis_timeseries':0, 'rbspa_mageis_alpha':1,
        'rbspa_emfisis_wfr':2, 'rbspa_emfisis_wfr_burst':-1,
        'rbspb_mageis_timeseries':3, 'rbspb_mageis_alpha':4,
        'rbspb_emfisis_wfr':5, 'rbspa_emfisis_wfr_burst':-1}   
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
    if bool(panelDict['rbspa_mageis_alpha']):
        cax_mageis_rbspa = fig.add_subplot(gs[panelDict['rbspa_mageis_alpha'], -1])
    if bool(panelDict['rbspb_mageis_alpha']):
        cax_mageis_rbspb = fig.add_subplot(gs[panelDict['rbspb_mageis_alpha'], -1])
    if bool(panelDict['rbspa_emfisis_wfr']+1) == True:
        cax_emfisis_rbspa = fig.add_subplot(gs[panelDict['rbspa_emfisis_wfr'], -1])
    if bool(panelDict['rbspb_emfisis_wfr']+1) == True:
        cax_emfisis_rbspb = fig.add_subplot(gs[panelDict['rbspb_emfisis_wfr'], -1])
        
    ### PLOT MAGEIS DATA ###
    if bool(panelDict['rbspa_mageis_timeseries']+1) == True:
        axi = ax[panelDict['rbspa_mageis_timeseries']]
        plot_mageis('A', 't', tRange, True, ax = axi, Nsmooth = 1)
        axi.legend(loc = 1, bbox_to_anchor=(1.2, 1))
   
    if bool(panelDict['rbspb_mageis_timeseries']+1) == True:  
        axi = ax[panelDict['rbspb_mageis_timeseries']]
        plot_mageis('B', 't', tRange, False, ax = axi)
        axi.legend(loc = 1, bbox_to_anchor=(1.2, 1))
        
    if bool(panelDict['rbspa_mageis_alpha']+1) == True:
        axi = ax[panelDict['rbspa_mageis_alpha']]
        plot_mageis('A', 'a', tRange, True, ax = axi, channels = 1, 
            downSampleAlpha=10, cax = cax_mageis_rbspa, vmin=2E4, vmax=2E5)
#        axi.set_ylabel('rbsp-A MagEIS {} keV\nlocal pitch angle (deg)'.format(mageis_params['Elow'][ch], 
#                mageis_params['Ehigh'][ch]))
        axi.legend()
    if bool(panelDict['rbspb_mageis_alpha']+1) == True:
        axi = ax[panelDict['rbspb_mageis_alpha']]
        plot_mageis('B', 'a', tRange, False, ax = axi, channels = 1, 
            cax = cax_mageis_rbspb)
#        axi.set_ylabel('rbsp-B MagEIS {}\nlocal pitch angle (deg)')
        axi.legend()
        
    ### PLOT EMFISIS DATA ###
    if tKey == 'muBurst': # Adjust the time range to get the full f_ce curves.
        emfisisTRange = [tRange[0], tRange[1] + timedelta(minutes=1)]
    if bool(panelDict['rbspa_emfisis_wfr']+1) == True:
        axi = ax[panelDict['rbspa_emfisis_wfr']]
        zz, MLT_A, Lstar_A, MLAT_A=plot_emfisis('A', tRange[0], emfisisTRange, axi,
            cax_emfisis_rbspa, burst_plot=plotBurst, vmax=10**-2, vmin=10**-10)
        # Plot the EMFISIS burst data
        axi.set_ylabel('EMFISIS WFR \n frequency (Hz)')
        axi.legend_.remove() # Remove legend
    if bool(panelDict['rbspb_emfisis_wfr']+1) == True:
        axi = ax[panelDict['rbspb_emfisis_wfr']]
        zz, MLT_B, Lstar_B, MLAT_B=plot_emfisis('B', tRange[0], emfisisTRange, axi,
            cax_emfisis_rbspb)
        axi.set_ylabel('EMFISIS WFR \n frequency (Hz)')
        axi.legend_.remove() # Remove legen
        
    ax[panelDict['rbspa_mageis_timeseries']].set_ylim(bottom = 10**4)
    ax[panelDict['rbspb_mageis_timeseries']].set_ylim(bottom = 10**5)
    
    # Annotate the panels with position information
    rbspaText = 'L* = {}\nMLT = {}\nMLAT = {}'.format(round(np.mean(Lstar_A), 1), 
        round(np.mean(MLT_A), 1), round(np.mean(MLAT_A), 1))
    rbspbText = 'L* = {}\nMLT = {}\nMLAT = {}'.format(round(np.mean(Lstar_B), 1), 
        round(np.mean(MLT_B), 1), round(np.mean(MLAT_B), 1))
    ax[0].text(.01, 0.85, rbspaText, transform=ax[0].transAxes, va='top')
    ax[3].text(.01, 0.85, rbspbText, transform=ax[3].transAxes, va='top')
        
    # Turn off x-axis labels for all but last subplot
    for a in ax[:-1]:
        a.set_xlabel('')
    ax[-1].set_xlabel('UTC')
        
    # Draw vertical line indicating the time interval used
    # in the resonant diffusion analysis, and quiet period.
    for a in ax[0:2]:
        a.axvline(datetime(2017, 3, 31, 11, 15, 0), c='k')
        a.axvline(datetime(2017, 3, 31, 11, 16,50), c='k')
        a.axvline(datetime(2017, 3, 31, 11, 17, 2), c='k',
            ls='--')
        a.axvline(datetime(2017, 3, 31, 11, 17, 13), c='k',
            ls='--')        
    # Format the time axis
    #if zoomedT:
    abcLabels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
    abcColors = ['k', 'k', 'w', 'k', 'k', 'w']
    for i, a in enumerate(ax):
        #second5s = mdates.SecondLocator(bysecond = range(0, 60, 5), 
        #    interval = 1)
        #a.xaxis.set_major_locator(second5s)
        #a.xaxis.set_minor_locator(mdates.SecondLocator())
        a.xaxis.set_tick_params(which = 'both', width = 2, length = 5)
        # Add panel labels
        a.text(.01, 0.95, abcLabels[i], transform=a.transAxes, va='top', 
            color=abcColors[i])        
    # Label the panels
    ax[0].set_title('MagEIS and EMFISIS data from March 31st, '
        '2017 microburst event')

    if tKey == 'muBurst':
        if zoomedT:
            ax[0].set_xlim(datetime(2017, 3, 31, 11, 17, 1), 
                datetime(2017, 3, 31, 11, 17, 12))
#            ax[0].set_xlim(datetime(2017, 3, 31, 11, 17, 8), 
#                datetime(2017, 3, 31, 11, 17, 20))
        else:
            ax[0].set_xlim(datetime(2017, 3, 31, 11, 16, 0), datetime(2017, 3, 31, 11, 18, 10))
    
    #gs.tight_layout(fig)
    plt.subplots_adjust(left=0.15, bottom=0.04, right=None, top=0.97,
                wspace=0.1, hspace=0.06)
                
    # Add RBSP labels
    labelAx = plt.axes([0, 0, .05, 1], frameon=False)
    labelAx.text(0.75, 0.75, 'RBSP-A', va='center', ha='right',
        rotation='vertical', fontsize=15)
    labelAx.annotate("",
            xy=(0.75, 0.5), xycoords='data',
            xytext=(0.75, 1), textcoords='data',
            arrowprops=dict(arrowstyle="<->",
                            connectionstyle="arc3"),
            )
    labelAx.text(0.75, 0.25, 'RBSP-B', va='center', ha='right',
        rotation='vertical', fontsize=15)
    labelAx.annotate("",
            xy=(0.75, 0), xycoords='data',
            xytext=(0.75, 0.5), textcoords='data',
            arrowprops=dict(arrowstyle="<->",
                            connectionstyle="arc3"),
            )
    
    for a in ax[:-1]:
        plt.setp(a.get_xticklabels(), visible=False)
###    if savePlot:
######        if zoomedT:
######            plt.savefig('fig2.pdf')
######            plt.savefig('fig2.png', dpi = 200)
######        else:
###        plt.savefig('fig1.pdf')
###        plt.savefig('fig1.png', dpi = 200)
    plt.show()

