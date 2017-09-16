# Fitting functions and scripts for the phase space density of quiet times right
# before the March 31st, 2017 microburst observation

import numpy as np
import mageis_diffusion_curves
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime

def sinAlpha(alpha, A, n):
    """
    This function will return a value from the function A*sin(alpha)^n. 
    This is used for fitting the equatorial pitch angle distribution.
    """
    return A*np.sin(np.deg2rad(alpha))**n

if __name__ == '__main__':
    # A dictionary of times to analyze
    tBoundsDict = {'q':[datetime(2017, 3, 31, 11, 15, 0), datetime(2017, 3, 31, 11, 17, 0)], 'm1':[datetime(2017, 3, 31, 11, 17, 0), datetime(2017, 3, 31, 11, 17, 20)], 'm2':[datetime(2017, 3, 31, 11, 17, 10), datetime(2017, 3, 31, 11, 17, 20)],
    'bigOne':[datetime(2017, 3, 31, 11, 17, 13), datetime(2017, 3, 31, 11, 17, 18)],
    'smallOne':[datetime(2017, 3, 31, 11, 17, 9, 500000), datetime(2017, 3, 31, 11, 17, 10, 500000)]}
    tBounds = tBoundsDict['q']

    # Script parameters
    rb_id = 'A'
    instrument = 'LOW'
    alphaBins = np.arange(0, 180, 5)
    psdObj = mageis_diffusion_curves.PhaseSpaceDensity(rb_id, tBounds, instrument)
    psdObj.loadData()
    psdObj.calcPsd()
    psdObj.binPsdAlpha(alphaBins, psdErr = psdObj.psdErr)

    # Map the local pitch angle to the magnetic equator.
    alpha0Arr = psdObj.alpha0(psdObj.BeqOverB, alphaBins)

    # Filter the psd by valid pitch angles.
    ida = np.where(psdObj.meanPsd[0] != 0)[0]

    popt = np.nan*np.ones((psdObj.meanPsd.shape[0], 2), dtype = float)
    perr = np.nan*np.ones((psdObj.meanPsd.shape[0], 2), dtype = float)

    # Do the least squares fitting.
    for i in range(psdObj.meanPsd.shape[0]):
        popt[i, :], pcov = scipy.optimize.curve_fit(sinAlpha, alpha0Arr[ida], 
                    psdObj.meanPsd[i, ida], p0 = [10**-4, 1.5], sigma = None, 
                    absolute_sigma = False)
        perr[i, :] = np.sqrt(np.diag(pcov))
    print(popt, perr)

    # save the fit parameters and errors.
    np.save('../data/quiet_fit_popt', popt)
    np.save('../data/quiet_fit_perr', perr)

    fitAlphaArr = np.arange(10, 180)
    ### PLOTTING ###
    fig = plt.figure(figsize=(15, 10), dpi = 80, facecolor = 'white')
    plt.rcParams.update({'font.size': 15})
    gs = gridspec.GridSpec(1, 1)
    psdPlt = fig.add_subplot(gs[0, 0], facecolor='w')
    for i in range(psdObj.meanPsd.shape[0]):
        psdPlt.errorbar(alpha0Arr[ida], psdObj.meanPsd[i, ida], 
            label = '{} keV'.format(psdObj.Emid[i]), ls = 'None', marker = 'o',
            yerr = psdObj.meanPsdErr[i, ida])
        psdPlt.plot(fitAlphaArr, sinAlpha(fitAlphaArr, popt[i, 0], popt[i, 1]))
    psdPlt.set(yscale = 'log', xlabel = r'$\alpha_{eq}$', 
        ylabel = r'PSD $c^3/(cm \ MeV)^3$', xlim = (0, 180))

    psdPlt.set_title('RBSP-{} phase space density for {} \n {} - {} UT'.format(rb_id, 
            tBounds[0].date(), tBounds[0].strftime("%H:%M:%S"), 
            tBounds[1].strftime("%H:%M:%S")))
    psdPlt.legend()
    plt.show()
        
