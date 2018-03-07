# This script will fit the RBSPICE pitch angle distribution with a
# sine function.

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import spacepy.pycdf
from datetime import datetime, timedelta

def sin_fit(t, A, w, phi, b):
    """ Sine function to fit. """
    return A*np.sin(t*w + phi)+b

def load_esrhelt(fname, tRange=None):
    """
    This function will load in the ESRHELT data and filter the epoch,
    pitch angle, and electron flux arrays.
    """
    d = spacepy.pycdf.CDF(fname).copy()
    d['Epoch'] = np.array(d['Epoch'])
    # If user provides a time range, filter the relevant keys.
    if tRange is not None:
        idT = np.where((tRange[0] < d['Epoch'][:]) & 
                    (tRange[1] > d['Epoch'][:]))[0]
        d['Epoch'] = d['Epoch'][idT]
        d['FEDU_Alpha'] = np.array(d['FEDU_Alpha'])[idT, :]
        d['FEDU'] =  np.array(d['FEDU'])[idT, :, :]
    return d

def fit_alpha(t, a, p0=None):
    """
    This function is a wrapper for scipy.optimize.curve_fit to fit 
    a sine curve to the pitch angle evolution.
    """
    popt, pcov = scipy.optimize.curve_fit(sin_fit, t, a, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    return popt, pcov

if __name__ == '__main__':
    fName = ('/home/mike/research/rbsp/data/rbspice/rbspa/'
            'rbsp-a-rbspice_lev-3_ESRHELT_20170331_v1.1.9-01.cdf')
    START_TIME = datetime(2017, 3, 31, 11, 16, 30)
    END_TIME = datetime(2017, 3, 31, 11, 17, 30)
    # START_TIME = datetime(2017, 3, 31, 11, 17, 10)
    # END_TIME = datetime(2017, 3, 31, 11, 17, 21)

    d = load_esrhelt(fName, tRange=[START_TIME, END_TIME])
    # Proce
    t = [(t - d['Epoch'][0]).total_seconds() for t in d['Epoch']]
    fitT = [START_TIME + timedelta(seconds=i) for i in 
            np.arange(0, int((END_TIME-START_TIME).total_seconds()), 0.25)]
    fitTsec = np.array([(t - fitT[0]).total_seconds() for t in fitT])

    p0 = np.array([[20, 2*np.pi/11, np.pi/3, 100], 
                   [40, 2*np.pi/11, np.pi/3, 100],
                   [65, 2*np.pi/11, np.pi/2, 100],
                   [65, 2*np.pi/11, np.pi/2, 100],
                   [65, 2*np.pi/11, np.pi/2, 100],
                   [65, 2*np.pi/11, np.pi/2, 100]])
    # popt, pcov = fit_alpha(t, d['FEDU_Alpha'][:, 0], p0=p0)
    # print(popt)

    colors = ['r', 'b', 'g', 'c', 'k', 'y']
    for (i, aArr) in enumerate(d['FEDU_Alpha'].T):
        plt.scatter(d['Epoch'], aArr, c=colors[i])
        popt, pcov = fit_alpha(t, d['FEDU_Alpha'][:, i], p0=p0[i, :])
        plt.plot(fitT, sin_fit(fitTsec, *popt), color=colors[i])

    plt.xlim(START_TIME, END_TIME)
    plt.title('Sine fit to RBSPICE PA data')
    plt.xlabel('UTC')
    plt.ylabel(r'$\alpha_L$')
    plt.show()