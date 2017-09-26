# This script creates a matrix of Meredith style PSD plots. 
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import numpy as np
import os

import resonant_diffusion_curves
import mageis_diffusion_curves
import psd_fit

# Constants of nature
c = 3E8 # m/s

# Load MagEIS data, calculate the PSD, and bin it to equatorial pitch angle.
tBounds = [datetime(2017, 3, 31, 11, 17, 0), datetime(2017, 3, 31, 11, 17, 20)]
rb_id = 'A'
instrument = 'LOW'
alphaBins = np.arange(0, 181, 5)
psdObj = mageis_diffusion_curves.PhaseSpaceDensity(rb_id, tBounds, instrument)
psdObj.loadData()
psdObj.calcPsd()
#alphaBins = alphaBins[np.where(alphaBins < 180)[0]]
psdObj.binPsdAlpha(alphaBins, psdErr = psdObj.psdErr)
alpha0Arr = psdObj.alpha0(psdObj.BeqOverB, alphaBins)

# Get the extrapolated PSD, using the fit parameters.
dir_path = os.path.dirname(os.path.realpath(__file__)) # Script dir
popt = np.load(os.path.join(dir_path, 'data', 'quiet_fit_popt.npy'))
extPsd = np.nan*np.ones((7, len(alphaBins)), dtype = float) 

for e in range(7):
    extPsd[e, :] = psd_fit.sinAlpha(alphaBins, popt[e, 0], popt[e, 1])

# Find the index where the alpha0Arr jumps (largest pitch angle sampled)
ida = np.where(np.abs(np.convolve(alpha0Arr, [-1, 1])) > 10)[0][0]

# Draw straight lines representing the trapped edges of the real data
n_perp =  psdObj.p_perp(np.linspace(0, psdObj.Ehigh[-1]), alpha0Arr[ida-1]) 
n_parallel =  psdObj.p_parallel(np.linspace(0, psdObj.Ehigh[-1]), 
    alpha0Arr[ida-1])
s_perp =  psdObj.p_perp(np.linspace(0, psdObj.Ehigh[-1]), alpha0Arr[ida])
s_parallel =  psdObj.p_parallel(np.linspace(0, psdObj.Ehigh[-1]), 
    alpha0Arr[ida]) 

################
### PLOTTING ###
################
nRows = 2
nCols = 3
dCols = 10
fig = plt.figure(figsize=(11, 12), dpi = 80, facecolor = 'white')
plt.rcParams.update({'font.size': 15})
gs = gridspec.GridSpec(nRows, dCols*nCols+2)
gs.update(wspace=0, hspace=0.01)
axArr = np.nan*np.ones((nRows, nCols), dtype = object)
for (ir, ic), ax in np.ndenumerate(axArr):
    axArr[ir, ic] = fig.add_subplot(gs[ir, ic*dCols:(ic+1)*dCols], facecolor='k')
colorAx = fig.add_subplot(gs[:, -2:], facecolor='k')

nArr = np.array([[0.5E6, 0.5E6, 0.5E6], [5E6, 5E6, 5E6]])
mlatArr = np.array([[0, 10, 20], [0, 10, 20], [0, 10, 20]])
diffFraction = 0.1 # Fraction of the cyclotron frequency to draw the diffusion curves.

for i, ax in np.ndenumerate(axArr):
    # Draw the extrapolarted patches
    zzz, p = psdObj.drawPatches(alphaBins, ax = ax, psd = extPsd)
    # Draw the edges of the data and extrapolation
    ax.plot(n_perp, n_parallel, 'w:')
    ax.plot(s_perp, s_parallel, 'w:') 
    # Draw the PSD data.
    ax, p = psdObj.drawPatches(alpha0Arr, ax = ax)
    cb = plt.colorbar(p, label=r'PSD $c^3/(cm \ MeV)^3$', cax=colorAx) 

    # Resonant-diffusion parameters
    vParallel_res = c*np.linspace(0, -0.99, num = 1000)
    #mlat = 0
    L = 5.7
    #n0 = 0.5E6 # Density at the time

    n = nArr[i]
    mlat = mlatArr[i]

    # Draw resonance curves
    # w/w_ce = 0.1 
    vPerp_res = resonant_diffusion_curves.resCurveVperp(
        vParallel_res, 0.1*resonant_diffusion_curves.wce(mlat, L), n, mlat, L)
    pPerp_res, pParallel_res = resonant_diffusion_curves.p(
        vPerp_res, vParallel_res)
    ax.plot(pPerp_res, pParallel_res, 'g')
    ax.plot(pPerp_res, -pParallel_res, 'g')
    label01 = mlines.Line2D([], [], color='g', markersize=15, 
        label=r'$0.1 \ \Omega_{ce}$')

    # w/w_ce = 0.4 
    vPerp_res = resonant_diffusion_curves.resCurveVperp(
        vParallel_res, 0.4*resonant_diffusion_curves.wce(mlat, L), n, mlat, L)
    pPerp_res, pParallel_res = resonant_diffusion_curves.p(
        vPerp_res, vParallel_res)
    ax.plot(pPerp_res, pParallel_res, 'r')
    ax.plot(pPerp_res, -pParallel_res, 'r')
    label04 = mlines.Line2D([], [], color='r', markersize=15, ls='-',
        label=r'$0.4 \ \Omega_{ce}$')

    # w/w_ce = 0.6
    vPerp_res = resonant_diffusion_curves.resCurveVperp(
        vParallel_res, 0.6*resonant_diffusion_curves.wce(mlat, L), n, mlat, L)
    pPerp_res, pParallel_res = resonant_diffusion_curves.p(
        vPerp_res, vParallel_res)
    ax.plot(pPerp_res, pParallel_res, 'b')
    ax.plot(pPerp_res, -pParallel_res, 'b')
    label06 = mlines.Line2D([], [], color='b', markersize=15, ls='-',
        label=r'$0.6 \ \Omega_{ce}$')
        
    #ax.legend(loc=4, handles=[label01, label04, label06], fontsize=10)

    Earr = psdObj.Emid
    vParallel_diff = vParallel_res

    for e in Earr[::2]:
        vPerp_diff = resonant_diffusion_curves.diffCurveVperp(
            vParallel_diff, diffFraction*resonant_diffusion_curves.wce(mlat, L), 
            n, mlat, L, e)
        pPerp_diff, pParallel_diff = resonant_diffusion_curves.p(vPerp_diff, 
            vParallel_diff)
        ax.plot(pPerp_diff, pParallel_diff, 'c--')
        ax.plot(pPerp_diff, -pParallel_diff, 'c--')
            
    ax.set(aspect = 'equal', ylim = (-1.2, 1.2), xlim = (0, 1.2))
        

    paramsStr = r'$\lambda = ${}$^{{\circ}}$ $n = {} \ e^-/cm^{{3}}$'.format(
        round(mlat), round(n/1E6, 2))
    ax.text(0.5, 0.95, paramsStr, horizontalalignment='center', 
        verticalalignment='center', transform=ax.transAxes, color = 'w', fontsize = 15)
        
# Turn off the subplot's ticklabels
axShape = axArr.shape
for ax in axArr[:, 1]:
    ax.set_yticklabels([])
for ax in axArr[:, -1]: # This is weird why it does not work with the loop above.
    ax.set_yticklabels([])
for ax in axArr[0, :]:
    ax.set_xticklabels([])
        
fig.suptitle('RBSP-{} phase space density for {} \n {} - {} UT'.format(rb_id, 
    tBounds[0].date(), tBounds[0].strftime("%H:%M:%S"), 
    tBounds[1].strftime("%H:%M:%S")))
fig.text(0.5, 0.02, r'$p_{\perp}/m_e c$', ha='center')
fig.text(0.02, 0.5, r'$p_{\parallel}/m_e c$', va='center', rotation='vertical')
        
gs.tight_layout(fig, rect=[0.02, 0.03, 1, 0.95])
plt.show()
