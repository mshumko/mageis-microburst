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

# Script parameters
tBounds = [datetime(2017, 3, 31, 11, 17, 1), datetime(2017, 3, 31, 11, 17, 12)]
rb_id = 'A'
instrument = 'LOW'
mlat0 = -20 # Degrees
n0 = 0.5E6 # e-/cm^3
L = 5.7
a = 1 # Electron number density power law coefficient.
mlats = [0, 20]
# Fraction of the cyclotron frequency to draw the diffusion curves.
diffFraction = 0.4 

dataAlphaBins = np.arange(0, 180, 5)
extrapAlphaBins = np.arange(30, 150, 5)
vmax = 10**-1
vmin = 10**-4
poptFnameDict = {'fit':'psd_fit_extrapolation_popt_111500_111700.npy',
            '1':'psd_fit_fudged_n_1_popt_111500_111700.npy',
            '2':'psd_fit_fudged_n_2_popt_111500_111700.npy',
            '4':'psd_fit_fudged_n_4_popt_111500_111700.npy'}

# Load MagEIS data, calculate the PSD, and bin it to equatorial pitch angle.
psdObj = mageis_diffusion_curves.PhaseSpaceDensity(rb_id, tBounds, instrument)
psdObj.loadData()
psdObj.calcPsd()
psdObj.binPsdAlpha(dataAlphaBins, psdErr = psdObj.psdErr)
alpha0Arr = psdObj.alpha0(psdObj.BeqOverB, dataAlphaBins)

# Get the extrapolated PSD, using the fit parameters from files in the
# poptFnameDict
dir_path = os.path.dirname(os.path.realpath(__file__)) # Script dir
popDict = {}
for key in poptFnameDict.keys():
    popDict[key] = np.load(os.path.join(dir_path, 'data', poptFnameDict[key]))
    
extPsd = np.nan*np.ones((7, len(extrapAlphaBins), len(poptFnameDict.keys())),
    dtype=float) 

# Sort the dictionary keys so that the first key is the fit, and then the 
# extrapolations in ascending order
dictKeys = sorted(poptFnameDict.keys())
dictKeys.remove('fit')
dictKeys.insert(0, 'fit')

# Now save the extrapolated PSD values over the extrapolated pitch angles.
for i_n, nKey in enumerate(dictKeys):
    for e in range(7):
        extPsd[e, :, i_n] = psd_fit.sinAlpha(extrapAlphaBins, popDict[nKey][e, 0],
            popDict[nKey][e, 1])

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
nCols = 4
dCols = 10
fig = plt.figure(figsize=(13, 11), dpi = 80, facecolor = 'white')
plt.rcParams.update({'font.size': 15})
gs = gridspec.GridSpec(nRows, dCols*nCols+2)
gs.update(wspace=0, hspace=0.05)
axArr = np.nan*np.ones((nRows, nCols), dtype = object)
for (ir, ic), ax in np.ndenumerate(axArr):
    axArr[ir, ic] = fig.add_subplot(gs[ir, ic*dCols:(ic+1)*dCols], facecolor='k')
colorAx = fig.add_subplot(gs[:, -2:], facecolor='k')

mlatArr = np.meshgrid(np.ones(nCols), mlats)[1]
#mlatArr = np.array([[0, 0, 0, 0], [10, 10, 10, 10], [20, 20, 20, 20]])

for i, ax in np.ndenumerate(axArr):
    # Draw the extrapolarted patches
    zzz, p = psdObj.drawPatches(extrapAlphaBins, ax=ax, psd=extPsd[:, :, i[1]], 
        vmin=vmin, vmax=vmax, cMapLog=False)
    # Draw the edges of the data and extrapolation
    ax.plot(n_perp, n_parallel, 'w:')
    ax.plot(s_perp, s_parallel, 'w:') 
    # Draw the PSD data.
    ax, p = psdObj.drawPatches(alpha0Arr, ax = ax, vmin=vmin, vmax=vmax,
        cMapLog=False)
    cb = plt.colorbar(p, label=r'PSD $c^3/(cm \ MeV)^3$', cax=colorAx) 

    # Resonant-diffusion parameters
    vParallel_res = c*np.linspace(0, -0.99, num = 1000)
    mlat = mlatArr[i]
    nn = resonant_diffusion_curves.n_e(n0, mlat0, mlat, a)

    # Draw resonance curves
    # w/w_ce = 0.1 
    vPerp_res = resonant_diffusion_curves.resCurveVperp(
        vParallel_res, 0.1*resonant_diffusion_curves.wce(mlat, L), nn, mlat0, 
        mlat, L, a=a)
    pPerp_res, pParallel_res = resonant_diffusion_curves.p(
        vPerp_res, vParallel_res)
    ax.plot(pPerp_res, pParallel_res, 'g')
    ax.plot(pPerp_res, -pParallel_res, 'g')
    label01 = mlines.Line2D([], [], color='g', markersize=15, 
        label=r'$0.1 \ \Omega_{ce}$')

    # w/w_ce = 0.4 
    vPerp_res = resonant_diffusion_curves.resCurveVperp(
        vParallel_res, 0.4*resonant_diffusion_curves.wce(mlat, L), nn, mlat0,
        mlat, L, a=a)
    pPerp_res, pParallel_res = resonant_diffusion_curves.p(
        vPerp_res, vParallel_res)
    ax.plot(pPerp_res, pParallel_res, 'r')
    ax.plot(pPerp_res, -pParallel_res, 'r')
    label04 = mlines.Line2D([], [], color='r', markersize=15, ls='-',
        label=r'$0.4 \ \Omega_{ce}$')

    # w/w_ce = 0.6
    vPerp_res = resonant_diffusion_curves.resCurveVperp(
        vParallel_res, 0.6*resonant_diffusion_curves.wce(mlat, L), nn, 
        mlat0, mlat, L, a=a)
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
            n0, mlat0, mlat, L, e, a=a)
        pPerp_diff, pParallel_diff = resonant_diffusion_curves.p(vPerp_diff, 
            vParallel_diff)
        ax.plot(pPerp_diff, pParallel_diff, 'c--')
        ax.plot(pPerp_diff, -pParallel_diff, 'c--')
            
    ax.set(aspect='equal', ylim=(-1.2, 1.2), xlim=(0, 1.2))
        
    mlat_str = r'$\lambda={0}^{{\circ}} \ ({1} \ cm^{{-3}})$'.format(
        round(mlat), round(resonant_diffusion_curves.n_e(
            n0, mlat0, mlat, a=a)*1E-6, 2))
    n_extrap = r'$n_{{ext}}={0}$'.format(dictKeys[i[1]])
    n_e_str = r'$n_{{e}} = {0}$'.format(
        round(resonant_diffusion_curves.n_e(n0, mlat0, mlat, a=a)*1E-6, 2))
    ax.text(1, 1, mlat_str + ', '+ n_extrap, horizontalalignment='right', 
        verticalalignment='top', transform=ax.transAxes, color='k', 
        fontsize=12, bbox=dict(facecolor='w', boxstyle="round"))
        
# Turn off the subplot's ticklabels   
for ii, ax in np.ndenumerate(axArr):
    if ii[1] == 0: # Can't figure out how to slide the array correctly...
        continue
    axArr[ii].set_yticklabels([])
for ii, ax in np.ndenumerate(axArr[:-1, :]):
    axArr[ii].set_xticklabels([])
for ax in axArr[-1, :]: # Throw on more x ticks.
    ax.set_xticks([0, 0.5, 1])
        
fig.suptitle('RBSP-{} phase space density for {} \n {} - {} UT'.format(rb_id, 
    tBounds[0].date(), tBounds[0].strftime("%H:%M:%S"), 
    tBounds[1].strftime("%H:%M:%S")))
fig.text(0.5, 0.01, r'$p_{\perp}/m_e c$', ha='center')
fig.text(0.01, 0.5, r'$p_{\parallel}/m_e c$', va='center', rotation='vertical')

gs.tight_layout(fig, h_pad = 0, w_pad = 0, rect=[0.02, 0.03, 1, 0.95])
plt.savefig('resonance_diffusion_matrix_{0}_{1}_n0_{2}_a_{3}.png'.format(
    tBounds[0].strftime("%H:%M:%S"), tBounds[1].strftime("%H:%M:%S"), 
    str(n0*1E-6).replace('.', '-'), a), dpi=100)
plt.show()
