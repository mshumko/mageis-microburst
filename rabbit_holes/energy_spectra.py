# This script will calculate the energy spectra of the microbursts.

import sys
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp/')
sys.path.insert(0, '/home/mike/research/mission-tools/misc/')
sys.path.insert(0, '/home/mike/research/mageis-microburst-stats/detection/')
import plot_mageis as plot_mageis_lib
import obrienBaseline
import exp_fit

SC_ID = 'A'
SMOOTH = None
INSTRUMENT = 'LOW'
#A_WIDTH = 0.5
#BURST_THRESH = 100
DATE = datetime(2017, 3, 31)
T_BOUNDS = [datetime(2017, 3, 31, 11, 17, 0), datetime(2017, 3, 31, 11, 17, 16)]
T_LIM = [datetime(2017, 3, 31, 11, 17, 9), datetime(2017, 3, 31, 11, 17, 13)]
MU_BOUNDS = [[datetime(2017, 3, 31, 11, 17, 9, 707936),
            datetime(2017, 3, 31, 11, 17, 9, 864081)], 
            [datetime(2017, 3, 31, 11, 17, 10, 84116),
            datetime(2017, 3, 31, 11, 17, 10, 473371)], 
            [datetime(2017, 3, 31, 11, 17, 10, 750628),
            datetime(2017, 3, 31, 11, 17, 11, 117222)], 
            [datetime(2017, 3, 31, 11, 17, 11, 376817),
            datetime(2017, 3, 31, 11, 17, 11, 768979)]]
BASELINE = True
FIT_ENERGY = np.arange(4)

if BASELINE:
    BASELINE_STR = 'baseline subtracted'
else:
    BASELINE_STR = ''
fig, ax = plt.subplots(2, sharex=False, figsize=(8,8))
fluxObj = plot_mageis_lib.magEISspectra(
    SC_ID, DATE, dataLevel = 3)
fluxObj.tBounds = T_BOUNDS
fluxObj.loadMagEIS(instrument=INSTRUMENT, highrate=True)
fluxObj.get_resolved_flux(smooth=SMOOTH)

if BASELINE:
    baselineFlux = obrienBaseline.baselineMultiprocessWrapper(
        fluxObj.flux, 1, 10, cadence=11/1000)
    for i in baselineFlux: 
        fluxObj.flux[:, int(i)] -= baselineFlux[i]

# Calculate everaged flux over the microburst
jBar = np.nan*np.ones((7, 4), dtype=float) #nE by nMuBurst
E0 = np.nan*np.ones(4)
J0 = np.nan*np.ones(4)
for i in range(4):
    validInd = np.where((fluxObj.times > MU_BOUNDS[i][0]) & 
        (fluxObj.times <= MU_BOUNDS[i][1]))[0]
    jBar[:, i] = np.sum(fluxObj.flux[validInd, :], axis=0)/len(validInd)
    #Calculate the exponential fit
    fitDict = exp_fit.fit_exponent(
        np.array(fluxObj.Emid)[FIT_ENERGY], jBar[FIT_ENERGY, i])
    E0[i] = fitDict['E0']
    J0[i] = fitDict['J0']

# Plot flux
for ee in range(7):
    validF = np.where(fluxObj.flux[:, ee] > 0)[0]
    ax[0].plot(fluxObj.times[validF], fluxObj.flux[validF, ee])

# Plot jBar and fit, and labels to show which flux was used in the analysis.
colors = ['r', 'b', 'g', 'k']
Earr = np.arange(25, 150)
for e in range(4):
    ax[0].axvline(MU_BOUNDS[e][0], color=colors[e], ls='--')
    ax[0].axvline(MU_BOUNDS[e][1], color=colors[e], ls='--')

for i, j in enumerate(jBar.T):
    ax[1].scatter(fluxObj.Emid, j, c=colors[i]) 
    # Plot fit, and annotate plot
    ax[1].plot(Earr, J0[i]*np.exp(-Earr/E0[i]), c=colors[i])
    ax[1].text(0.4, 0.95-0.05*i, r'E0={} keV | J0={:.2E} $(cm^2 \ sr \ s \ keV)^{{-1}}$'.format(round(E0[i]), round(J0[i])), 
        color=colors[i], va='top',transform=ax[1].transAxes)
    
    
ax[0].set(yscale='log', ylabel=r'J $(cm^2 \ sr \ s \ keV)^{-1}$',
    xlim=T_LIM, ylim=(10**4, 5E5), xlabel='UTC', 
    title='MagEIS microburst energy spectra {}'.format(BASELINE_STR))
ax[1].set(yscale='log', ylabel=r'$\bar{J}$ $(cm^2 \ sr \ s \ keV)^{-1}$',
    xlabel='Energy (keV)')
plt.show()
