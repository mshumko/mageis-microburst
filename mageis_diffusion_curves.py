# Code to create the diffusion surfaces from the magEIS data.

### Imports ###
import numpy as np
import sys
import copy
from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import operator
from matplotlib.patches import Polygon
import matplotlib.collections

sys.path.insert(0, '/home/mike/Dropbox/0_grad_work/mission_tools')
import plot_emfisis_spectrogram
import plot_mageis_spectra
import psd_fit

### Constants of nature ###
c = 3E10 # cm/s
mu_0 = 1.26E-6 # H·m−1 or N·A−2
eps0 = 8.85E-12 # F/m
m_p = 1.7E-27 # kg
m_e = 9.1E-31 # kg
q_e = -1.6E-19 # C
Erest = 511 # keV

# Define relativistic $\beta$ and $\gamma$, functions of kinetic energy for electrons
beta = lambda Ek: np.sqrt(1-(Ek/511+1)**(-2))
gamma = lambda Ek:np.sqrt(1-beta(Ek)**2)**(-1/2)

# Define the dipole magnetic field
B0 = 31.2E-6 # Tesla from Schultz and Lanzerotti, MagB is from Eq. 1.23
magB = lambda λ, L: (B0/L**3)*np.sqrt(1 + 3*np.power(np.sin(np.deg2rad(λ)), 2))/np.cos(np.deg2rad(λ))**6

# Frequency, magnitude of k, and number density definitions
wce = lambda λ, L: np.abs(q_e)*magB(λ, L)/m_e
n_e = lambda n0, λ = None: n0 # Electron number density. Currently constant, but can assume a complex function.
wpe = lambda n0, λ = None: np.sqrt(4*np.pi*n_e(n0, λ)*q_e**2/(m_e*eps0))
magk = lambda w, n0, λ, L: (w/c)*np.sqrt(1 - wpe(n0, λ)**2/(w*(w - wce(λ, L))))

# Define the resonance curve
def resCurveVperp(vParallel, w, n0, λ, L, n = 1):
    """
    This function defines the perpendicular velocity of a resonant particle with n = 1
    """
    A = -(c*(w + magk(w, n0, λ, L)*vParallel)/(wce(λ, L)*n))**2
    return np.sqrt(A + c**2 - vParallel**2)

# Define Alfven speed, alpha parameter from Summers 1998 paper, and phase velocity.
# omega must be normalized to the elctron gyrofrequency, and the B field in 
# units of nT and electron number density, n in #/cc
v_a = lambda B, n: 1e-9*B/np.sqrt((1e6)*mu_0*n*m_p)
summersAlpha = lambda B, n: (m_p/m_e)*(v_a(B, n)/(c/1e2))**2
u_ph = lambda omega, B, n: np.sqrt(summersAlpha(B, n)*omega*(1 - omega))

class PhaseSpaceDensity(plot_mageis_spectra.magEISspectra): # Utilize inheritance
    def __init__(self, rb_id, tRange, instrument, **kwargs):
        """
        NAME:    PhaseSpaceDensity(rb_id, times, instrument)
        USE:     This class calculates the relativistic phase space density.
        INPUT:   REQUIRED:
                   rb_id: either 'A' or 'B', 
                   times: a time range to look for a spacecraft rotation. Use 
                   times = [t0, t0 + 12 s] for normal VAP spin rate.
                   instrument: The magEIS instrument, 'LOW', 'M35', 'M75', 'HIGH'
                 OPTIONAL:
                    ### FINISH WRITING THIS LATER ###

        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-08-08
        """
        # Save the physicist some time in case he messes up or tries to utilize
        # functionality that has not been implemented.
        assert rb_id.upper() in ['A', 'B'], 'ERROR: Incorrect RBSP id!'
        assert instrument.upper() == 'LOW', ('ERROR: Only the low instrument'
            ' detector parameters accounted for!')
        assert tRange[0] > datetime(2013, 7, 12), ('ERROR: Only the low instrument'
            ' detector parameters accounted for!')
        
        # Define magEIS detector constants
        if rb_id.upper() == 'A':
            self.Emid = [34, 54, 78, 108, 143, 182, 223] # keV
            self.Elow = [29, 46, 68, 95, 126, 164, 206] # keV
            self.Ehigh = [41, 66, 92, 126, 164, 204, 247] # keV
            # Units of (keV cm^2 sr)
            self.G0dE = [4.13E-2, 5.73E-2, 6.056E-2, 6.88E-2, 7.35E-2, 
                6.90E-2, 5.98E-2]
            self.Ebins = [29, 41, 66, 92, 126, 164, 204, 247]
            
        if rb_id.upper() == 'B':
            self.Emid = [32, 51, 74, 101, 132, 168, 208] # keV
            self.Elow = [27, 43, 63, 88, 117, 152, 193] # keV
            self.Ehigh = [39, 63, 88, 117, 150, 188] # keV
            # Units of (keV cm^2 sr)
            self.G0dE = [4.33E-2, 5.41E-2, 5.926E-2, 6.605E-2, 6.460E-2,
                6.23E-2, 5.96E-2]
            self.Ebins = [27, 39, 63, 88, 117, 150, 188]
                
        self.tSpin = 10.9 # Spin period in seconds.
        self.rb_id = rb_id
        self.tRange = tRange
        self.instrument = instrument

        # Load the magEIS plotting object
        dataLevel = kwargs.get('dataLevel', 3)
        super().__init__(rb_id, self.tRange[0], dataLevel = dataLevel)
        self.tBounds = self.tRange
        
        return
                
    def loadData(self, dataLevel = 3):
        """
        NAME:    loadMagEIS
        USE:     This function loads in magEIS data and 
                 magnetic ephemeris.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    remapAlpha: fold pitch angles from 0 to 360 to
                                0 to 180 degrees.
                    datalevel: Which level of data to load.

        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-09-07
        """
        # Widen the user time range so that a) there will be at least one magEphem
        # data point and b) so that it can be flattened, without any part of the
        # user's time range cut off.        
        self.tBounds = [self.tRange[0] - timedelta(minutes = 2), 
            self.tRange[1] + timedelta(minutes = 2)]
        
        self.loadMagEIS(instrument = self.instrument, highrate = True)


        
        # Calculate the flattened time and pitch angle arrays.
        self.resolveSpinTimes(flattenTime = True)
        self.flatAlpha = self.magEISdata['HighRate_Alpha360'].flatten()

        self.loadMagEphem()        
        self.BeqOverB = 1/np.mean(self.magEphem['BoverBeq'])
        
        # Lastly, filter the flattened times and pitch angles.
        validIdt = np.where((self.times >= self.tRange[0]) & 
            (self.times <= self.tRange[1]) & (self.flatAlpha != -1E31))[0]
        self.times = self.times[validIdt]
        self.flatAlpha = self.flatAlpha[validIdt]

        # Equatorial pitch angle array
        self.alpha0Arr = self.alpha0(self.BeqOverB, self.flatAlpha)
     
        # Convert the count data into a 2d array with nEnergy, nTime
        self.magEISdata['HighRate'] = np.array([
            self.magEISdata['HighRate'][:, :, E_ch].flatten()[validIdt]
            for E_ch in range(8)])
        return
        
    def calcPsd(self, lowerAlpha = None, upperAlpha = None):
        """
        This function calculates the phase space density, and it filters by
        pitch angle and invalid indicies
        """
        #alpha = self.magEISdata['HighRate_Alpha360']
        #validIda = np.where((alpha != -1E31))[0]
        
        if lowerAlpha is not None:
            lowerIda = np.where(self.flatAlpha >= lowerAlpha)[0]
        else:
            lowerIda = range(len(self.flatAlpha))
            
        if upperAlpha is not None:
            upperIda = np.where(self.flatAlpha <= upperAlpha)[0]
        else:
            upperIda = range(len(self.flatAlpha))
            
        validIda = np.array(list(set(lowerIda) & set(upperIda)))

        # Create data phase space density in units of c^3/(cm^6 s^3 meV^3)
        # In addition, calculate error, assuming Poisson statistics.
        psd = 1E9*np.array([self.f(self.j(i)[0][validIda], self.Emid[i], 
            jErr = self.j(i)[1][validIda]) for i in range(7)]) 
        self.psd = psd[:, 0, :]
        self.psdErr = psd[:, 1, :]
        return self.psd, self.psdErr

    def calcDiffusionCurveConstantU(self, v_p,  u_0, v_0):
        """
        This function returns the perpendicular velocity from the diffusion
        equation given in Eq. 6 in Summers et al. 1998. 

        Input velocity units must be divided by c!

        THIS FUNCTION ASSUMES THAT THE WAVE PHASE VELOCITY IS CONSTANT!
        """
        numerator = ( -(1 - (u_0*v_0)**2)*v_p**2 + 
            2*u_0*(1 - v_0**2)*v_p + v_0**2 + u_0**2)
        denomimator = 1 - u_0**2
        return np.sqrt(numerator/denomimator)

    def calcDiffusionCurveConstantUWrapper(self, v_parallel, E, B, n, omega):
        """
        Energy must be in keV!
        """
        u = u_ph(omega, B, n)
        v_0 = beta(E)
        v_perp = self.calcDiffusionCurveConstantU(v_parallel,  u, v_0)
        
        #filter out the nan's
        validV = np.where(np.logical_not(np.isnan(v_perp)))[0]
        return v_perp[validV], v_parallel[validV]

    def calcResonanceCurves(self, ):
        """
        
        """

        return
        
    def binPsdAlpha(self, binEdges, psdErr = None, zeroPsdFill = 0, 
            remapAlpha = True):
        """
        This function will bin the phase space density using the alpha bin edges
        array.
        
        remapAlpha kwarg will bin the pitch angles assuming the relevant angles
        are from 0 to 180.
        """
        if remapAlpha:
            # Fold binned and data alpha from 0 to 360 to 0 to 180.
            binEdges = binEdges[np.where((binEdges >= 0) & 
                (binEdges <= 180))[0]]
            self.flatAlphaFolded = copy.copy(self.flatAlpha)
            self.flatAlphaFolded[self.flatAlpha > 180] = (
                -self.flatAlpha[self.flatAlpha > 180] % 180)
        self.meanPsd = np.zeros((self.psd.shape[0], len(binEdges) - 1))
        if psdErr is not None:
            self.meanPsdErr = np.zeros((self.psd.shape[0], len(binEdges) - 1))
        self.alphaBinMid = np.convolve([0.5, 0.5], binEdges, mode = 'valid')
        
        # Loop over pitch angle bins.
        for a in range(len(binEdges[:-1])):
            validAlpha = np.where((self.flatAlphaFolded >= binEdges[a]) & 
                                  (self.flatAlphaFolded < binEdges[a+1]))[0] 

            # If no data points found at that equatorial pitch angle, fill in
            # the psd with a dummy variable zeroPsdFill (This avoids nan's for
            # later plotting).
            if len(validAlpha) == 0:
                self.meanPsd[:, a] = zeroPsdFill*np.ones((self.psd.shape[0]))
                continue
            # Calculate mean phase space density in that bin.
            self.meanPsd[:, a] = np.mean(self.psd[:, validAlpha], axis = 1)
            if psdErr is not None:
                self.meanPsdErr[:, a] = np.sqrt(np.sum(
                    psdErr[:, validAlpha]**2, axis = 1))/len(validAlpha)
                #np.mean(self.psd[:, validAlpha], axis = 1)
        
        return self.alphaBinMid, self.meanPsd
        
        
    def drawPatches(self, eqPitchAngleBins, energyBins = None, ax = None, 
            psd = None):
        # Create the verticies from the pitch angle bins.
        self.__makePatchVerticies__(eqPitchAngleBins, energyBins, psd = psd)   
        # Use lists here to filter out the pitch angle bins with no samples.   
        patches = []
        c = []

        for i in range(self.verticies.shape[2]):
            c.append(self.c[i])
            patches.append(Polygon(self.verticies[:, :, i]))

        p = matplotlib.collections.PatchCollection(patches)
#        p.set_cmap('rainbow')
        p.set_cmap('plasma')
        p.set_array(np.array(c))
        p.autoscale()
        p.set_norm(matplotlib.colors.LogNorm())
        p.set_clim(vmin=10**-4, vmax=10**-1)
        ax.add_collection(p)
        #fig.colorbar(p, ax=ax)
        return ax, p
        
    def __makePatchVerticies__(self, eqPitchAngleBins, energyBins = None,
            ignoreBoxHeght = 0.2, psd = None):
        """
        
        """
        if energyBins is None:
            energyBins = self.Ebins
            
        if psd is None:
            psd = self.meanPsd
            
        aa, ee = np.meshgrid(eqPitchAngleBins, energyBins)
        # Create a 3d array with dimentions (nVerticies, nCoords, nQuadralaterals)
        #self.verticies = np.nan*np.ones((4, 2, np.product(aa.shape)-2))
        self.verticies = np.nan*np.ones((4, 2, 0))
        
        # Make a 1d array of colors, analogy to the 3rd dimention of self.verticies
        self.c = np.array([])
        pPerp = self.p_perp(ee, aa)
        pParallel = self.p_parallel(ee, aa)
        
        # i index runs over energy, j index runs over pitch angle
        for i in range(pPerp.shape[0]-1):
            for j in range(pPerp.shape[1]-1):
                if np.isnan(psd[i, j]):
                    continue
                if (pParallel[i, j] - pParallel[i+1, j+1]) > ignoreBoxHeght:
                    # This is to avoid drawing large rectangles at local
                    # pitch angle of ~90 (Ambiguous if it maps to norward or
                    # southward pitch angles.
                    continue
                # Create the verticies
                vertex = np.array(
                    [[pPerp[i, j], pParallel[i, j]], 
                    [pPerp[i+1, j], pParallel[i+1, j]], 
                    [pPerp[i+1, j+1], pParallel[i+1, j+1]], 
                    [pPerp[i, j+1], pParallel[i, j+1]]] )
                vertex = np.expand_dims(vertex, 2)
                self.verticies = np.concatenate((self.verticies, vertex), 
                    axis = -1)
                self.c = np.append(self.c, psd[i, j])
        return self.verticies

    
    ### Define phase space density (PSD, f), relativistic momentum, and e flux.
    # @staticmethod decorator does not work with lambda functions! Need to do
    # staticmethod(lambda ...)
    
    # Drew's derivation of the phase space density
    @staticmethod
    def f(j, Ek, Erest = Erest, jErr = None):
        psd = np.divide(j, c*Ek*(Ek + 2*Erest))
        if jErr is None:
            return psd
        else:
            psdErr = np.divide(jErr, c*Ek*(Ek + 2*Erest))
            return psd, psdErr

    # Calculate the electron flux for a given energy channel.
    def j(self, E_ch):
        j = self.magEISdata['HighRate'][E_ch, :]/self.G0dE[E_ch]
        err = np.sqrt(self.magEISdata['HighRate'][E_ch, :])/self.G0dE[E_ch]
        return j, err
    
    @staticmethod 
    def p_perp(E, alpha0, Erest = Erest):
        return np.sqrt(E*(E+2*Erest))*np.sin(np.deg2rad(alpha0)
            )*np.sign(alpha0)/Erest
            
    @staticmethod 
    def p_parallel(E, alpha0, Erest = Erest):
        return np.sqrt(E*(E+2*Erest))*np.cos(np.deg2rad(alpha0)
            )*np.sign(alpha0)/Erest

    # Relativistic momentum from the kinetic energy output units of keV*s/m
    p = staticmethod(lambda T, Erest = Erest: np.sqrt(T**2 + 2*T*Erest)/c)

    @staticmethod
    def alpha0(*args):
        """
        Equatorial pitch angle
        
        First argument is B_eq, second is B_sc and third is alpha_sc. Or can provide
        B_eq/B_sc directly and alpha_sc
        """
        assert len(args) in [3, 2], "Number of arguments must be 2 or 3!"
        
        if len(args) == 3:
            a0 = np.arcsin(np.sqrt(np.divide(args[0], args[1]))*np.sin(np.deg2rad(args[2])))
        elif len(args) == 2:
            a0 = np.arcsin(np.sqrt(args[0])*np.sin(np.deg2rad(args[1])))
            
        ida = np.where(args[-1] > 90)[0]
        a0[ida] = np.pi - a0[ida] #(-a0[ida] % np.pi/2) 
        return np.rad2deg(a0)

if __name__ == '__main__':
    tBoundsDict = {'q':[datetime(2017, 3, 31, 11, 15, 0), datetime(2017, 3, 31, 11, 17, 0)], 'm1':[datetime(2017, 3, 31, 11, 17, 0), datetime(2017, 3, 31, 11, 17, 20)], 'm2':[datetime(2017, 3, 31, 11, 17, 10), datetime(2017, 3, 31, 11, 17, 20)],
    'bigOne':[datetime(2017, 3, 31, 11, 17, 13), datetime(2017, 3, 31, 11, 17, 18)],
    'smallOne':[datetime(2017, 3, 31, 11, 17, 9, 500000), datetime(2017, 3, 31, 11, 17, 10, 500000)]}
    tPeriod = 'm1'

    tBounds = tBoundsDict[tPeriod]
    
    if tPeriod == 'q' or tPeriod == 'm1': # Only draw the interpolated PSD for the quiet time
        extPSD = True
        save_description = '_fit_extrapolated'
    else:
        extPSD = False
        save_description = ''
        
    rb_id = 'A'
    instrument = 'LOW'
    alphaBins = np.arange(0, 181, 5)
    psdObj = PhaseSpaceDensity(rb_id, tBounds, instrument)
    psdObj.loadData()
    psdObj.calcPsd()
    alphaBins = alphaBins[np.where(alphaBins < 180)[0]]
    psdObj.binPsdAlpha(alphaBins, psdErr = psdObj.psdErr)
    plotPitchAngles = False
    plotMeredithPlot = True
    saveMeredithPlt = True
    drawDiffusionCurves = False
    drawResonanceCurves = False
    drawEqualE = False
    vmin = 1E-4
    vmax = 1E-140*7
    
    da = 10 # smooth every da points

    # Get equatorial alpha0array
    alpha0Arr = psdObj.alpha0(psdObj.BeqOverB, alphaBins)

    ### PLOTS ###
    fig = plt.figure(figsize=(15, 10), dpi = 80, facecolor = 'white')
    plt.rcParams.update({'font.size': 15})

    if plotPitchAngles:
        gs = gridspec.GridSpec(1, 1)
        psdPlt = fig.add_subplot(gs[0, 0], facecolor='k')
        for i in range(psdObj.meanPsd.shape[0]):
        #    plt.plot(psdObj.alphaBinMid, psdObj.meanPsd[i])
            psdPlt.errorbar(psdObj.alphaBinMid, psdObj.meanPsd[i], 
                yerr = psdObj.meanPsdErr[i], label = '{} keV'.format(psdObj.Emid[i]))
        psdPlt.set_yscale('log')
        psdPlt.legend(loc = 1)

        psdPlt.set_title('RBSP-{} '.format(rb_id) + r'PSD($\alpha_0$, E)' + 
            '\n {}, {} - {} UT'.format(tBounds[0].date(), tBounds[0].strftime("%H:%M:%S"), 
            tBounds[1].strftime("%H:%M:%S")))
        psdPlt.set_ylabel(r'PSD $c^3/(cm \ MeV)^3$')
        psdPlt.set_xlabel(r'$\alpha_0$')
        plt.show()

    # Now create the arrays for momentum space energy plots
    p_e_perp = np.array([psdObj.p_perp(psdObj.Emid[i], np.linspace(1, 180)) for i in range(7)])
    p_e_parallel = np.array([psdObj.p_parallel(psdObj.Emid[i], np.linspace(1, 180)) for i in range(7)])

    # Define diffusion and resonance curve parameters
    #B = 180 # nT, at scattering location
    n = 1 # electrons/cm^3
    omega = 0.1 # Normalized chorus gyrofrequency
    v_parallel = np.linspace(0, 1)
    harmonic = 1
    λ = 10
    L = 6

    ###### PLOTS ######

    if plotMeredithPlot:
        gs = gridspec.GridSpec(2, 2)
        psdPlt = fig.add_subplot(gs[0, 0], facecolor='k')
        alphaPlt = fig.add_subplot(gs[1, 0], facecolor='k', sharex = psdPlt)
        polarPsdPlt = fig.add_subplot(gs[:, 1], facecolor='k')
    
    psd = [np.convolve(np.ones(da)/da, np.transpose(psdObj.psd)[:, e], mode = 'same')
            for e in range(7)]
    psdPlt.plot(psdObj.times, np.transpose(psd))
    psdPlt.set(xlabel = 'UTC', ylabel = r'PSD $c^3/(cm \ MeV)^3$')
    if da > 1:
        psdPlt.text(0.2, 0.9,'{} pt smooth'.format(da), 
            horizontalalignment='center', verticalalignment='center',
            transform=psdPlt.transAxes, color = 'w')
    alphaPlt.set(xlabel = 'UTC', ylabel = 'Local Pitch angle', ylim = (0, 360))
    alphaPlt.scatter(psdObj.times, psdObj.flatAlpha)

    psdPlt.set_yscale('log')
    psdPlt.axes.get_xaxis().set_visible(False)

    # Put ticks on seconds, and rotate them.
    plt.setp(alphaPlt.xaxis.get_majorticklabels(), rotation=30, ha='right')

    for i in range(7):
        # Plot energy contours in momenum space.
        if drawEqualE:
            polarPsdPlt.plot(p_e_perp[i], p_e_parallel[i], 'w--')
            
        # Draw diffusion curves
        if drawDiffusionCurves and i % 2 == 0:
            v_perpFlt, v_parallelFlt = psdObj.calcDiffusionCurveConstantUWrapper(
                v_parallel, psdObj.Emid[i], magB(λ, L)*1E9, n, omega)
            polarPsdPlt.plot(v_perpFlt, v_parallelFlt, 'w-')
            polarPsdPlt.plot(v_perpFlt, -v_parallelFlt, 'w-')
            paramStr = ('L = {}\nmlat = {}\n|B| = {} nT \nomega = {}'
                ' \nn = {} cm^-3'.format(
                L, λ, round(magB(λ, L)*1E9), omega, n))
            polarPsdPlt.text(0.6, 0.8, paramStr, color = 'w',
                transform=polarPsdPlt.transAxes)    
                
                
    # Draw extrapolated PSD.
    if extPSD:
        # Load the fit parameters and bin the psd to the pitch angle bins
        popt = np.load('../data/quiet_fit_popt.npy')
        extPsd = np.nan*np.ones((7, len(alphaBins)), dtype = float) 
        
        for e in range(7):
            extPsd[e, :] = psd_fit.sinAlpha(alphaBins, popt[e, 0], popt[e, 1])
            
        # Draw the extrapolarted patches
        ax, p = psdObj.drawPatches(alphaBins, ax = polarPsdPlt, psd = extPsd)
        
        # Find the index where the alpha0Arr jumps (largest pitch angle sampled)
        ida = np.where(np.abs(np.convolve(alpha0Arr, [-1, 1])) > 10)[0][0]
        
        # Draw straight lines representing the trapped edges of the real data
        n_perp =  psdObj.p_perp(np.linspace(0, psdObj.Ehigh[-1]), alpha0Arr[ida-1]) 
        n_parallel =  psdObj.p_parallel(np.linspace(0, psdObj.Ehigh[-1]), 
            alpha0Arr[ida-1])
        s_perp =  psdObj.p_perp(np.linspace(0, psdObj.Ehigh[-1]), alpha0Arr[ida])
        s_parallel =  psdObj.p_parallel(np.linspace(0, psdObj.Ehigh[-1]), 
            alpha0Arr[ida])
        
        polarPsdPlt.plot(n_perp, n_parallel, 'w:')
        polarPsdPlt.plot(s_perp, s_parallel, 'w:') 
    ax, p = psdObj.drawPatches(alpha0Arr, ax = polarPsdPlt)
    cb = plt.colorbar(p, label = r'PSD $c^3/(cm \ MeV)^3$')    
    
    
    # Draw resonance curves
    if drawResonanceCurves:
        
        v_perp_res = resCurveVperp(c*v_parallel, wce(λ, L)*0.1, 1E6*n, λ, L, 
            n = harmonic)/c
#        polarPsdPlt.plot(v_perp_res, v_parallel, 'k--', 
#            v_perp_res, -v_parallel, 'k--', 
#            label = r'$\omega = 0.4$ resonance curve')
        polarPsdPlt.plot(v_perp_res, v_parallel, 'k--')
        polarPsdPlt.plot(v_perp_res, -v_parallel, 'k--', 
            label = r'$\omega = 0.1$ resonance curve')
        v_perp_res = resCurveVperp(c*v_parallel, wce(λ, L)*0.4, 1E6*n, λ, L, 
            n = harmonic)/c
        polarPsdPlt.plot(v_perp_res, v_parallel, 'k:')
        polarPsdPlt.plot(v_perp_res, -v_parallel, 'k:', 
            label = r'$\omega = 0.4$ resonance curve')
        polarPsdPlt.legend(loc = 4)
        
    polarPsdPlt.set(xlabel = r'$p_{\perp}/m_e c$', 
        ylabel = r'$p_{\parallel}/m_e c$')
        
    fig.suptitle('RBSP-{} phase space density for {} \n {} - {} UT'.format(rb_id, 
        tBounds[0].date(), tBounds[0].strftime("%H:%M:%S"), 
        tBounds[1].strftime("%H:%M:%S")))

    #ax.set_rlim(3, 1000)
    #ax.set_yscale('log')
#    ax.set(title = 'Phase space density from RBSP-{} \n Spin at {}'.format(rb_id, psdObj.magEISdata['Epoch'][0].isoformat()),
#    xlabel = r'$p_{\perp}/m_e c$', ylabel = r'$p_{\parallel}/m_e c$')
#    ax.set_xlim(right = ax.get_ylim()[1])
    #polarPsdPlt.set_aspect('equal')
    polarPsdPlt.set(ylim = (-1.2, 1.2), xlim = (0, 1.2))
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.95])
    
    if saveMeredithPlt:
        plt.savefig('rbsp{}_psd_{}_{}_to_{}{}.pdf'.format(rb_id.lower(), 
            tBounds[0].date(), tBounds[0].strftime("%H%M%S"), 
            tBounds[1].strftime("%H%M%S"), save_description))
        plt.savefig('rbsp{}_psd_{}_{}_to_{}{}.png'.format(rb_id.lower(), 
            tBounds[0].date(), tBounds[0].strftime("%H%M%S"), 
            tBounds[1].strftime("%H%M%S"), save_description))
    else:    
        plt.show()

