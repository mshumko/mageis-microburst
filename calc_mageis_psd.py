# Code to create the diffusion surfaces from the magEIS data.

### Imports ###
import numpy as np
import sys
import os
import copy
from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.lines as mlines
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import operator
from matplotlib.patches import Polygon
import matplotlib.collections

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp')
import plot_emfisis_spectra
import plot_mageis
import psd_fit
import resonant_diffusion_curves

### Constants of nature ###
c = 3E8 # m/s
mu_0 = 1.26E-6 # H·m−1 or N·A−2
eps0 = 8.85E-12 # F/m
m_p = 1.7E-27 # kg
m_e = 9.1E-31 # kg
q_e = -1.6E-19 # C
Erest = 511 # keV

class PhaseSpaceDensity(plot_mageis.PlotMageis): # Utilize inheritance
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
        plot_mageis.PlotMageis.__init__(self, rb_id, self.tRange[0], 'highrate', 
            tRange=self.tRange, instrument='low')
        return
                
    def loadData(self):
        """
        NAME:    loadMagEIS
        USE:     This function loads in magEIS data and 
                 magnetic ephemeris.
        INPUT:   REQUIRED:
                    None
                 OPTIONAL:
                    None

        AUTHOR:  Mykhaylo Shumko
        RETURNS: None
        MOD:     2017-09-07
        """
        self.times, self.j = self.getFluxTimeseries()
        self.flatAlpha = self.magEISdata['HighRate_Alpha'][:, :1000].flatten()

        # Refilter the times since initially it was calculated from the
        # unflattened arrays and now they are flattened and could be 
        # narrowed down more.
        validInd = np.where((self.times > self.tRange[0]) & 
                            (self.times < self.tRange[1]))[0]
        self.times = self.times[validInd]
        self.flatAlpha = self.flatAlpha[validInd]
        self.j = self.j[validInd, :]

        # Widen the user time range so that a) there will be at least one magEphem
        # data point and b) so that it can be flattened, without any part of the
        # user's time range cut off.      
        self.tRange = [self.tRange[0] - timedelta(minutes=2), 
            self.tRange[1] + timedelta(minutes=2)]  
        self.loadMagEphem()        
        self.BeqOverB = 1/np.mean(self.magEphem['BoverBeq'])        
        # Equatorial pitch angle array
        self.alpha0Arr = self.alpha0(self.BeqOverB, self.flatAlpha)
        return
        
    def calcPsd(self, lowerAlpha = None, upperAlpha = None):
        """
        This function calculates the phase space density, and it filters by
        pitch angle and invalid indicies
        """
        # Logic to filter calculate PSD using the correct indicies.
        if lowerAlpha is not None:
            lowerIda = np.where(self.flatAlpha >= lowerAlpha)[0]
        else:
            lowerIda = range(len(self.flatAlpha))
        if upperAlpha is not None:
            upperIda = np.where(self.flatAlpha <= upperAlpha)[0]
        else:
            upperIda = range(len(self.flatAlpha))
        validIda = np.array(list(set(lowerIda) & set(upperIda)))

        # Create data phase space density in units of c^3/(cm^6 s^3 meV^3).
        self.psd = 1E9*np.array([self.f(self.j[validIda, ee], self.Emid[ee]) for ee in range(7)]) 
        return self.psd

    def binPsdAlpha(self, binEdges, psdErr = None, zeroPsdFill = 0, 
            remapAlpha=False):
        """
        NAME:    binPsdAlpha(self, binEdges, psdErr = None, zeroPsdFill = 0, 
                    remapAlpha = True)
        USE:     This function will bin the phase space density using the 
                 pitch angle binEdges.
        INPUT:   REQUIRED:
                    binEdges - Local pitch angle bins.
                 OPTIONAL:
                    remapAlpha = True: Fold the binEdges array and the observed
                        pitch angles into 0-180 degrees.
                    psdErr = None: <<< THIS MAY BE WRONG! >>> Phase space 
                        density errors. Will claulcuate error due to Poisson
                        statistics if None. 
                    zeroPsdFill = 0: If no phase space density is found in the 
                        bin, assign it a value of 0 instead of Nan.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: self.alphaBinMid - Middle pitch angles from the binEdges,
                 self.meanPsd - (nE, len(binEdges)-1) array where nE is the 
                    number of energy channels.
        MOD:     2017-11-27
        """
        self.meanPsd = np.zeros((self.psd.shape[0], len(binEdges) - 1))
        if psdErr is not None:
            self.meanPsdErr = np.zeros((self.psd.shape[0], len(binEdges) - 1))
        self.alphaBinMid = np.convolve([0.5, 0.5], binEdges, mode = 'valid')
        
        # Loop over pitch angle bins.
        for a in range(len(binEdges[:-1])):
            validAlpha = np.where((self.flatAlpha >= binEdges[a]) & 
                                  (self.flatAlpha < binEdges[a+1]))[0] 

            # If no data points found at that equatorial pitch angle, fill in
            # the psd with a dummy variable zeroPsdFill (This avoids nan's for
            # later plotting).
            if len(validAlpha) == 0:
                self.meanPsd[:, a] = zeroPsdFill*np.ones((self.psd.shape[0]))
                continue
            # Calculate mean phase space density in that bin.
            self.meanPsd[:, a] = np.mean(self.psd[:, validAlpha], axis = 1)
        return self.alphaBinMid, self.meanPsd
        
    def drawPatches(self, eqPitchAngleBins, **kwargs):
        """
        NAME:    drawPatches(self, eqPitchAngleBins, **kwargs)
        USE:     This function takes a set of pitch angle bins and plots
                 patches in momentum space. The eqPitchAngleBins and 
                 energyBins are used to create quadralateralls and plot them.
        INPUT:   REQUIRED:
                    eqPitchAngleBins - Equatorial pitch angle bins/
                 OPTIONAL:
                    energyBins = None - Energy bins. If None will use 
                        self.Ebins
                    ax = None: Subplot object
                    psd = None: Phase space density for the bins. It None will
                        use self.meanPsd.
                    vmin = 10**-4: Sets the minimum value to display with the 
                        colormap.
                    vmax = 10**-1: Sets the maximum value to display with the 
                        colormap.
                    cmap = 'plasma': Sets the colormap for the patches.
                    cMapLog = True: Sets the colorbar scale to log or linear 
                        scale.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: ax - subplot object, and p - patches objects.
        MOD:     2017-10-03
        """
        energyBins = kwargs.get('energyBins', None)
        ax = kwargs.get('ax', None)
        psd = kwargs.get('psd', None)
        vmin = kwargs.get('vmin', 10**-4)
        vmax = kwargs.get('vmax', 10**-1)
        cmap = kwargs.get('cmap', 'plasma')
        cMapLog = kwargs.get('cMapLog', True)
        
        # Create the verticies from the pitch angle bins.
        self.__makePatchVerticies__(eqPitchAngleBins, energyBins, psd = psd)   
        # Use lists here to filter out the pitch angle bins with no samples.   
        patches = []
        c = []

        for i in range(self.verticies.shape[2]):
            c.append(self.c[i])
            patches.append(Polygon(self.verticies[:, :, i]))

        p = matplotlib.collections.PatchCollection(patches)
        p.set_cmap(cmap)
        p.set_array(np.array(c))
        p.autoscale()
        if cMapLog:
            p.set_norm(matplotlib.colors.LogNorm())
        p.set_clim(vmin=vmin, vmax=vmax)
        ax.add_collection(p)
        return ax, p
        
    def __makePatchVerticies__(self, eqPitchAngleBins, energyBins = None,
            ignoreBoxHeght = 0.2, psd = None):
        """
        NAME:    __makePatchVerticies__(self, eqPitchAngleBins, 
                    energyBins = None, ignoreBoxHeght = 0.2, psd = None)
        USE:     This helper function for drawPatches calculates the locations
                 of the patch verticies from pitch angle and energy bins (It 
                 maps them to momentum space).
        INPUT:   REQUIRED:
                    eqPitchAngleBins - Equatorial pitch angle bins/
                 OPTIONAL:
                    energyBins = None - Energy bins. If None will use 
                        self.Ebins
                    psd = None: Phase space density for the bins. It None will
                        use self.meanPsd.
                    ignoreBoxHeght = 0.2: If there is a gap in pitch angles 
                        (due to partial equatorial picth angle coverage) do not
                        try to draw boxes over the non-sampled pitch angles. 
                        This compares the distance between the lower left to
                        upper right verticies.
        AUTHOR:  Mykhaylo Shumko
        RETURNS: self.verticies a (4, 2, nV) array where the first dimention is
                 over the verticies, the second dimention is the p_perp, 
                 p_parallel coordinates of those verticies, and the third 
                 dimention of size nV are the patches.
        MOD:     2017-10-03
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
    def f(j, Ek, Erest=Erest, jErr=None):
        psd = np.divide(j, 100*c*Ek*(Ek + 2*Erest))
        if jErr is None:
            return psd
        else:
            psdErr = np.divide(jErr, 100*c*Ek*(Ek + 2*Erest))
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
    p = staticmethod(lambda T, Erest = Erest: np.sqrt(T**2 + 2*T*Erest)/(100*c))

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
    tBoundsDict = {'q':[datetime(2017, 3, 31, 11, 15, 0), datetime(2017, 3, 31, 11, 17, 0)], 
    'm1':[datetime(2017, 3, 31, 11, 17, 0), datetime(2017, 3, 31, 11, 17, 20)], 
    'm2':[datetime(2017, 3, 31, 11, 17), datetime(2017, 3, 31, 11, 17, 11)],
    'bigOne':[datetime(2017, 3, 31, 11, 17, 13), datetime(2017, 3, 31, 11, 17, 18)],
    'smallOne':[datetime(2017, 3, 31, 11, 17, 9, 500000), datetime(2017, 3, 31, 11, 17, 10, 500000)]}
    tPeriod = 'm2'

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
    psdObj.binPsdAlpha(alphaBins)
    plotPitchAngles = False
    plotMeredithPlot = True
    saveMeredithPlt = False
    drawDiffusionCurves = True
    drawResonanceCurves = True
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
    alphaPlt.set(xlabel = 'UTC', ylabel = 'Local Pitch angle', ylim = (0, 180))
    alphaPlt.scatter(psdObj.times, psdObj.flatAlpha)

    psdPlt.set_yscale('log')
    psdPlt.axes.get_xaxis().set_visible(False)

    # Put ticks on seconds, and rotate them.
    plt.setp(alphaPlt.xaxis.get_majorticklabels(), rotation=30, ha='right')  
                
    # Draw extrapolated PSD.
    if extPSD:
        # Load the fit parameters and bin the psd to the pitch angle bins
        dir_path = os.path.dirname(os.path.realpath(__file__)) # Script dir
        popt = np.load(os.path.join(dir_path, 'data', 'quiet_fit_popt.npy'))
        extPsd = np.nan*np.ones((7, len(alphaBins)), dtype = float) 
        
        for e in range(7):
            extPsd[e, :] = psd_fit.sinAlpha(alphaBins, popt[e, 0], popt[e, 1])
            
        # Draw the extrapolarted patches
        ax, p = psdObj.drawPatches(alphaBins, ax=polarPsdPlt, psd=extPsd,
            cMapLog=False)
        
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
    ax, p = psdObj.drawPatches(alpha0Arr, ax=polarPsdPlt, cMapLog=False)
    cb = plt.colorbar(p, label = r'PSD $c^3/(cm \ MeV)^3$')    
    

    # Plot energy contours in momenum space.
    if drawEqualE:
        for i in range(7):
            polarPsdPlt.plot(p_e_perp[i], p_e_parallel[i], 'w--')

    # Resonant-diffusion parameters
    vParallel_res = c*np.linspace(0, -0.99, num = 1000)
    mlat = 0
    mlat0 = 20
    L = 5.7
    n0 = 0.5E6 # Density at the time

    # Draw resonance curves
    if drawResonanceCurves:
        # w/w_ce = 0.1 
        vPerp_res = resonant_diffusion_curves.resCurveVperp(
            vParallel_res, 0.1*resonant_diffusion_curves.wce(mlat, L), n0,
            mlat0, mlat, L)
        pPerp_res, pParallel_res = resonant_diffusion_curves.p(
            vPerp_res, vParallel_res)
        polarPsdPlt.plot(pPerp_res, pParallel_res, 'g')
        polarPsdPlt.plot(pPerp_res, -pParallel_res, 'g')
        label01 = mlines.Line2D([], [], color='g', markersize=15, 
            label=r'$0.1 \ \Omega_{ce}$')

        # w/w_ce = 0.4 
        vPerp_res = resonant_diffusion_curves.resCurveVperp(
            vParallel_res, 0.4*resonant_diffusion_curves.wce(mlat, L), n0, 
            mlat0, mlat, L)
        pPerp_res, pParallel_res = resonant_diffusion_curves.p(
            vPerp_res, vParallel_res)
        polarPsdPlt.plot(pPerp_res, pParallel_res, 'r')
        polarPsdPlt.plot(pPerp_res, -pParallel_res, 'r')
        label04 = mlines.Line2D([], [], color='r', markersize=15, ls='-',
            label=r'$0.4 \ \Omega_{ce}$')

        # w/w_ce = 0.6
        vPerp_res = resonant_diffusion_curves.resCurveVperp(
            vParallel_res, 0.6*resonant_diffusion_curves.wce(mlat, L), n0, 
            mlat0, mlat, L)
        pPerp_res, pParallel_res = resonant_diffusion_curves.p(
            vPerp_res, vParallel_res)
        polarPsdPlt.plot(pPerp_res, pParallel_res, 'b')
        polarPsdPlt.plot(pPerp_res, -pParallel_res, 'b')
        label06 = mlines.Line2D([], [], color='b', markersize=15, ls='-',
            label=r'$0.6 \ \Omega_{ce}$')
            
        polarPsdPlt.legend(loc=4, handles=[label01, label04, label06], fontsize=10)
        
    if drawDiffusionCurves:
        Earr = psdObj.Emid
        vParallel_diff = vParallel_res
        
        for e in Earr[::2]:
            vPerp_diff = resonant_diffusion_curves.diffCurveVperp(
                vParallel_diff, 0.4*resonant_diffusion_curves.wce(mlat, L), 
                n0, mlat0, mlat, L, e)
            pPerp_diff, pParallel_diff = resonant_diffusion_curves.p(vPerp_diff, 
                vParallel_diff)
            polarPsdPlt.plot(pPerp_diff, pParallel_diff, 'c--')
            polarPsdPlt.plot(pPerp_diff, -pParallel_diff, 'c--')
            
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
    polarPsdPlt.set_aspect('equal')
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

