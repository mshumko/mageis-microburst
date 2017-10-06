# This script calculates the return time.
import numpy as np
from datetime import datetime
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.colors

import IRBEM

R_e = 6.371E6 # meters

beta = lambda Ek, Erest = 511: np.sqrt(1-((Ek/Erest)+1)**(-2)) # Ek,a dErest must
vparalel = lambda Ek, Bm, B, alpha:c*beta(Ek)*np.sqrt(
    1 - np.sin(np.deg2rad(alpha))*np.abs(B/Bm) )
    
    
def tReturn(X, maginput, alpha, hemi='same'):
    """
    
    """
    model = IRBEM.MagFields()
    # Trace field line
    fLine = model._interpolate_field_line(X, maginput)
    
    B_m = fLine['inputB']/np.sin(np.deg2rad(alpha))**2 # Mirror point |B|
    
    # Tranformed coordinates of the sc location.
    coords = IRBEM.Coords()
    scLoc = coords.coords_transform(X['dateTime'], 
        [X['x1'], X['x2'], X['x3']], 0, 1)
    
    if scLoc[0, 2] > 0:
        # Find index (and location of spacecraft)
        S_sc = scipy.optimize.brentq(fLine['fB'], 0, len(fLine['S'])/2)
        
        # Find the index of the mirror point.
        if hemi == 'same':
            S_m = scipy.optimize.brentq(fLine['fB']+fLine['inputB']-B_m, 0, 
                len(fLine['S'])/2)
    else:
        S_sc = scipy.optimize.brentq(fLine['fB'], len(fLine['S'])/2, 
                                             len(fLine['S'])-1)
        # Find the index of the mirror point.                                    
        if hemi == 'same':
            S_m = scipy.optimize.brentq(fLine['fB']+fLine['inputB']-B_m, 
                len(fLine['S'])/2, len(fLine['S'])-1)
                
    #plt.plot(fLine['fx'](fLine['S']), fLine['fz'](fLine['S']))
    plt.scatter(fLine['fx'](fLine['S']), fLine['fz'](fLine['S']), 
        c=fLine['fB'](fLine['S'])+fLine['inputB'], norm=matplotlib.colors.LogNorm())
    #plt.scatter(scLoc[0, 0], scLoc[0, -1], c='r', marker='*', s=500)
    plt.scatter(fLine['fx'](S_sc), fLine['fz'](S_sc), c='r', marker='*', s=500)
    plt.scatter(fLine['fx'](S_m), fLine['fz'](S_m), c='b', marker='*', s=500)
    plt.show()
    return 

# FIREBIRD location
#X = {'x1':651, 'x2':63, 'x3':15.9, 'dateTime':datetime(2015, 2, 2, 6, 12, 43)}
# RBSP-A location at the time of the microburst
X = {'x1':30409, 'x2':-8.8, 'x3':115.1, 'dateTime':datetime(2017, 3, 31, 11, 18)}
maginput = {'Kp':20}
tReturn(X, maginput, 320)

###model = IRBEM.MagFields()
###fLine = model._interpolate_field_line(X, maginput)
###B_m = fLine['inputB']/np.sin(a)**2

###t = [2*np.sum(np.divide(ds[1:-1], vparalel(Ei, fLine['inputB'], dB, 
###                                              Erest = Erest)[1:-1])) for Ei in E]

