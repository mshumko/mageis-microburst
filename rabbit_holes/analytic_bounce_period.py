# Calculate electron bounce periods at VAP
import sys
sys.path.insert(0, '/home/ms30715/research_code/auxiliary')
import numpy as np
import plot_mageis_spectra   
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from datetime import datetime, timedelta

Re = 6371 #km
c = 3.0E8 # m/s

### Relativistic functions for velocity and dialation effect ###
beta = lambda Ek, Erest = 511: np.sqrt(1-((Ek/Erest)+1)**(-2)) # Ek,a dErest must be in keV
gamma = lambda Ek, Erest = 511:np.sqrt(1-beta(Ek, Erest = 511)**2)**(-1/2)
Tsl = lambda L, alpha0, v: 4*6.371E6*np.divide(L, v) * \
       (1.3802 - 0.3198*(np.sin(np.deg2rad(alpha0)) + \
       np.sqrt(np.sin(np.deg2rad(alpha0)))))
       
def alpha0(*args):
    """
    First argument is B_eq, second is B_sc and third is alpha_sc. Or can provide
    B_eq/B_sc directly and alpha_sc
    """
    assert len(args) in [3, 2], "Number of arguments must be 2 or 3!"
    
    if len(args) == 3:
        a0 = np.arcsin(np.sqrt(np.divide(args[0], args[1]))*np.sin(np.deg2rad(args[2])))
    elif len(args) == 2:
        a0 = np.arcsin(np.sqrt(args[0])*np.sin(np.deg2rad(args[1])))
    return np.rad2deg(a0)
    
#print(alpha0(0.5, 40))

### Script starts here ###
if __name__ == '__main__':
    rb_id = 'A'
    date = datetime(2017, 3, 31)
    tBounds = [datetime(2017, 3, 31, 11, 15), datetime(2017, 3, 31, 11, 25)]
    fluxObj = plot_mageis_spectra.magEISspectra(rb_id, date, dataLevel = 3)
    fluxObj.tBounds = tBounds
    fluxObj.loadMagEphem()
    
    BeqOverB = 1/np.mean(fluxObj.magEphem['BoverBeq']) # Asumme its fairly constant
    L = np.mean(fluxObj.magEphem['Lstar'])
    
    # Create a meshgrid of pitch angles and energies
    alphaArr = np.arange(0, 180, 1)
    alpha0Arr = alpha0(BeqOverB, alphaArr)
    Earr = np.arange(20, 100, 1)

    aa, ee = np.meshgrid(alphaArr, Earr)
    a0a0, ee = np.meshgrid(alpha0Arr, Earr)
    tb = Tsl(L, a0a0, c*beta(ee))
    
    cflevels = np.arange(np.min(tb), np.max(tb), 0.01)
    clevels = np.arange(0.5, 3, 0.25)
    
    cf = plt.contourf(aa, ee, tb, c = tb, levels = cflevels, s = 40)#, 
    #    vmin=clevels[0], vmax=clevels[-1])    
    c = plt.contour(aa, ee, tb, colors = 'k', lw = 5, levels = clevels)
    plt.clabel(c, inline=True, colors = 'k' ,
                  fontsize=10, manual = False, fmt = '%.2f')
    
    plt.colorbar(cf, label = r'$T_b$')
    plt.title('Schulz & Lanzerotti electron bounce period at RBSP-A \n Time: {}'.format(tBounds[0].isoformat()))
    plt.xlabel(r'$\alpha_{rbsp}$ (local pitch angle)')
    plt.ylabel('Energy (keV)')
    plt.tight_layout()
    plt.show()
