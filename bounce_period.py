# This script calculates the return time.
import numpy as np
from datetime import datetime
#import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0, '/home/mike/research/mission-tools/vap')
import plot_mageis
import IRBEM

# FIREBIRD location
#X = {'x1':651, 'x2':63, 'x3':15.9, 'dateTime':datetime(2015, 2, 2, 6, 12, 43)}
# RBSP-A location at the time of the microburst
X = {'x1':30409, 'x2':-8.8, 'x3':115.1, 'dateTime':datetime(2017, 3, 31, 11, 18)}
maginput = {'Kp':40}
#tReturn(X, maginput, 320)

model = IRBEM.MagFields()

alphaArr = np.arange(10, 90, 5)
Earr = np.arange(30, 200, 5)
aa, ee = np.meshgrid(alphaArr, Earr) # Meshgrid for plotting
tArr = np.nan*np.ones((len(Earr), len(alphaArr))) # Bounce period array.

for i_a, a in enumerate(alphaArr):
    tArr[:, i_a] = model.bounce_period(X, maginput, Earr, alpha = a)
    
# Plot the bounce periods
fig = plt.figure(figsize=(15, 10), dpi=80, facecolor='white')
gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0, 0], facecolor='k')
p = ax.pcolormesh(ee, aa, tArr)
plt.colorbar(p, label='Bounce period (s)')
ax.set(xlabel='Energy (keV)', ylabel=r'$\alpha_{sc}$', 
    title='Electron bounce period in the T89 model, Kp = {}'.format(
    maginput['Kp']/10))
#plt.show()

### PLOT the HighRate data ###
fig2 = plt.figure(figsize=(15, 10), dpi=80, facecolor='white')
gs = gridspec.GridSpec(2,1)
bx = fig2.add_subplot(gs[0, 0])
cx = fig2.add_subplot(gs[1, 0], sharex=bx)

tKey = 'muBurst'
date = datetime(2017, 3, 31)
times = {'muBurst':[datetime(2017, 3, 31, 11, 17), 
                    datetime(2017, 3, 31, 11, 17, 30)],
           'quiet0':[datetime(2017, 3, 31, 2), 
                    datetime(2017, 3, 31, 2, 56)],
           'quiet1':[datetime(2017, 3, 31, 11),
                    datetime(2017, 3, 31, 11, 17)],
           'quiet2':[datetime(2017, 3, 31, 19, 45), 
                    datetime(2017, 3, 31, 20, 54)]            
        }
rb_id = 'A'
fluxObj = plot_mageis.magEISspectra(rb_id, date, dataLevel = 3)
fluxObj.tBounds = times[tKey]
fluxObj.loadMagEIS(instrument='LOW', highrate=True)
fluxObj.plotHighRateTimeSeries(smooth=10, ax=bx)
bx.set_ylim(bottom=5*10**2)
bx.set(title='MagEIS timeseries from {}'.format(date.isoformat()))

print(fluxObj.magEISdata[fluxObj.alphaKey].shape)
print(fluxObj.times.shape)
validIda = np.where(fluxObj.magEISdata[fluxObj.alphaKey] != -1E31)
cx.fill_between(fluxObj.times[validIda], 
    fluxObj.magEISdata[fluxObj.alphaKey][validIda] + 10,
    fluxObj.magEISdata[fluxObj.alphaKey][validIda] - 10)
cx.plot(fluxObj.times[validIda], fluxObj.magEISdata[fluxObj.alphaKey][validIda], 
    'r')

plt.show()


# Now calculate the 

###R_e = 6.371E6 # meters

###beta = lambda Ek, Erest = 511: np.sqrt(1-((Ek/Erest)+1)**(-2)) # Ek,a dErest must
###vparalel = lambda Ek, Bm, B, alpha:c*beta(Ek)*np.sqrt(
###    1 - np.sin(np.deg2rad(alpha))*np.abs(B/Bm) )
    
    
###def tReturn(X, maginput, alpha, hemi='same'):
###    """
###    
###    """
###    model = IRBEM.MagFields()
###    # Trace field line
###    fLine = model._interpolate_field_line(X, maginput)
###    
###    B_m = fLine['inputB']/np.sin(np.deg2rad(alpha))**2 # Mirror point |B|
###    
###    # Tranformed coordinates of the sc location.
###    coords = IRBEM.Coords()
###    scLoc = coords.coords_transform(X['dateTime'], 
###        [X['x1'], X['x2'], X['x3']], 0, 1)
###    
###    if scLoc[0, 2] > 0:
###        # Find index (and location of spacecraft)
###        S_sc = scipy.optimize.brentq(fLine['fB'], 0, len(fLine['S'])/2)
###        
###        # Find the index of the mirror point.
###        if hemi == 'same':
###            S_m = scipy.optimize.brentq(fLine['fB']+fLine['inputB']-B_m, 0, 
###                len(fLine['S'])/2)
###    else:
###        S_sc = scipy.optimize.brentq(fLine['fB'], len(fLine['S'])/2, 
###                                             len(fLine['S'])-1)
###        # Find the index of the mirror point.                                    
###        if hemi == 'same':
###            S_m = scipy.optimize.brentq(fLine['fB']+fLine['inputB']-B_m, 
###                len(fLine['S'])/2, len(fLine['S'])-1)
###                
###    #plt.plot(fLine['fx'](fLine['S']), fLine['fz'](fLine['S']))
###    plt.scatter(fLine['fx'](fLine['S']), fLine['fz'](fLine['S']), 
###        c=fLine['fB'](fLine['S'])+fLine['inputB'], norm=matplotlib.colors.LogNorm())
###    #plt.scatter(scLoc[0, 0], scLoc[0, -1], c='r', marker='*', s=500)
###    plt.scatter(fLine['fx'](S_sc), fLine['fz'](S_sc), c='r', marker='*', s=500)
###    plt.scatter(fLine['fx'](S_m), fLine['fz'](S_m), c='b', marker='*', s=500)
###    plt.show()
###    return 

###model = IRBEM.MagFields()
###fLine = model._interpolate_field_line(X, maginput)
###B_m = fLine['inputB']/np.sin(a)**2

###t = [2*np.sum(np.divide(ds[1:-1], vparalel(Ei, fLine['inputB'], dB, 
###                                              Erest = Erest)[1:-1])) for Ei in E]

