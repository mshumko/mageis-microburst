# This script calculates the return time.
import numpy as np
from datetime import datetime, timedelta
#import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib.gridspec as gridspec
plt.rcParams.update({'font.size': 12})

import sys
sys.path.insert(0, '/home/mike/research/mission-tools/vap')
import plot_mageis
import IRBEM

# FIREBIRD location
#X = {'x1':651, 'x2':63, 'x3':15.9, 'dateTime':datetime(2015, 2, 2, 6, 12, 43)}
# RBSP-A location at the time of the microburst
X = {'x1':30409, 'x2':-8.8, 'x3':115.1, 'dateTime':datetime(2017, 3, 31, 11, 18)}
maginput = {'Kp':40}

model = IRBEM.MagFields()

alphaArr = np.arange(10, 90, 5)
Earr = np.arange(30, 100, 5)
aa, ee = np.meshgrid(alphaArr, Earr) # Meshgrid for plotting
tArr = np.nan*np.ones((len(Earr), len(alphaArr))) # Bounce period array.

for i_a, a in enumerate(alphaArr):
    tArr[:, i_a] = model.bounce_period(X, maginput, Earr, alpha = a)
    
# Plot the bounce periods
fig = plt.figure(figsize=(15, 8), dpi=80, facecolor='white')
gs = gridspec.GridSpec(2,2)
ax = fig.add_subplot(gs[:, 0], facecolor='k')
p = ax.pcolormesh(ee, aa, tArr)
CS = ax.contour(ee, aa, tArr, 6, colors='w')
ax.clabel(CS, inline=1)
plt.colorbar(p, label='Bounce period (s)')
ax.set(xlabel='Energy (keV)', ylabel=r'$\alpha_{sc}$', 
    title='Electron bounce period in the T89 model, Kp = {}'.format(
    maginput['Kp']/10))

# Draw boxes in E-alpha space for the microbursts

### PLOT the HighRate data ###
#fig2 = plt.figure(figsize=(15, 10), dpi=80, facecolor='white')
#gs = gridspec.GridSpec(2,1)
bx = fig.add_subplot(gs[0, 1])
cx = fig.add_subplot(gs[1, 1], sharex=bx)
#dx = fig2.add_subplot(gs[3, 0], sharex=bx)

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

### NOW PLOT PITCH ANGLES ###
validIda = np.where(fluxObj.magEISdata[fluxObj.alphaKey] != -1E31)
alpha = fluxObj.magEISdata[fluxObj.alphaKey][validIda]

# Now map alpha into a 180 degree space
flipAlpha = np.where(alpha > 180)[0]
alpha[flipAlpha] = -alpha[flipAlpha] % 360

cx.fill_between(fluxObj.times[validIda], alpha + 10, alpha - 10, color='b')
cx.plot(fluxObj.times[validIda], alpha, 'r')
    
cx.set(ylim=(0, 60), xlim=(datetime(2017, 3, 31, 11, 17, 8), 
                    datetime(2017, 3, 31, 11, 17, 13)) )
# Draw guides
t = [None]*4
t[0] = datetime(2017, 3, 31, 11, 17, 9, int(78E4))
t[1] = datetime(2017, 3, 31, 11, 17, 10, int(24E4))
t[2] = datetime(2017, 3, 31, 11, 17, 10, int(85E4))
t[3] = datetime(2017, 3, 31, 11, 17, 11, int(48E4))

tt = [None]*4
for i in range(4):
    tt[i] = md.date2num(t[i])

for a in [bx, cx]:    
    a.axvline(x=t[0], c='k')
    a.axvline(x=t[1], c='k')
    a.axvline(x=t[2], c='k')
    a.axvline(x=t[3], c='k')

plt.setp(bx.get_xticklabels(), visible=False)

bx.set_ylim(bottom=5*10**2)
bx.set(title='MagEIS LOW timeseries from {}'.format(date.isoformat()))
bx.legend()
bx.set_ylabel(r'Flux $(keV \ cm^2 \ sr \ s)^{-1}$')
cx.set_title('MagEIS LOW pitch angle evolution')
cx.set_ylabel(r'$\alpha_{eq}$')

#dx.set_ylabel(r'$\alpha_{eq}$')
cx.set_xlabel('UTC')
gs.tight_layout(fig)
plt.show()
