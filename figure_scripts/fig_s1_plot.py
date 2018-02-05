# This script generates the Figure 2 plot with AC-6 data.
from datetime import datetime, timedelta
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/mike/research/mission-tools/ac6')
import read_ac_data

date = datetime(2017, 3, 31)
tRange = [datetime(2017, 3, 31, 11, 18), datetime(2017, 3, 31, 11, 20, 21)]
dataA = read_ac_data.read_ac_data_wrapper('a', date, dType = '10Hz', tRange = tRange)
dataB = read_ac_data.read_ac_data_wrapper('b', date, dType = '10Hz', tRange = tRange)
    
tLag = 10.4 # s
dLag = 78.1 # km

fig = plt.figure(figsize=(15, 10), dpi = 80, facecolor = 'white')
plt.rcParams.update({'font.size': 15})
gs = matplotlib.gridspec.GridSpec(2, 1)
sTplt = fig.add_subplot(gs[0, 0], facecolor='w')
sPplt = fig.add_subplot(gs[1, 0], facecolor='w') #, sharex = sTplt
posPlt = sPplt.twinx()

# Plot the unshifted dos rate data.
validIda = np.where(dataA['dos1rate'] != -1E31)[0]
validIdb = np.where(dataB['dos1rate'] != -1E31)[0]
sTplt.plot(dataA['dateTime'][validIda], dataA['dos1rate'][validIda], 'r',
    label = 'AC6-A dos1rate')
sTplt.plot(dataB['dateTime'][validIdb], dataB['dos1rate'][validIdb], 'b',
    label = 'AC6-B dos1rate')
###sTplt.plot(dataA['dateTime'][validIda], dataA['dos2rate'][validIda], 'r:',
###    label = 'AC6-A dos2rate')
###sTplt.plot(dataB['dateTime'][validIdb], dataB['dos2rate'][validIdb], 'b:',
###    label = 'AC6-B dos2rate')
###sTplt.plot(dataA['dateTime'][validIda], dataA['dos3rate'][validIda], 'r--',
###    label = 'AC6-A dos2rate')
###sTplt.plot(dataB['dateTime'][validIdb], dataB['dos3rate'][validIdb], 'b--',
###    label = 'AC6-B dos2rate')
    
# Plot the time shifted dos rates.
shiftedTimes = [t + timedelta(seconds = tLag) for t in dataA['dateTime'][validIda]]

sPplt.plot(shiftedTimes, dataA['dos1rate'][validIda], 'r',
    label = 'AC6-A dos1rate')
sPplt.plot(dataB['dateTime'][validIdb], dataB['dos1rate'][validIdb], 'b',
    label = 'AC6-B dos1rate')
    
# Plot spacecraft positions
validPos = np.where(dataB['Lm_OPQ'] != -1E31)[0]
posPlt.plot(dataB['dateTime'][validPos], dataB['Lm_OPQ'][validPos], 'k', 
    label = 'AC6-B Lm_OPQ')
    
# Set all of the subplot parameters here
sTplt.set(yscale = 'log', ylabel = 'counts/s', 
    xlim = (datetime(2017, 3, 31, 11, 19), datetime(2017, 3, 31, 11, 20, 20)) )
sPplt.set(yscale = 'log', ylabel = 'counts/s', xlabel = 'UTC', 
    xlim = (datetime(2017, 3, 31, 11, 19), datetime(2017, 3, 31, 11, 20, 20)) )
posPlt.set(ylim = (4, 7.5), ylabel = 'McIlwain L (black)')

sTplt.text(0.05, 0.9,'MLT = {}'.format(round(np.mean(dataA['MLT_OPQ']))), 
    horizontalalignment='left', verticalalignment='top',
    transform = sTplt.transAxes, color = 'k')

sPplt.text(0.05, 0.9,'AC6-A shifted by {} s\nin-track separation {} km'.format(tLag, dLag), 
    horizontalalignment='left', verticalalignment='top',
    transform = sPplt.transAxes, color = 'k')

abcLabels = ['(a)', '(b)']
abcColors = ['k', 'k']
for i, a in enumerate([sTplt, sPplt]):
    a.xaxis.set_minor_locator(matplotlib.dates.SecondLocator())
    a.xaxis.set_tick_params(which = 'minor', width = 2, length = 5)
    a.xaxis.set_tick_params(which = 'major', width = 2, length = 15)
    a.text(0.05, 0.95, abcLabels[i], transform=a.transAxes, va='top', 
            color=abcColors[i])     

plt.setp(sPplt.xaxis.get_majorticklabels(), rotation=30, ha='right')
plt.setp(sTplt.get_xticklabels(), visible=False)

sTplt.legend(loc = 1)
posPlt.set_ylim(top=8)
gs.tight_layout(fig)
plt.savefig('fig2.pdf')
plt.show()
