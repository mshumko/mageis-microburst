# This script will plot the HOPE data and convert it to PSD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import statistics
from datetime import datetime

from spacepy import pycdf
# import PSD library
sys.path.insert(0, '/home/mike/research/mageis-microburst')
import calc_mageis_psd

START_TIME = datetime(2017, 3, 31, 11, 16, 0)
END_TIME = datetime(2017, 3, 31, 11, 18, 0)
START_ENERGY_IND = 50
END_ENERGY_IND = 70
PA_PLOT_IND = 2

# Load in data 
# Usefull keys: PITCH_ANGLE [11], Epoch_Ele [nT], FEDU [nT, 11, 72], 
# FEDO [nT, 72], HOPE_ENERGY_Ele [nT, 72]
dPath = ('/home/mike/research/rbsp/data/hope/rbspa/'
        'rbspa_rel03_ect-hope-PA-L3_20170331_v6.1.0.cdf')
d = pycdf.CDF(dPath).copy()
# Filter times and a few usefull variables
validIdt = np.where((d['Epoch_Ele'] > START_TIME)
            & (d['Epoch_Ele'] < END_TIME))[0]
d['Epoch_Ele'] = d['Epoch_Ele'][validIdt]
d['FEDU'] = d['FEDU'][validIdt, :, :]

d['FEDO'] = d['FEDO'][validIdt, :]
d['HOPE_ENERGY_Ele'] = d['HOPE_ENERGY_Ele'][validIdt, :]
d['e_energy_ch'] = [statistics.mode(i) for i in np.transpose(d['HOPE_ENERGY_Ele'])]
# Calculate PSD

for a in range(11):
    d['FEDU_PSD'] = np.nan*np.ones_like(d['HOPE_ENERGY_Ele'])
    for i, (j, Ek) in enumerate(zip(np.transpose(d['FEDU'][:, a, :]), d['e_energy_ch'])):
        d['FEDU_PSD'][:, i] = 1E9*calc_mageis_psd.PhaseSpaceDensity.f(j, Ek/1000)

    # Plot the PSD
    fig, ax = plt.subplots()
    ax.axhline(0.1)
    colors = cm.rainbow(np.linspace(0, 1, END_ENERGY_IND - START_ENERGY_IND))
    for (i, ii) in enumerate(range(START_ENERGY_IND, END_ENERGY_IND)):
        ax.scatter(d['Epoch_Ele'], d['FEDU_PSD'][:, ii], 
                    label='{}'.format(round(d['e_energy_ch'][ii]/1000)), c=colors[i])
    ax.legend(loc=1, bbox_to_anchor=(1.15, 1.1), title='Energy [keV]',
                fontsize=10)
    
    ax.set_yscale('log')
    ax.set_xlabel('UTC')
    ax.set_ylabel(r'PSD $c^3/(cm \ MeV)^3$')
    ax.set_xlim(START_TIME, END_TIME)
    ax.set_ylim(bottom=10**-2, top=10)
    fig.autofmt_xdate()
    
    ax.set_title(r'2017-03-31 | RBSP-A | HOPE Electrons | $\alpha_L$' + 
                ' = {}'.format(d['PITCH_ANGLE'][a]))
    plt.savefig(('/home/mike/Dropbox/0_grad_work/mageis_microburst/plots/'
                'HOPE/20170331_hope_psd_a_{}.png'.format(int(d['PITCH_ANGLE'][a]))))
    plt.close()