"""
This script and supporting functions are made to produce figure 2 in the
MagEIS microburst paper

Mykhaylo Shumko
Last modified: 2017-10-18
"""

import numpy as np
import sys
import os
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')))
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#import matplotlib.dates as mdates
from datetime import datetime, timedelta

#import operator

sys.path.insert(0, '/home/mike/Dropbox/0_grad_work/mission_tools/rbsp')
import plot_mageis_spectra
import plot_rbspice

# "Interactive time range selection."
tKey = 'muBurst'
times = {'muBurst':[datetime(2017, 3, 31, 11, 15, 0), 
                    datetime(2017, 3, 31, 11, 18, 10)],
            'later':[datetime(2017, 3, 31, 11, 35, 0), 
                    datetime(2017, 3, 31, 11, 38)],
            'all':[datetime(2017, 3, 31, 11, 15), 
                    datetime(2017, 3, 31, 11, 20)]}
tRange = times[tKey]
