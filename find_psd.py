import numpy as np
import mageis_diffusion_curves
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime

def sinAlpha(alpha, A, n):
    """
    This function will return a value from the function A*sin(alpha)^n. 
    This is used for fitting the equatorial pitch angle distribution.
    """
    return A*np.sin(np.deg2rad(alpha))**n


tBoundsDict = {'q':[datetime(2017, 3, 31, 11, 15, 0), datetime(2017, 3, 31, 11, 17, 0)], 'm1':[datetime(2017, 3, 31, 11, 17, 0), datetime(2017, 3, 31, 11, 17, 20)], 'm2':[datetime(2017, 3, 31, 11, 17, 10), datetime(2017, 3, 31, 11, 17, 20)],
    'bigOne':[datetime(2017, 3, 31, 11, 17, 13), datetime(2017, 3, 31, 11, 17, 18)],
    'smallOne':[datetime(2017, 3, 31, 11, 17, 9, 500000), datetime(2017, 3, 31, 11, 17, 10, 500000)], 'None':None}
    tBounds = tBoundsDict['None']
