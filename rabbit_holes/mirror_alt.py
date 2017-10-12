# This script will calculate the mirror point altitude for 
# the MagEIS microburst electrons
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

import IRBEM

Re = 6371 #km

X = {'x1':30409, 'x2':-8.8, 'x3':115.1, 'dateTime':datetime(2017, 3, 31, 11, 18)}
maginput = {'Kp':40}
alpha = 40

model = IRBEM.MagFields()
mirror_dict = model.find_mirror_point(X, maginput, alpha)
alt = Re*(np.sqrt(mirror_dict['POSIT'][0]**2+mirror_dict['POSIT'][1]**2+mirror_dict['POSIT'][2]**2) - 1)
print(alt)
