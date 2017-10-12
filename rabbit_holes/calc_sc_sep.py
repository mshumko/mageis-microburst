import spacepy.datamodel 
from datetime import datetime
import dateutil.parser
import numpy as np

# Magnetic epehemeris paths
ephemPathA = '/home/mike/research/rbsp/magephem/rbspa/rbspa_def_MagEphem_TS04D_20170331_v1.0.0.txt'
ephemPathB = '/home/mike/research/rbsp/magephem/rbspb/rbspb_def_MagEphem_TS04D_20170331_v1.0.0.txt'

Re = 6371 # km

# Times of interest
tRange = [datetime(2017, 3, 31, 11, 10), datetime(2017, 3, 31, 11, 20)]

mEphemA = spacepy.datamodel.readJSONheadedASCII(ephemPathA)
mEphemB = spacepy.datamodel.readJSONheadedASCII(ephemPathB)

# Parse and filter by times
timeA = np.array(list(map(lambda t: dateutil.parser.parse(t, ignoretz = True), 
    mEphemA['DateTime'])))
timeB = np.array(list(map(lambda t: dateutil.parser.parse(t, ignoretz = True),
    mEphemB['DateTime'])))
idtA = np.where((timeA > tRange[0]) & (timeA < tRange[1]))[0]
idtB = np.where((timeB > tRange[0]) & (timeB < tRange[1]))[0]

# The GEO positions
posA = mEphemA['Rgeo'][idtA]
posB = mEphemB['Rgeo'][idtB]

d = Re*np.sqrt((posA[:, 0] - posB[:, 0])**2 + (posA[:, 1] - posB[:, 1])**2 + 
    (posA[:, 2] - posB[:, 2])**2)

print('Time               Separation (km)')
for i, t in enumerate(timeA[idtA]):
    print(t, d[i])

