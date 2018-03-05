# This script will 

from spacepy import pycdf
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

class Proc_L0_RBSPICE:
    def __init__(self, fPath, tRange=None):
        """
        This class will load in the L0 RBSPICE data. It has methods to
        convert the count rates to flux, and spin sector to pitch angle.
        """
        self.rawD = pycdf.CDF(fPath)
        #print(self.d)
        self.tRange = tRange

        if self.tRange is not None:
            self.d = {}
            self._filterTimes()
        else:
            self.d = self.rawD.copy()


        return

    def _filterTimes(self):
        """
        This function will filter by times.
        """
        print(self.tRange)
        idT = np.where((self.tRange[0] > np.array(self.rawD['Epoch'][:])) & 
                    (self.tRange[1] < np.array(self.rawD['Epoch'][:])))[0]
        #print(self.rawD['Epoch'][:100])
        print(idT)
        # Filter data
        for key in filter(lambda x: ('Epoch' in x or 
                        ('Counts' in x and x[-1] == 's')), self.rawD.keys()):
            self.d[key] = self.rawD[key].copy()[idT]
        return

if __name__ == '__main__':
    tRange = [datetime(2017, 3, 31, 11), datetime(2017, 3, 31, 11, 18)]
    obj = Proc_L0_RBSPICE(('/home/mike/research/rbsp/data/rbspice/rbspa/'
                        'rbsp-a-rbspice_lev-1_EBR_20170331_v1.1.2-01'),
                        tRange=None)

    #print(obj.d.keys())

    for SSD in obj.d['EBR'].T:
        plt.plot(obj.d['Epoch'], SSD)
    plt.show()