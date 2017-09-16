# These functions read in the burst data from EMFSIS
import sys, os
import numpy as np
import dateutil.parser
from datetime import datetime

try:
    sys.path.insert(0, '/home/mike/Dropbox/0_grad_work')
    import directories
except SystemError:
    print('Could not import the directories.py file, '
        'please supply data directories manualy!')

def loadFile(f):
    """
    This generator function will yeild lines from the file one by one.
    It will parse time with the parseTime() function and for each line output,
    attach a data type string indicator of 'time', 'frequencies', or 'spectra'.
    """
    dataType = 'frequencies'

    with open(f) as handle:
        for i, line in enumerate(handle):
            splitLine = line.split()

            # If the line is a time, return it.
            if '2017' in splitLine:
                yield 'time', parseTime(splitLine)
                continue
            if '2000' in splitLine: # From this point, the data are all Nan's
                raise ValueError('Valid data ends here')
                
            
            # Try to parse the data into floats and yeild it along with the
            # data type that starts with 
            try:
                yield dataType, list(map(float, splitLine))
            except ValueError as err:
                # After this line in the data, spectra starts.
                if 'Time' in str(err): 
                    dataType = 'spectra'
                continue                

def parseTime(timeStr):
    """
    This is the helper function for loadFile() that parses Dave's time
    stamps and yeilds them as datetime objects. 
    """
    timeStr[:-1] = list(map(int, timeStr[:-1]))
                
    # Convert decimal seconds seconds and milliseconds
    timeStr[-1] = float(timeStr[-1])
    timeStr.append(int(1E6*(timeStr[-1] - int(timeStr[-1])))) 
    # Seconds
    timeStr[-2] = int(timeStr[-2])
    return datetime(*timeStr[1:])

if __name__ == '__main__':

    sc_id = 'A'
    times = np.nan*np.ones(10200, dtype = object)
    spectra = np.nan*np.ones((10200, 512), dtype = float) 
    frequencies = np.nan*np.ones(512, dtype = float)
    idf = 0; idt = 0; # Keep track of which row we are saving the data to.

    fname = 'RBSP-A_20170331T11FFT1024overlap000_calibratedwaveformburst_v1.0.0.txt'

    fileGenerator = loadFile(os.path.join(directories.emfisis_dir(sc_id), fname))

    # The try statement will attemp to loop over the file, and stop when it 
    # raises a ValueError indicating the following data is nans
    try: 
        for i, line in enumerate(fileGenerator):

            # The if statements catch the generator output which yeilds the 
            # type of data as the first tag and redirects the output to the 
            # appropriate array.
            if line[0] == 'time': # Save the times
                times[idt] = line[1]
                idt += 1
                ids = 0 # Reset (or start) the spectra column index.

            elif line[0] == 'frequencies': # Save the frequencies
                frequencies[idf : idf+len(line[1])] = line[1]
                idf += len(line[1]) # Increment the frequency index.

            # Save the spectra. This is tricky since its a 2d array.
            elif line[0] == 'spectra': 
                #print(line[1])
                spectra[idt, ids : ids+len(line[1])] = line[1]
                ids += len(line[1])
    except ValueError:
        # Filter the data by nans
        validt = np.where(~np.isnan(spectra[:, 0]))[0]
        times = times[validt]
        spectra = spectra[validt, :]

        # Save the arrays to binary numoy format
        np.save('../data/emfisis_burst_times', times)
        np.save('../data/emfisis_burst_spectra', spectra)
        np.save('../data/emfisis_burst_frequencies', frequencies)

