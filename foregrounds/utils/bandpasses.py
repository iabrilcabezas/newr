'''
bandpasses
'''

import numpy as np
import utils.sed as used
from utils.params import PATH_DICT


band_names = used.get_band_names()



def write_bandpass_experiment(freq):

    '''
    writes delta bandpass in the absence of experiment bandpass
    '''

    freq_array = np.array([freq - 0.5, freq, freq + 0.5])
    bps_array  = np.array([0.,1.,0.])

    full_array = np.vstack((freq_array, bps_array)).T
    np.savetxt(PATH_DICT['bbpower_path'] +\
                             f'examples/data/bandpasses/{freq}.txt', full_array)

def import_bandpasses():

    '''
    Returns dictionary of bandpass class for each channel in EXPERIMENT
    '''

   # if EXPERIMENT == 'S4_SAT':
   #     freqs = used.S4_frequencies()
   #     bpss = {band_names[i]: used.Bpass_delta(n) for i, n in enumerate(freqs)}
    
    #if EXPERIMENT == 'so':
    bpss = {n: used.Bpass_band(n,PATH_DICT['bbpower_path'] +\
                             f'examples/data/bandpasses/{n}.txt') for n in band_names}


    return bpss
