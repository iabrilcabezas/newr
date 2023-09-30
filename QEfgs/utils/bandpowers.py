''''
bandpowers utility
'''

import numpy as np
import binning as bins
from QEfgs.utils.params import LMAX_OUT, LMIN_OUT, DELL_OUT

def get_ell_arrays(lmax):

    '''
    defines ell array and ell factor for given lmax
    '''

    L = np.linspace(0, lmax, lmax+1)
    Lfac = (L * (L+1))**2 / (2*np.pi)

    return L, Lfac

def get_ell_bins():

    '''
    defines binning for bandpowers
    '''

    bin_scheme = bins.multipole_binning(DELL_OUT, lmin = LMIN_OUT, lmax = LMAX_OUT)

    return bin_scheme