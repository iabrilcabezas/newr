'''
write.py
    write/read to folders
'''

import numpy as np
from QEfgs.utils.params import OUTPUT_PATH

def write_complex(name, array):
    '''
    write complex array to file 
    '''

    # assert np.all(np.iscomplex(array)), 'this funtion is to save complex arrays only'
    np.savetxt(OUTPUT_PATH + f'{name}_real.txt', array.real)
    np.savetxt(OUTPUT_PATH + f'{name}_imag.txt', array.imag)


def read_complex(name):
    '''
    read array that is saved as complex
    '''

    real_array = np.loadtxt(OUTPUT_PATH + f'{name}_real.txt')
    imag_array = np.loadtxt(OUTPUT_PATH + f'{name}_imag.txt')

    array = np.array(real_array + 1j * imag_array, dtype=complex)

    return array

def write_cell(name, array):

    '''
    simple write 1D array to .txt file
    '''

    np.savetxt(OUTPUT_PATH + f'{name}_cell.txt', array)

def read_cell(name):
    '''
    read cell array
    '''
    array = np.loadtxt(name)
    return array