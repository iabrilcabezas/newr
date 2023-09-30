'''
theory_cell
    functions from theory
'''
import numpy as np
import camb
import cmb
from QEfgs.utils.params import LMAX
from QEfgs.utils.params import ALPHA_PHI, ALPHA_PSI, R

def prepare_cls(lmax=4000,Alens=1.):

    '''
    Borrowed from T. Namikawa

    Computes unlensed, lensed and tensor cells in the desired lmax range from CAMB

    NB. CAMB Cl arrays are always in (TT, EE, BB, TE)
    '''

    #Set up a new set of parameters for CAMB
    pars = camb.CAMBparams()
    # This function sets up CosmoMC-like settings
    # One massive neutrino and helium set using BBN consistency
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06, Alens=Alens)
    pars.InitPower.set_params(As=2e-9, ns=0.965, r=1)
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    pars.WantTensors = True
    pars.NonLinear = camb.model.NonLinear_both
    results = camb.get_transfer_functions(pars)
    powers =results.get_cmb_power_spectra(pars, CMB_unit='muK')

    # unlensed Cls
    ls = np.arange(powers['unlensed_scalar'].shape[0])
    fac = ls*(ls+1)/(2*np.pi)*cmb.Tcmb**2
    ucl = np.zeros((4,lmax+1))
    ucl[0,2:] = powers['unlensed_scalar'][2:lmax+1,0]/fac[2:lmax+1] # TT
    ucl[1,2:] = powers['unlensed_scalar'][2:lmax+1,1]/fac[2:lmax+1] # EE
    ucl[2,2:] = powers['unlensed_scalar'][2:lmax+1,3]/fac[2:lmax+1] # TE
    ucl[3,2:] = powers['lens_potential'][2:lmax+1,0]/(ls[2:lmax+1]**2*(ls[2:lmax+1]+1)**2)*(2*np.pi)

    # lensed Cls
    ls = np.arange(powers['lensed_scalar'].shape[0])
    fac = ls*(ls+1)/(2*np.pi)*cmb.Tcmb**2
    lcl = np.zeros((4,lmax+1))
    lcl[0,2:] = powers['lensed_scalar'][2:lmax+1,0]/fac[2:lmax+1] # TT
    lcl[1,2:] = powers['lensed_scalar'][2:lmax+1,1]/fac[2:lmax+1] # EE
    lcl[2,2:] = powers['lensed_scalar'][2:lmax+1,2]/fac[2:lmax+1] # BB
    lcl[3,2:] = powers['lensed_scalar'][2:lmax+1,3]/fac[2:lmax+1] # EB

    # tensor
    tCL = powers['tensor']
    ls  = np.arange(tCL.shape[0])
    fac = ls*(ls+1)/(2*np.pi)*cmb.Tcmb**2
    rcl = np.zeros((4,lmax+1))
    rcl[1,2:] = tCL[2:lmax+1,1]/fac[2:lmax+1] # EE
    rcl[2,2:] = tCL[2:lmax+1,2]/fac[2:lmax+1] # BB
    # ucl[TT, EE, TE], lcl[TT, EE, BB, EB], rcl[0, EE, BB, 0]
    return ucl, lcl, rcl

def theory_camb():
    ''' returns theory cells '''
    ucl, lcl, rcl = prepare_cls(lmax=LMAX)
    return ucl, lcl, rcl


def alens_alpha():

    ''' computes Alens from alpha'''

    alens_psi = (1 - ALPHA_PSI)**2
    alens_phi = (1 - ALPHA_PHI)**2
    return alens_psi, alens_phi

def get_observed_spectrum():

    ''''
    Observed spectrum:
    Cell = r * tensor modes + alens * scalar modes
    '''

    alens_psi, alens_phi = alens_alpha()
    lcl, rcl = theory_camb()[1:]

    ocl_phi = R*rcl + lcl*alens_phi
    ocl_psi = R*rcl + lcl*alens_psi

    return {'phi': ocl_phi, 'psi': ocl_psi}
