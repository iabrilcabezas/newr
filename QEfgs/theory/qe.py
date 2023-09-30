'''
qe.py
    calling quadratic estimators
'''

import numpy as np
import curvedsky as cs
from QEfgs.theory.theory_cell import get_observed_spectrum, theory_camb
from QEfgs.fgs.fgsreadin import almEB_maskdust_raw_ellrange
from QEfgs.utils.params import LMAX_OUT, LMIN_PHI, LMAX_PHI, LMIN_PSI, LMAX_PSI
from QEfgs.utils.params import NSIDE, RECONS_PSI, RECONS_PHI

# filtering

def get_Fl_phi():

    '''
    filtering for phi = 1/ observed phi
    '''

    ocl_phi = get_observed_spectrum()['phi']

    # do 1/Cell filtering:
    Fl_phi = np.zeros((3, LMAX_PHI+1, LMAX_PHI + 1)) # (T, E, B)
    for l in range(LMIN_PHI, LMAX_PHI+1):
        Fl_phi[:,l,:l+1] = 1./ocl_phi[:3, l, None] # TT, EE, BB

    return Fl_phi

def get_Fl_psi():

    '''
    filtering for psi = 1/ observed psi
    '''

    ocl_psi = get_observed_spectrum()['psi']

    Fl_psi = np.zeros((3, LMAX_PSI+1, LMAX_PSI + 1)) # (T, E, B)
    for l in range(LMIN_PSI, LMAX_PSI+1):
        Fl_psi[l,:l+1] = 1./ocl_psi[:3, l, None] # TT, EE, BB

    return Fl_psi

def get_upsi():

    '''
    reconstruct psi (unnormalized)
    rec_lens.qbb(BB spectrum, inverse-variance filtered B-mode alm)
    '''

    # filtering
    Fl_psi = get_Fl_psi()

    if RECONS_PSI == 'BB':

        Blm = almEB_maskdust_raw_ellrange()['psi']['BB'] # (psi, BB) component
        fBlm   = Blm * Fl_psi[2, :, :] ## BB component

        rcl = theory_camb()[-1] # tensor output
        # grad component
        psi_bb = cs.rec_lens.qbb(LMAX_OUT, LMIN_PSI, LMAX_PSI,rcl[2,:LMAX_PSI+1],fBlm,fBlm)[0]
        return psi_bb

    return None

def get_uphi():

    '''
    reconstruct phi (unnormalized)
    rec_lens.qeb(EE spectrum, inverse-variance filtered E-mode alm, inverse B)
    '''

    Fl_phi = get_Fl_phi()

    if RECONS_PHI == 'EB':

        Elm = almEB_maskdust_raw_ellrange()['phi']['EE']
        Blm = almEB_maskdust_raw_ellrange()['phi']['BB']

        fElm = Elm * Fl_phi[1,:,:] # EE
        fBlm = Blm * Fl_phi[2,:,:] # BB

        lcl = theory_camb()[1] # lensed cl output
        # grad component
        phi_qeb = cs.rec_lens.qeb(LMAX_OUT, LMIN_PHI, LMAX_PHI, lcl[1,:LMAX_PHI+1], \
                                  fElm, fBlm, nside_t = NSIDE)[0]

        return phi_qeb

    return None

def norm_psi():

    '''
    normalization QE for psi
    norm_quad.qbb(theory BB spectrum, observed bb spectrum)
    '''

    ocl_psi = get_observed_spectrum()['psi']

    if RECONS_PSI == 'BB':
        rcl = theory_camb()[-1] # tensor
        Aglm_bb = cs.norm_quad.qbb('lens',LMAX_OUT, LMIN_PSI, LMAX_PSI,\
                                   rcl[2,:LMAX_PSI+1],ocl_psi[2,:LMAX_PSI+1])[0]

        return Aglm_bb

    return None

def norm_phi():

    '''
    normalization QE for psi
    norm_quad.qeb(theory EE spectrum, observed EE spectrum, observed BB spectrum)
        ignore theory BB spectrum as it is much smaller than EE
    '''

    ocl_phi = get_observed_spectrum()['phi']

    if RECONS_PHI == 'EB':
        lcl, _ = theory_camb()[1:]

        Aglm_eb = cs.norm_quad.qeb('lens', LMAX_OUT, LMIN_PHI, LMAX_PHI, lcl[1, :LMAX_PHI+1],\
                                    ocl_phi[1,:LMAX_PHI + 1], ocl_phi[2,:LMAX_PHI+1])[0]

        return Aglm_eb

    return None

def get_rphi():

    '''
    reconstruct phi = norm * QE
    '''
    A = norm_phi()
    phi = get_uphi()

    return A[:LMAX_OUT+1, None] * phi

def get_rpsi():

    '''
    reconstruct psi = norm * QE
    '''

    A = norm_psi()
    psi = get_upsi()

    return A[:LMAX_OUT+1, None] * psi
