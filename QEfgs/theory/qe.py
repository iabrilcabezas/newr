'''
qe.py
    calling quadratic estimators
'''

import numpy as np
import curvedsky as cs
from QEfgs.theory.theory_cell import get_observed_spectrum, theory_camb, get_EBlm_plm_theo_dellrange
from QEfgs.fgs.fgsreadin import almEB_maskdust_raw_ellrange
# from QEfgs.utils.write import read_complex
from QEfgs.utils.params import LMAX_OUT, LMIN_PHI, LMAX_PHI, LMIN_PSI, LMAX_PSI
from QEfgs.utils.params import NSIDE, RECONS_PSI, RECONS_PHI, BASE_NAME_THEO, OUTPUT_PATH

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
        Fl_psi[:, l,:l+1] = 1./ocl_psi[:3, l, None] # TT, EE, BB

    return Fl_psi

def get_upsi(recons_data):

    '''
    reconstruct psi (unnormalized)
    rec_lens.qbb(BB spectrum, inverse-variance filtered B-mode alm)
    recons_data alms
    '''

    # filtering
    Fl_psi = get_Fl_psi()

    if RECONS_PSI == 'BB':

        Blm = recons_data['psi']['B'] # (psi, BB) component
        fBlm   = Blm * Fl_psi[2, :, :] ## BB component

        rcl = theory_camb()[-1] # tensor output
        # grad component
        psi_bb = cs.rec_lens.qbb(LMAX_OUT, LMIN_PSI, LMAX_PSI,rcl[2,:LMAX_PSI+1],fBlm,fBlm)[0]
        return psi_bb

    return None

def get_uphi(recons_data):

    '''
    reconstruct phi (unnormalized)
    rec_lens.qeb(EE spectrum, inverse-variance filtered E-mode alm, inverse B)

    recons_data alms
    '''

    Fl_phi = get_Fl_phi()

    if RECONS_PHI == 'EB':

        Elm = recons_data['phi']['E']
        Blm = recons_data['phi']['B']

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

def get_rphi(recons_data):

    '''
    reconstruct phi = norm * QE
    '''
    A = norm_phi()
    phi = get_uphi(recons_data)

    return A[:LMAX_OUT+1, None] * phi

def get_rpsi(recons_data):

    '''
    reconstruct psi = norm * QE
    '''

    A = norm_psi()
    psi = get_upsi(recons_data)

    return A[:LMAX_OUT+1, None] * psi

def reconstruct_on(recons_type):

    '''
    define which alms are being used for reconstruction
    recons_type ['dust', 'cmb']
    '''

    if recons_type == 'dust':
        reconstruct_alms = almEB_maskdust_raw_ellrange()
    elif recons_type == 'cmb':
        reconstruct_alms = get_EBlm_plm_theo_dellrange()
    else:
        return None
    return reconstruct_alms

def rec_psipsi(recons_data):

    '''
    reconstruct psi and obtain psi x psi cell
    '''

    psilm = get_rpsi(recons_data)
    gg = cs.utils.alm2cl(LMAX_OUT, psilm)

    return psilm, gg

def rec_phiphi(recons_data):
    '''
    reconstruct phi and obtain phi x phi cell
    '''
    philm = get_rphi(recons_data)
    pp_recon = cs.utils.alm2cl(LMAX_OUT, philm)

    return philm, pp_recon

def recpsi_theory(recons_data):

    '''
    reconstruct psi and obtain psi x phi theory cell
    '''
    psilm = rec_psipsi(recons_data)[0]
    # palm = read_complex(f'{BASE_NAME_THEO}_palm')
    palm = np.load(f'{OUTPUT_PATH}{BASE_NAME_THEO}_palm')
    gp_theo = cs.utils.alm2cl(LMAX_OUT,psilm,palm[:LMAX_OUT+1,:LMAX_OUT+1])

    return gp_theo

def recphi_theory(recons_data):

    '''
    reconstruct phi and obtain psi x phi theory cell
    '''

    philm = rec_phiphi(recons_data)[0]
    # palm = read_complex(f'{BASE_NAME_THEO}_palm')
    palm = np.load(f'{OUTPUT_PATH}{BASE_NAME_THEO}_palm')
    pg_ptheo = cs.utils.alm2cl(LMAX_OUT, philm, palm[:LMAX_OUT+1,:LMAX_OUT+1])

    return pg_ptheo
