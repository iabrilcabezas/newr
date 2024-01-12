'''
fgsreadin.py
    reads-in foregrounds
'''
from pathlib import Path
import numpy as np
import healpy as hp
import pymaster as nmt
import pysm3
import pysm3.units as u
import curvedsky as cs
import cmb
#from QEfgs.utils.write import read_complex, write_complex
from QEfgs.utils.sed import fdust
from QEfgs.utils.params import DUST_PATH, COORD_DUST, DUST_TYPE, DUST_SUBTYPE
from QEfgs.utils.params import NSIDE, COORD_OUT, LMAX, DUST_FREQ
from QEfgs.utils.params import FOOTPRINT_PATH, COORD_MASK, FOOTPRINT
from QEfgs.utils.params import OUTPUT_PATH, NAME_RUN
from QEfgs.utils.params import LMAX_PHI, LMAX_PSI

# constants
TCMB = cmb.Tcmb

base_name = 'fgs/' + NAME_RUN
cell_base_name = 'fgs_cell/' + NAME_RUN

def hp_rotate(map_hp, coord):
    """Rotate healpix map between coordinate systems

    :param map_hp: A healpix map in RING ordering
    :param coord: A len(2) list of either 'G', 'C', 'E'
    Galactic, equatorial, ecliptic, eg ['G', 'C'] converts
    galactic to equatorial coordinates
    :returns: A rotated healpix map
    """
    if map_hp is None:
        return None
    if coord[0] == coord[1]:
        return map_hp
    rotator_func = hp.rotator.Rotator(coord=coord)
    new_map = rotator_func.rotate_map_pixel(map_hp)
    return new_map

def read_dustmap():

    '''
    reads-in dust map
    '''

    if DUST_TYPE == 'pysm':

        sky = pysm3.Sky(nside = NSIDE, preset_strings = [DUST_SUBTYPE])
        dustmap = sky.get_emission(DUST_FREQ * u.GHz)

    else:
        dustmap = hp.read_map(DUST_PATH, field = None )
        if DUST_TYPE == 'planck':
            # SCALE FREQ DOWN:
            dustmap *= fdust(DUST_FREQ) * 1e6 # from 353GHz (where dust = 1) to any other frequency, and K to muK

    dustmap_T, dustmap_Q, dustmap_U = dustmap[0], dustmap[1], dustmap[2]

    # downgrade and normalize:
    dustmap_T = hp.ud_grade(dustmap_T, NSIDE) / TCMB
    dustmap_Q = hp.ud_grade(dustmap_Q, NSIDE) / TCMB
    dustmap_U = hp.ud_grade(dustmap_U, NSIDE) / TCMB

    # rotate:
    dustmap_T = hp_rotate(dustmap_T, coord = [COORD_DUST, COORD_OUT])
    dustmap_Q = hp_rotate(dustmap_Q, coord = [COORD_DUST, COORD_OUT])
    dustmap_U = hp_rotate(dustmap_U, coord = [COORD_DUST, COORD_OUT])

    return dustmap_T, dustmap_Q, dustmap_U

def read_mask():

    '''
    reads-in mask and also returns w2, w4 factors

    Standard procedure is:
    1. Read mask
    2. Rotate mask
    3. Downgrade/upgrade mask
    4. Transform to binary mask (in case it goes with e.g. number hit count)
    '''

    if FOOTPRINT.split('_')[0] == 'planck':

        fsky_dict = {'20': 0, '40': 1, '60': 2,'70':3, '80':4 ,'90':5, '97':6}

        w_mask = hp.read_map(FOOTPRINT_PATH, field = fsky_dict[FOOTPRINT.split('_')[1]])
        
        w_mask = hp_rotate(w_mask, coord = [COORD_MASK, COORD_OUT])
        mask = hp.ud_grade(w_mask, NSIDE)

    elif FOOTPRINT == 'so':

        w_mask = hp.read_map(FOOTPRINT_PATH)
        # rotate
        w_mask = hp_rotate(w_mask, coord = [COORD_MASK, COORD_OUT])
        w_mask = hp.ud_grade(w_mask, NSIDE)
        
        so_patch = w_mask > 0.35
        mask_so_now = np.zeros_like(w_mask)
        mask_so_now[so_patch] = 1
        mask = nmt.mask_apodization(mask_so_now, 5.0 , apotype="C1") # apo_deg # apodized


    elif FOOTPRINT.split('_')[0] == 'ACT':

        fsky = FOOTPRINT_PATH.split("GAL")[1][:3]
        assert FOOTPRINT.split('_')[1] == fsky, 'fsky does not match PATH to ACT mask'

        w_mask = hp.read_map(FOOTPRINT_PATH)
        w_mask = hp_rotate(w_mask, coord =[COORD_MASK, COORD_OUT])
        mask = hp.ud_grade(w_mask, NSIDE)

    
    elif FOOTPRINT == 'BICEP':
        bicep3_mask = hp.read_map(FOOTPRINT_PATH)
        bicep3_mask = np.nan_to_num(bicep3_mask)
        bicep3_mask = hp_rotate(bicep3_mask, coord = [COORD_MASK, COORD_OUT])

        # make binary mask
        binary = bicep3_mask > 0.001
        bicep3_mask[~binary] = 0.
        bicep3_mask[binary] = 1.
        bicep_apo = nmt.mask_apodization(bicep3_mask, 5., apotype="C1")

        mask = hp.ud_grade(bicep_apo, NSIDE)

    elif FOOTPRINT == 'none':
        mask = np.ones(hp.nside2npix(NSIDE))

    else:
        print('unrecognized footprint')
        return None
    print(f'{FOOTPRINT} fsky = {np.mean(mask):.2f}')

    return mask, np.mean(mask**2), np.mean(mask**4)

def almEB_maskdust_raw():

    '''
    measures a_lm polarizaiton coefficients of masked dust map
    '''

    dustT, dustQ, dustU = read_dustmap()

    mask = read_mask()[0]

    dustT_mask = dustT * mask
    dustQ_mask = dustQ * mask
    dustU_mask = dustU * mask

    almT = cs.utils.hp_map2alm(NSIDE, LMAX, mmax = LMAX, map = dustT_mask)
    almE, almB = cs.utils.hp_map2alm_spin( NSIDE, LMAX, mmax = LMAX, spin=2, \
                                          map0 = dustQ_mask, map1 = dustU_mask)

    return almT, almE, almB

def cellfromalm_mask():
    '''
    measures power spectra from alms. corrects for mask with w2 factor
    '''

    if Path(f'{OUTPUT_PATH}{base_name}_almT.npy').exists():

        almT = np.load(f'{OUTPUT_PATH}{base_name}_almT.npy')
        almE = np.load(f'{OUTPUT_PATH}{base_name}_almE.npy')
        almB = np.load(f'{OUTPUT_PATH}{base_name}_almB.npy')
    
    else:
        almT, almE, almB = almEB_maskdust_raw()

    w2factor = read_mask()[1]

    cl_TT = cs.utils.alm2cl(LMAX, almT, almT)
    cl_EE = cs.utils.alm2cl(LMAX, almE, almE)
    cl_EB = cs.utils.alm2cl(LMAX, almE, almB)
    cl_BB = cs.utils.alm2cl(LMAX, almB, almB)
    # w2 correction
    cl_TT /= w2factor
    cl_EE /= w2factor
    cl_EB /= w2factor
    cl_BB /= w2factor

    np.savetxt(OUTPUT_PATH + cell_base_name + '_cl_TT.txt', cl_TT)
    np.savetxt(OUTPUT_PATH + cell_base_name + '_cl_EE.txt', cl_EE)
    np.savetxt(OUTPUT_PATH + cell_base_name + '_cl_EB.txt', cl_EB)
    np.savetxt(OUTPUT_PATH + cell_base_name + '_cl_BB.txt', cl_BB)

def get_EBlm_write():
    '''
    computes and writes foreground alms
    '''

    almT, almE, almB = almEB_maskdust_raw()
    # save it:
    # write_complex(f'{base_name}_almT', almT)
    # write_complex(f'{base_name}_almE', almE)
    # write_complex(f'{base_name}_almB', almB)
    np.save(f'{OUTPUT_PATH}{base_name}_almT', almT)
    np.save(f'{OUTPUT_PATH}{base_name}_almE', almE)
    np.save(f'{OUTPUT_PATH}{base_name}_almB', almB)


def almEB_maskdust_raw_ellrange():

    '''
    returns needed alms for reconstruction in the required ell range
    '''

    # almT = read_complex(f'{base_name}_almT')
    # almE = read_complex(f'{base_name}_almE')
    # almB = read_complex(f'{base_name}_almB')
    almT = np.load(f'{OUTPUT_PATH}{base_name}_almT.npy')
    almE = np.load(f'{OUTPUT_PATH}{base_name}_almE.npy')
    almB = np.load(f'{OUTPUT_PATH}{base_name}_almB.npy')

    # phi:
    Tlm_phi = almT[:LMAX_PHI+1, :LMAX_PHI+1]
    Elm_phi = almE[:LMAX_PHI+1, :LMAX_PHI+1]
    Blm_phi = almB[:LMAX_PHI+1, :LMAX_PHI+1]
    # psi
    Tlm_psi = almT[:LMAX_PSI+1, :LMAX_PSI+1]
    Elm_psi = almE[:LMAX_PSI+1, :LMAX_PSI+1]
    Blm_psi = almB[:LMAX_PSI+1, :LMAX_PSI+1]

    return {'phi': {'T' :Tlm_phi, 'E': Elm_phi, 'B': Blm_phi},\
             'psi': {'T': Tlm_psi, 'E': Elm_psi, 'B': Blm_psi}}
