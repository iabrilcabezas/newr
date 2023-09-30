'''
fgsreadin.py
    reads-in foregrounds
'''
import numpy as np
import healpy as hp
import pymaster as nmt
import pysm3
import pysm3.units as u
import curvedsky as cs
from QEfgs.utils.params import DUST_PATH, COORD_DUST, DUST_TYPE, DUST_SUBTYPE
from QEfgs.utils.params import NSIDE, TCMB, COORD_OUT, LMAX, DUST_FREQ
from QEfgs.utils.params import FOOTPRINT_PATH, COORD_MASK, FOOTPRINT, APO_DEG
from QEfgs.utils.params import OUTPUT_PATH, NAME_RUN


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
    reads-in mask and also returns w2 factor
    '''
    # raw map
    w_mask = hp.read_map(FOOTPRINT_PATH)
    # rotate
    w_mask = hp_rotate(w_mask, coord = [COORD_MASK, COORD_OUT])
    w_mask = hp.ud_grade(w_mask, NSIDE)

    if FOOTPRINT == 'so':
        print('fsky = 0.1')
        so_patch = w_mask > 0.35
        mask_so_now = np.zeros_like(w_mask)
        mask_so_now[so_patch] = 1
        mask_so_apo = nmt.mask_apodization(mask_so_now, APO_DEG , apotype="C1") # apo_deg

        return mask_so_apo, np.mean(mask_so_apo**2)

    return None

def almEB_maskdust_raw():

    '''
    measures a_lm polarizaiton coefficients of masked dust map
    '''

    dustQ, dustU = read_dustmap()[1:]

    mask = read_mask()[0]

    dustQ_mask = dustQ * mask
    dustU_mask = dustU * mask

    almE, almB = cs.utils.hp_map2alm_spin( NSIDE, LMAX, mmax = LMAX, spin=2, \
                                          map0 = dustQ_mask, map1 = dustU_mask)

    return almE, almB

def cellfromalm_mask():
    '''
    measures power spectra from alms. corrects for mask with w2 factor
    '''
    almE, almB = almEB_maskdust_raw()
    w2factor = read_mask()[1]

    cl_EE = cs.utils.alm2cl(LMAX, almE, almE)
    cl_EB = cs.utils.alm2cl(LMAX, almE, almB)
    cl_BB = cs.utils.alm2cl(LMAX, almB, almB)
    # w2 correction
    cl_EE /= w2factor
    cl_EB /= w2factor
    cl_BB /= w2factor

    np.savetxt(OUTPUT_PATH + NAME_RUN + '_cl_EE.txt', cl_EE)
    np.savetxt(OUTPUT_PATH + NAME_RUN + '_cl_EB.txt', cl_EB)
    np.savetxt(OUTPUT_PATH + NAME_RUN + '_cl_BB.txt', cl_BB)
