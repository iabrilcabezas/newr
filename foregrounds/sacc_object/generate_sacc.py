'''
generate_sacc.py
'''
import numpy as np
import sacc
import utils.noise_calc as nc
from utils.params import EXPERIMENT, PATH_DICT, POLARIZATION
from utils.params import LMAX, NBANDS, LBANDS, LARR_ALL, DELL, LEFF, LMIN
import utils.sed as used

band_names = used.get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)
ncombs = len(indices_tr[0])

def get_windows(weight):

    '''
    Returns window with binning according to weights

    ** Parameters **
    weight: 'Dl' or 'Cl'
        weight as Dell [l * (l + 1)/ 2pi] or 'Cell' [equal weights]

    ** Returns **
    windows: np.array([NBANDS, LMAX + 1])
    '''

    weight_types = ['Dl', 'Cl']
    assert weight in weight_types, 'not a type of weight!'

    windows = np.zeros([NBANDS,LMAX+1])

    cl_weights = np.ones_like(LARR_ALL)

    for b_i,(b_l0,b_lf) in enumerate(zip(LBANDS[:-1],LBANDS[1:])):

        if weight == 'Dl':
            windows[b_i,b_l0:b_lf] = (LARR_ALL * (LARR_ALL + 1)/(2*np.pi))[b_l0:b_lf]
        if weight == 'Cl':
            windows[b_i, b_l0:b_lf] = cl_weights[b_l0:b_lf]

        windows[b_i,:] /= DELL

    return windows

def import_beams(ell_array):

    '''
    Returns dictionary of beams for each channel in EXPERIMENT, evaluated at each ell in ell_array

    ** Parameters **
    ell_array: np.array()
        array where to evaluate the beams
    '''

    if EXPERIMENT == 'S4_SAT': 
        beams ={band_names[i]: b for i, b in \
                    enumerate(nc.S4_SAT_beams(ell_array))}
    if EXPERIMENT  == 'so':
        beams ={band_names[i]: b for i, b in \
                    enumerate(nc.Simons_Observatory_V3_SA_beams(ell_array))}


    return beams

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

def add_tracers(ell_array):

    '''
    Creates sacc object and add tracers according to EXPERIMENT
    Beam of tracer is evaluated at ell = ell_array

    ** Parameters **
    ell_array: np.array()
        ell array where beams are calculated

    ** Returns **
    s_d: sacc.Sacc()
        sacc object with added tracers
    '''

    s_d = sacc.Sacc()

     # Bandpasses:
    bpss = import_bandpasses()
    # Beams
    beams = import_beams(ell_array)

    for i_band, name_b in enumerate(band_names):
        bandpass = bpss[name_b]
        beam = beams[name_b]
        s_d.add_tracer('NuMap', f'band{i_band+1}',
                        quantity='cmb_polarization',
                        spin=2,
                        nu=bandpass.nu,
                        bandpass=bandpass.bnu,
                        ell=ell_array,
                        beam=beam,
                        nu_unit='GHz',
                        map_unit='uK_CMB')

    return s_d


def add_powerspectra(s_d):

    '''
    Adds power spectra to Sacc object

    ** Parameters **
    s_d: sacc object
        object to add P(k)
    bpw_freq_sig: np.array
        power spectra that will be added
    leff: np.array
        ell each power spectra band corresponds to
    do_bin: bool
        if the power spectra is binned, provide window too
    weight: 'Cl' or 'Dl'
        binnin of window if do_bin
    '''

    windows = get_windows('Cl')
    s_wins = sacc.BandpowerWindow(LARR_ALL, windows.T)

    map_names=[]
    for i_freq in range(nfreqs):
        map_names.append(f'band{i_freq+1}_B')

    for (i_tr1, i_tr2) in zip(indices_tr[0], indices_tr[1]):
        band12 = np.array([map_names[i_tr1][:-2], map_names[i_tr2][:-2]])
        pol12  = np.array([map_names[i_tr1][-1].lower(), map_names[i_tr2][-1].lower()])
        cl_type = f'cl_{pol12[0]}{pol12[1]}'

        s_d.add_ell_cl(cl_type, band12[0], band12[1], LEFF, np.ones_like(LEFF),
                        window = s_wins)

    return s_d

def write_tracer():

    '''
    write fake sacc object with all the info we need to generate good object
    '''

    s_d = add_tracers(LARR_ALL)
    # fake data
    s_d = add_powerspectra(s_d)

    fake_cov = np.diag(np.ones(ncombs * NBANDS))
    s_d.add_covariance(fake_cov)
    #print(s_d.get_ell_cl('cl_bb', 'band1', 'band1'))
    #print(s_d.tracers.items())

    s_d.save_fits(PATH_DICT['output_path'] + f'sacc_{EXPERIMENT}_{LMIN}_{LMAX}_{DELL}.fits', overwrite = True)

