
import numpy as np
from utils.params import POLARIZATION
import utils.sed as used
from utils.bandpowers import dell2cell_lmax
from sacc_object.generate_sacc import import_bandpasses, get_windows

band_names = used.get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)
ncombs = len(indices_tr[0])

ctype_dict_ncomp = {'d00': 1, 'dc0': 2, 'dcs': 3,'00s':1}

def get_bpw_freq_sig(ctype, lmax, do_bin, weight = 'Cl'):

    '''
    Computes SED of all components
    Convolves SED with instrument (bands, CMB units)
    Computes total signal in each bandpower

    ** Parameters **
    ctype: 'd00', 'dc0', 'dcs'
        type of model (dust only, dust + CMB, dust + CMB + sync)
    lmax: int
        max ell to compute SED
    do_bin: bool
        bin signal with window?
    weight: 'Cl' or 'Dl'
        binning of window
    '''

    assert POLARIZATION == 'B', 'reading B components but you have specified otherwise'

    bpss = import_bandpasses()

    dl2cl = dell2cell_lmax(lmax)

    ncomp = ctype_dict_ncomp[ctype]

    dls_comp = np.zeros([ncomp,nmodes,ncomp,nmodes,lmax+1]) #[ncomp,np,ncomp,np,nl]
    dls_sed  = used.get_component_spectra(lmax)

    if ctype == 'd00':
        dls_comp[0,0,0,0,:] = dls_sed[0]
    if ctype == '00s':
        dls_comp[0,0,0,0,:] = dls_sed[2]

    if ctype == 'dc0':
        (dls_comp[1,0,1,0,:],
        dls_comp[0,0,0,0,:]) = dls_sed[0], dls_sed[1]

    if ctype == 'dcs':
        (dls_comp[1,0,1,0,:],
        dls_comp[0,0,0,0,:],
        dls_comp[2,0,2,0,:]) = dls_sed[0], dls_sed[1], dls_sed[2]

    dls_comp *= dl2cl[None, None, None, None, :]

    if do_bin:
        windows = get_windows(weight)
        bpw_comp=np.sum(dls_comp[:,:,:,:,None,:]*windows[None,None,None, None, :,:],axis=5)
        to_bpw_freq = bpw_comp
    else:
        to_bpw_freq = dls_comp

    # Convolve with bandpasses
    seds = used.get_convolved_seds(band_names, bpss)

    if ctype == 'd00':
        seds = np.array([seds[1,:]])
    elif ctype == '00s':
        seds = np.array([seds[2,:]])
    else:
        print('no sed for ctype defined')
        return None

    bpw_freq_sig = np.einsum('ik,jm,iljno', seds, seds, to_bpw_freq)

    return bpw_freq_sig

