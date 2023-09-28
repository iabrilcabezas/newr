
import numpy as np
from astropy.io import fits
from sed_byhand.compute_freq_sig import get_bpw_freq_sig
import utils.sed as used
from utils.params import POLARIZATION, LMAX, PATH_DICT
from utils.params import EXPERIMENT, LMIN, NSIDE
from utils.params import amp_s_bb, alpha_s_bb, amp_d_bb, alpha_d_bb

band_names = used.get_band_names()
nfreqs = len(band_names)
nmodes = len(POLARIZATION)
nmaps=nmodes*nfreqs
indices_tr=np.triu_indices(nmaps)
ncombs = len(indices_tr[0])

print(float(abs(alpha_d_bb)))
for ctype in ['d00', '00s']:

    bpw_freq_sig = get_bpw_freq_sig(ctype = ctype, lmax = LMAX, do_bin = True, weight = 'Cl')
    nells = bpw_freq_sig.shape[-1]
    bpw_freq_sig = bpw_freq_sig.reshape([nfreqs*nmodes,nfreqs*nmodes, nells])

    hdu_w = fits.PrimaryHDU(bpw_freq_sig)
    hdu_w.writeto(PATH_DICT['output_path'] + f'{EXPERIMENT}_{NSIDE}_{LMIN}_{LMAX}_{nells}_{ctype}_' + \
                   f'{float(abs(amp_s_bb)):.2f}_{float(abs(alpha_s_bb)):.2f}_{float(abs(amp_d_bb)):.2f}_{float(abs(alpha_d_bb)):.2f}' + \
                    '_cell.fits', overwrite = True)
