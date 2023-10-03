'''
once palm, theory is measured, and alm maps too: reconstruct
'''
import numpy as np
from QEfgs.theory.qe import reconstruct_on
from QEfgs.theory.qe import rec_psipsi, rec_phiphi, recphi_theory, recpsi_theory
# from QEfgs.utils.write import write_cell, write_complex
from QEfgs.utils.params import NAME_RUN, CMB_SEED, OUTPUT_PATH

def run_reconstruction(recons_data, name):

    '''
    runs reconstruction, saves to file
    '''

    philm, pp_recon = rec_phiphi(recons_data)
    psilm, gg       = rec_psipsi(recons_data)
    gp_theo         = recpsi_theory(recons_data)
    pg_theo         = recphi_theory(recons_data)

    # write_complex(f'{name}_philm', philm)
    # write_complex(f'{name}_psilm', psilm)
    # write_cell(f'{name}_gg', gg)
    # write_cell(f'{name}_pp_recon', pp_recon)
    # write_cell(f'{name}_gp_theo', gp_theo)
    # write_cell(f'{name}_pg_theo', pg_theo)
    np.save(f'{OUTPUT_PATH}{name}_philm', philm)
    np.save(f'{OUTPUT_PATH}{name}_psilm', psilm)
    np.save(f'{OUTPUT_PATH}{name}_gg', gg)
    np.save(f'{OUTPUT_PATH}{name}_pp_recon', pp_recon)
    np.save(f'{OUTPUT_PATH}{name}_gp_theo', gp_theo)
    np.save(f'{OUTPUT_PATH}{name}_pg_theo', pg_theo)

dust_data = reconstruct_on('dust')
cmb_data  = reconstruct_on('cmb')

cmb_name = f'recons/cmb_{int(CMB_SEED):03}'
dust_name = f'recons/dust_{NAME_RUN}'

run_reconstruction(dust_data, dust_name)
run_reconstruction(cmb_data, cmb_name)
