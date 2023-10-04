'''
plot results
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import binning as bins
import curvedsky as cs
from QEfgs.utils.bandpowers import get_ell_arrays, get_ell_bins
from QEfgs.utils.write import read_complex, read_cell
from QEfgs.utils.params import LMAX_OUT, BASE_NAME_THEO, R, CMB_SEED, NAME_RUN, OUTPUT_PATH

rcParams['font.size']=24
rcParams['axes.linewidth']=1.5
rcParams['xtick.major.width']=1.5
rcParams['xtick.minor.width']=1
rcParams['ytick.major.width']=1.5
rcParams['ytick.minor.width']=1
rcParams['xtick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['ytick.direction'] = 'in'
rcParams['xtick.minor.visible'] = True
rcParams['ytick.minor.visible'] = True
rcParams['ytick.right'] = True
rcParams['text.usetex'] = True
rcParams['font.family'] = 'Helvetica'

label_rpp = r'$rC_\ell^{\phi\phi}$'
label_pp  = r'$C_\ell^{\phi\phi}$'

label_dict = {'gg':0, 'pp':1, 'gp':2, 'pg':3, 'pg_recons':4}

label_theory = [label_rpp, label_pp, label_rpp, label_pp, label_rpp]
label_cmb = [r'$C_\ell^{\hat{\psi}\hat{\psi}}$', r'$C_\ell^{\hat{\phi}\hat{\phi}}$', \
             r'$C_\ell^{\hat{\psi}{\phi}}$', r'$C_\ell^{\hat{\phi}{\phi}}$',\
             r'$C_\ell^{\hat{\psi}\hat{\phi}}$']
label_dust = [r'$C_\ell^{\hat{\psi}_d\hat{\psi}_d}$', r'$C_\ell^{\hat{\phi}_d\hat{\phi}_d}$', \
             r'$C_\ell^{\hat{\psi}_d{\phi}}$', r'$C_\ell^{\hat{\phi}_d{\phi}}$', \
             r'$C_\ell^{\hat{\psi}_d\hat{\phi}_d}$']

# {'gg': [r'$C_\ell^{{\hat{\psi}}_d{\hat{\psi}}_d}$', label_rpp] ,
#                 'pp': [r'$C_\ell^{{\hat{\phi}}_d{\hat{\phi}}_d}$', label_pp],
#                 'gp': [r'$C_\ell^{{\hat{\psi}}_d{\phi}}$', label_rpp],
#                 'pg': [r'$C_\ell^{{\hat{\phi}}_d{\phi}}$', label_pp],
#                 'pg_recons': [r'$C_\ell^{{\hat{\psi}}_d{\hat{\phi}}_d}$', label_rpp]}

cmb_name = f'recons/cmb_{int(CMB_SEED):03}'
dust_name = f'recons/dust_{NAME_RUN}'

colours = ['forestgreen','crimson','cornflowerblue']

def plotting_routine(plot_type, include_dust, data_plot_dust, data_plot_cmb, theory):

    '''
    plotting routine
    '''

    _, Lfac = get_ell_arrays(LMAX_OUT)
    mb = get_ell_bins()

    # binned:
    fig, ax = plt.subplots(dpi = 144)
    if include_dust:
        ydata_plot_dust = bins.binning(Lfac * data_plot_dust, mb)
        negative = ydata_plot_dust < 0
        ax.semilogy(mb.bc[negative], np.abs(ydata_plot_dust[negative]), linestyle = 'dashed', color = colours[0])
        ax.semilogy(mb.bc[~negative], ydata_plot_dust[~negative], label = label_dust[label_dict[plot_type]], color = colours[0])

    ydata_plot_cmb = bins.binning(Lfac * data_plot_cmb, mb)
    negative = ydata_plot_cmb < 0
    ax.semilogy(mb.bc[negative], np.abs(ydata_plot_cmb[negative]), linestyle = 'dashed', color = colours[1])
    ax.semilogy(mb.bc[~negative], ydata_plot_cmb[~negative], label = label_cmb[label_dict[plot_type]], color = colours[1])

    ax.semilogy(mb.bc, bins.binning(Lfac * theory, mb), label = label_theory[label_dict[plot_type]], color = colours[2])
    ax.legend(prop = {'size':15})
    ax.set_xlabel(r'$L$')
    ax.set_ylabel(r'$L^2 \left(L+1\right)^2 / 2\pi \cdot$' + r'$C_\ell$')
    name_fig = f'{plot_type}_{CMB_SEED}'
    if include_dust:
        name_fig += NAME_RUN
    fig.savefig(f'./QEfgs/Figures/{name_fig}.png', bbox_inches = 'tight')
    plt.show()




def get_plot_data(plot_type):

    '''
    plot reconstruction
    '''

    # palm      =  read_complex(f'{BASE_NAME_THEO}_palm')
    # pp_theory = cs.utils.alm2cl(LMAX_OUT,palm)
    pp_theory   = np.load(f'{OUTPUT_PATH}{BASE_NAME_THEO}_pp_theory.npy')[:LMAX_OUT+1]
    rpp_theory = R * pp_theory

    if plot_type == 'pp':

        # data_plot_cmb = read_cell(f'{cmb_name}_pp_recon')
        # data_plot_dust = read_cell(f'{dust_name}_pp_recon')
        data_plot_cmb  = np.load(f'{OUTPUT_PATH}{cmb_name}_pp_recon.npy')
        data_plot_dust = np.load(f'{OUTPUT_PATH}{dust_name}_pp_recon.npy')
        theory = pp_theory

    if plot_type == 'gg':

        # data_plot_cmb = read_cell(f'{cmb_name}_gg')
        # data_plot_dust = read_cell(f'{dust_name}_gg')
        data_plot_cmb  = np.load(f'{OUTPUT_PATH}{cmb_name}_gg.npy')
        data_plot_dust = np.load(f'{OUTPUT_PATH}{dust_name}_gg.npy')
        theory = rpp_theory

    if plot_type == 'gp':

        # data_plot_cmb = read_cell(f'{cmb_name}_gp_theo')
        # data_plot_dust = read_cell(f'{dust_name}_gp_theo')
        data_plot_cmb  = np.load(f'{OUTPUT_PATH}{cmb_name}_gp_theo.npy')
        data_plot_dust = np.load(f'{OUTPUT_PATH}{dust_name}_gp_theo.npy')
        theory = rpp_theory

    if plot_type == 'pg':

        # data_plot_cmb = read_cell(f'{cmb_name}_pg_theo')
        # data_plot_dust = read_cell(f'{dust_name}_pg_theo')
        data_plot_cmb  = np.load(f'{OUTPUT_PATH}{cmb_name}_pg_theo.npy')
        data_plot_dust = np.load(f'{OUTPUT_PATH}{dust_name}_pg_theo.npy')
        theory = pp_theory

    if plot_type == 'pg_recons':

        # philm_cmb = read_complex(f'{cmb_name}_philm')
        # psilm_cmb = read_complex(f'{cmb_name}_psilm')
        # data_plot_cmb = cs.utils.alm2cl(LMAX_OUT,psilm_cmb,philm_cmb)
        # philm_dust = read_complex(f'{dust_name}_philm')
        # psilm_dust = read_complex(f'{dust_name}_psilm')
        # data_plot_dust = cs.utils.alm2cl(LMAX_OUT,psilm_dust,philm_dust)
        data_plot_cmb  = np.load(f'{OUTPUT_PATH}{cmb_name}_gp_recon.npy')
        data_plot_dust = np.load(f'{OUTPUT_PATH}{dust_name}_gp_recon.npy')
        theory = rpp_theory

    return data_plot_cmb, data_plot_dust, theory


for a in ['pp','gg','pg','gp','pg_recons']:
    for dustbool in [True, False]:
        data_plot_cmb, data_plot_dust, theory = get_plot_data(a)
        plotting_routine(a, dustbool, data_plot_dust, data_plot_cmb, theory)
