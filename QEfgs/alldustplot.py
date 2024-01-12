'''
plot results
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.ticker as mtick
from matplotlib import rcParams
import binning as bins
from QEfgs.fgs.fgsreadin import read_mask
from QEfgs.utils.bandpowers import get_ell_arrays, get_ell_bins
from QEfgs.utils.params import LMAX_OUT, BASE_NAME_THEO, R, NAME_THEO, NAME_RUN, OUTPUT_PATH, NAME_LRANGES

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
# rcParams['font.family'] = 'Helvetica'

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

cmb_name = f'recons/{NAME_THEO}_{NAME_LRANGES}'
dust_name = f'recons/dust_{NAME_RUN}_{NAME_LRANGES}'

n = 3
colours = pl.cm.viridis(np.linspace(0,1,n))
#colours = ['forestgreen','crimson','cornflowerblue']

def get_plot_data():

    '''
    plot reconstruction
    '''

    # palm      =  read_complex(f'{BASE_NAME_THEO}_palm')
    # pp_theory = cs.utils.alm2cl(LMAX_OUT,palm)
    pp_theory   = np.load(f'{OUTPUT_PATH}{BASE_NAME_THEO}_pp_theory.npy')[:LMAX_OUT+1]
    rpp_theory = R * pp_theory

    w2_factor, w4_factor = read_mask()[1:]

    cmbname = f'{OUTPUT_PATH}recons/{NAME_THEO}_{NAME_LRANGES}'
    dustname1 = f'{OUTPUT_PATH}recons/dust_3050_so_pysm_d1_90_{NAME_LRANGES}'
    dustname2 = f'{OUTPUT_PATH}recons/dust_3050_so_DF_0_95_{NAME_LRANGES}'
    dustname3 = f'{OUTPUT_PATH}recons/dust_3050_so_van_0_90_{NAME_LRANGES}'

    plot_data = {}
    
    plot_type = 'pp'
    data_plot_cmb  = np.load(f'{cmbname}_pp_recon.npy')
    data_plot_dust1 = np.load(f'{dustname1}_pp_recon.npy')
    data_plot_dust2 = np.load(f'{dustname2}_pp_recon.npy')
    data_plot_dust3 = np.load(f'{dustname3}_pp_recon.npy')
    # apply correction factor:
    data_plot_dust1 /= w4_factor
    data_plot_dust2 /= w4_factor
    data_plot_dust3 /= w4_factor

    theory = pp_theory

    plot_data[plot_type] = {'cmb': data_plot_cmb, 'theory': theory, 'dust1': data_plot_dust1, 'dust2': data_plot_dust2,'dust3':data_plot_dust3}


    plot_type = 'gg'

    data_plot_cmb  = np.load(f'{cmbname}_gg.npy')
    data_plot_dust1 = np.load(f'{dustname1}_gg.npy')
    data_plot_dust2 = np.load(f'{dustname2}_gg.npy')
    data_plot_dust3 = np.load(f'{dustname3}_gg.npy')
    # apply correction factor:
    data_plot_dust1 /= w4_factor
    data_plot_dust2 /= w4_factor
    data_plot_dust3 /= w4_factor

    theory = rpp_theory

    plot_data[plot_type] = {'cmb': data_plot_cmb, 'theory': theory, 'dust1': data_plot_dust1, 'dust2': data_plot_dust2,'dust3':data_plot_dust3}

    plot_type = 'gp'

    data_plot_cmb  = np.load(f'{cmbname}_gp_theo.npy')
    data_plot_dust1 = np.load(f'{dustname1}_gp_theo.npy')
    data_plot_dust2 = np.load(f'{dustname2}_gp_theo.npy')
    data_plot_dust3 = np.load(f'{dustname3}_gp_theo.npy')
    # apply correction factor:
    data_plot_dust1 /= w2_factor
    data_plot_dust2 /= w2_factor
    data_plot_dust3 /= w2_factor

    theory = rpp_theory

    plot_data[plot_type] = {'cmb': data_plot_cmb, 'theory': theory, 'dust1': data_plot_dust1, 'dust2': data_plot_dust2,'dust3':data_plot_dust3}

    plot_type = 'pg'
    data_plot_cmb  = np.load(f'{cmbname}_pg_theo.npy')
    data_plot_dust1 = np.load(f'{dustname1}_pg_theo.npy')
    data_plot_dust2 = np.load(f'{dustname2}_pg_theo.npy')
    data_plot_dust3 = np.load(f'{dustname3}_pg_theo.npy')
    # apply correction factor:
    data_plot_dust1 /= w2_factor
    data_plot_dust2 /= w2_factor
    data_plot_dust3 /= w2_factor

    theory = pp_theory

    plot_data[plot_type] = {'cmb': data_plot_cmb, 'theory': theory, 'dust1': data_plot_dust1, 'dust2': data_plot_dust2,'dust3':data_plot_dust3}

    plot_type = 'pg_recons'
    data_plot_cmb  = np.load(f'{cmbname}_gp_recon.npy')
    data_plot_dust1 = np.load(f'{dustname1}_gp_recon.npy')
    data_plot_dust2 = np.load(f'{dustname2}_gp_recon.npy')
    data_plot_dust3 = np.load(f'{dustname3}_gp_recon.npy')
    # apply correction factor:
    data_plot_dust1 /= w4_factor
    data_plot_dust2 /= w4_factor
    data_plot_dust3 /= w4_factor

    theory = rpp_theory

    plot_data[plot_type] = {'cmb': data_plot_cmb, 'theory': theory, 'dust1': data_plot_dust1, 'dust2': data_plot_dust2,'dust3':data_plot_dust3}

    return plot_data

def plotting_routine(plot_type, plot_data, save = True):

    '''
    plotting routine
    '''

    _, Lfac = get_ell_arrays(LMAX_OUT)
    mb = get_ell_bins()

    # binned:
    fig, ax = plt.subplots(dpi = 144)



    ydata_plot_dust1 = bins.binning(Lfac * plot_data[plot_type]['dust1'], mb)
    ydata_plot_dust2 = bins.binning(Lfac * plot_data[plot_type]['dust2'], mb)
    ydata_plot_dust3 = bins.binning(Lfac * plot_data[plot_type]['dust3'], mb)
        
    ax.plot(mb.bc, ydata_plot_dust1, label = label_dust[label_dict[plot_type]] + r'$\textrm{PySM}$', color = colours[0], linestyle = 'solid')
    ax.plot(mb.bc, ydata_plot_dust2, label = label_dust[label_dict[plot_type]]+ r'$\textrm{DF}$', color = colours[0], linestyle = 'dotted')
    ax.plot(mb.bc, ydata_plot_dust3, label = label_dust[label_dict[plot_type]]+ r'$\textrm{van}$', color = colours[0], linestyle = 'dashed')

    ydata_plot_cmb = bins.binning(Lfac * plot_data[plot_type]['cmb'], mb)
 
    ax.plot(mb.bc, ydata_plot_cmb, label = label_cmb[label_dict[plot_type]], color = colours[1])
    ax.plot(mb.bc, bins.binning(Lfac * plot_data[plot_type]['theory'], mb), label = label_theory[label_dict[plot_type]], color = colours[2])

    if (np.any(ydata_plot_dust3<0) or np.any(ydata_plot_cmb<0)):
        ax.axhline(0, color = 'gray', linestyle = 'dashed')
        ax.set_yscale('symlog')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter(r'${%.1e}$'))
    else:
        ax.set_yscale('log')
    #ax.locator_params(axis = 'y', nbins = 3)
  

    ax.legend(prop = {'size':15}, ncol= 2)
    ax.set_xlabel(r'$L$')
    ax.set_ylabel(r'$L^2 \left(L+1\right)^2 C_\ell / 2\pi$')
    name_fig = f'{plot_type}_{NAME_THEO}'
    name_fig += f'_{NAME_LRANGES}'
    if save:
        fig.savefig(f'./QEfgs/Figures/all_dust_{name_fig}.pdf', bbox_inches = 'tight')
    plt.show()

def make_plots():
    data_plot = get_plot_data()
    for a in ['pg_recons', 'pp']: #['pp','gg','pg','gp','pg_recons']
        plotting_routine(a, data_plot)

make_plots()


