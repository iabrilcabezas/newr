'''
plot and compare power spectra of foreground maps
'''
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import rcParams
from QEfgs.fgs.fgsreadin import read_dustmap, read_mask
from QEfgs.fgs.fgsreadin import cellfromalm_mask, get_EBlm_write
from QEfgs.utils.params import OUTPUT_PATH, LMAX, FOOTPRINT, COORD_OUT
from QEfgs.utils.bandpowers import get_ell_arrays
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
rcParams['font.family'] = 'DejaVu Sans'

TASK = 'compute_alm'

TASK_TYPES = ['compute_alm', 'measure_cell', 'plot_cell', 'plot_map']
assert TASK in TASK_TYPES, 'task undefined'

# plot cell results:
freq_dict = {'DF_0': 95, 'pysm_d1': 90, 'van_0':90}
colour = ['cornflowerblue', 'gold','forestgreen']
labels_cell = [r'$\rm{DUSTFILAMENTS}$', r'$\rm{PySM - d1}$', r'$\rm{Vansyngel}$']
labels_map = ['T', 'Q', 'U']
L_array = get_ell_arrays(LMAX)[0]

# measure maps:
def plot_cell():
    '''plots cell of fgs'''
    for pol in ['EE','EB','BB']:

        fig, ax = plt.subplots()

        for d, DUST in enumerate(['DF_0', 'pysm_d1','van_0']):
            cell = np.loadtxt(OUTPUT_PATH + f'fgs/{LMAX}_{FOOTPRINT}_{DUST}_{freq_dict[DUST]}_cl_{pol}.txt')

            if pol == 'EB':
                negative = cell < 0
                ax.semilogy(L_array[negative], np.abs(cell[negative]), linestyle = 'dashed', color = colour[d])
                ax.semilogy(L_array[~negative], cell[~negative], label = labels_cell[d], color = colour[d])

            else:
                ax.semilogy(L_array, cell, label = labels_cell[d], color = colour[d])
            ax.legend(prop={'size':15})
            ax.set_xlabel(r'$\ell$')
            ax.set_ylabel(r'$C_\ell$' + f'[{pol}]')
        fig.savefig(f'./QEfgs/Figures/cl_{pol}.png', bbox_inches = 'tight')
        plt.show()

def plot_maps():
    '''plot maps'''

    dust_t, dust_q, dust_u = read_dustmap()

    mask = read_mask()[0]

    dust_t_mask = dust_t * mask
    dust_q_mask = dust_q* mask
    dust_u_mask = dust_u * mask

    for i, dustcomp in enumerate([dust_t_mask, dust_q_mask, dust_u_mask]):
        plt.figure(1, dpi = 200, facecolor='w')
        hp.mollview(dustcomp, coord = [COORD_OUT, 'C'],bgcolor='w', title = labels_map[i] + ' map')
        plt.savefig(f'./QEfgs/Figures/dust_map_{labels_map[i]}.png', bbox_inches = 'tight')
        plt.show()

if TASK == 'measure_cell':
    # change config.yml to measure spectra from different dust types
    cellfromalm_mask()

if TASK == 'plot_cell':
    plot_cell()

if TASK == 'plot_map':
    plot_maps()

if TASK == 'compute_alm':
    get_EBlm_write()