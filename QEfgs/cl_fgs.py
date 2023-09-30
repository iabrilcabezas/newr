'''
plot and compare power spectra of foreground maps
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
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


# # change config.yml to measure spectra from different dust types
# from QEfgs.fgs.fgsreadin import cellfromalm_mask
# cellfromalm_mask()



# # plot cell results:
# freq_dict = {'DF_0': 95, 'pysm_d1': 90, 'van_0':90}
# colour = ['cornflowerblue', 'gold','forestgreen']
# labels = [r'$\rm{DUSTFILAMENTS}$', r'$\rm{PySM - d1}$', r'$\rm{Vansyngel}$']

# L_array = get_ell_arrays(LMAX)[0]

# for pol in ['EE','EB','BB']:

#     fig, ax = plt.subplots()

#     for d, DUST in enumerate(['DF_0', 'pysm_d1','van_0']):
#         cell = np.loadtxt(OUTPUT_PATH + f'{LMAX}_{FOOTPRINT}_{DUST}_{freq_dict[DUST]}_cl_{pol}.txt')

#         if pol == 'EB':
#             negative = cell < 0
#             ax.semilogy(L_array[negative], np.abs(cell[negative]), linestyle = 'dashed', color = colour[d])
#             ax.semilogy(L_array[~negative], cell[~negative], label = labels[d], color = colour[d])

#         else:
#             ax.semilogy(L_array, cell, label = labels[d], color = colour[d])
#         ax.legend(prop={'size':15})
#         ax.set_xlabel(r'$\ell$')
#         ax.set_ylabel(r'$C_\ell$' + f'[{pol}]')
#     fig.savefig(f'./QEfgs/Figures/cl_{pol}.png', bbox_inches = 'tight')
#     plt.show()

# measure maps:
import healpy as hp
from QEfgs.fgs.fgsreadin import read_dustmap, read_mask

labels = ['T', 'Q', 'U']

def plot_maps():
    '''plot maps'''

    dustT, dustQ, dustU = read_dustmap()

    mask = read_mask()[0]

    dustT_mask = dustT * mask
    dustQ_mask = dustQ * mask
    dustU_mask = dustU * mask

    for i, dustcomp in enumerate([dustT_mask, dustQ_mask, dustU_mask]):
        fig = plt.figure(1, dpi = 200, facecolor='w')
        hp.mollview(dustcomp, coord = [COORD_OUT, 'C'],bgcolor='w', title = labels[i] + ' map') 
        plt.savefig(f'./QEfgs/Figures/dust_map_{labels[i]}.png', bbox_inches = 'tight')
        plt.show()

plot_maps()
