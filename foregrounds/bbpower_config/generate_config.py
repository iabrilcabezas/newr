'''
generate_config.py

generates config file for BBCompSep
'''
import yaml

from utils.params import amp_d_bb, alpha_d_bb, l0_d_bb, beta_d, temp_d, nu0_d
from utils.params import r_tensor, A_lens
from utils.params import amp_s_bb, alpha_s_bb, l0_s_bb, beta_s, nu0_s
from utils.params import epsilon_ds
from utils.params import decorr_amp_d, decorr_amp_s, decorr_nu01_d, decorr_nu01_s, decorr_nu02_d, decorr_nu02_s

from utils.params import PATH_DICT
from utils.params import LMIN, LMAX
from utils.params import NSIDE

def name_config(cmb, dust, sync, decorr, cross):

    name = ''
    if dust:
        name += 'd'
    else:
        name += '0'
    if cmb:
        name += 'c'
    else:
        name += '0'
    if sync:
        name += 's'
    else:
        name += '0'
    if decorr:
        name += '_D'
    else:
        name += '_0'
    if cross:
        name += '_C'
    else:
        name += '_0'

    name += f'_{NSIDE}_{LMIN}_{LMAX}'

    return name

def generate_cell_model_dict(cmb, dust, sync, decorr, cross):

    if not cmb:

        rtensor = 0.
        Alens   = 0.

    else:
        rtensor = r_tensor
        Alens = A_lens

    dict_cmbmodel = {'cmb_templates':   [PATH_DICT['bbpower_path'] + 'examples/data/camb_lens_nobb.dat',
                                         PATH_DICT['bbpower_path'] + 'examples/data/camb_lens_r1.dat'],
                    'params': { 'r_tensor': ['r_tensor','tophat',  [-1, float(rtensor), 1]],
                                'A_lens':   ['A_lens',  'tophat',  [0.0, float(Alens), 2.0]]}}

    dict_comp1 = {'name': 'Dust',
                  'sed': 'Dust',
                  'cl': {'BB': 'ClPowerLaw'},
                  'sed_parameters': {'beta_d': ['beta_d', 'Gaussian', [ float(beta_d), 0.5]],
                                     'temp_d': ['temp', 'fixed', [ float(temp_d)]],
                                     'nu0_d':  ['nu0',  'fixed', [ float(nu0_d)]]},
                   'cl_parameters': {'BB': {'amp_d_bb': ['amp', 'tophat', [0.0, float(amp_d_bb), 'inf']],
                                            'alpha_d_bb': ['alpha', 'tophat', [-1.0, float(alpha_d_bb), 0.0]],
                                            'l0_d_bb': ['ell0', 'fixed', [float(l0_d_bb)]]}}}

    dict_comp2 = {  'name': 'Synchrotron',
                    'sed': 'Synchrotron',
                    'cl': { 'BB': 'ClPowerLaw'},
                    'sed_parameters': { 'beta_s': ['beta_pl', 'Gaussian', [ float(beta_s), 0.6]],
                                        'nu0_s' : ['nu0', 'fixed', [float(nu0_s)]]},
                    'cl_parameters': {'BB': {   'amp_s_bb': ['amp', 'tophat', [0., float(amp_s_bb), "inf"]],
                                                'alpha_s_bb': ['alpha', 'tophat', [-1., float(alpha_s_bb), 0.]],
                                                'l0_s_bb': ['ell0', 'fixed', [float(l0_s_bb)]] }}
                }
    
    
    if decorr:
        dict_comp1['decorr'] = {'decorr_amp_d': ['decorr_amp', 'fixed', [decorr_amp_d]],
                               'decorr_nu01_d': ['decorr_nu01', 'fixed', [decorr_nu01_d]],
                               'decorr_nu02_d': ['decorr_nu02', 'fixed', [decorr_nu02_d]]}
        dict_comp2['decorr'] = {'decorr_amp_s': ['decorr_amp', 'fixed', [decorr_amp_s]],
                               'decorr_nu01_s': ['decorr_nu01', 'fixed', [decorr_nu01_s]],
                               'decorr_nu02_s': ['decorr_nu02', 'fixed', [decorr_nu02_s]]}
    if cross:
        assert (dust is True) & (sync is True), 'correlation with only one component present'
        dict_comp1['cross'] = { 'epsilon_ds': ['component_2', 'tophat', [-1., float(epsilon_ds), 1.]]}

    if dust:
        dict_fgmodel = {'component_1' : dict_comp1}
    else:
        dict_fgmodel = {'component_1': dict_comp2}
    if (dust & sync):
        dict_fgmodel['component_2'] = dict_comp2

    return dict_cmbmodel, dict_fgmodel


def generate_bbcompsep_dict(cmb, dust, sync, decorr, cross):

    dict_cmbmodel, dict_fgmodel = generate_cell_model_dict(cmb, dust, sync, decorr, cross)

    path_config  = PATH_DICT['output_path'] + name_config(cmb, dust, sync, decorr, cross) + '.yml'

    dict_bbcompsep = {  'global':  {'nside': int(NSIDE), 'compute_dell': False},
                        'modules': 'bbpower',
                        'launcher': 'local',
                        'stages': [{'name': 'BBCompSep', 'nprocess': 1}],
                        'config': path_config,
                        }
    
    dict_bbcompsep['BBCompSep'] = {'sampler': 'predicted_spectra',
                                    'predict_at_minimum': False,
                                    'predict_to_sacc': True,
                                    'pol_channels': ['B'],
                                    'l_min': int(LMIN),
                                    'l_max': int(LMAX),
                                    'bands': 'all',
                                    'cmb_model': dict_cmbmodel,
                                    'fg_model': dict_fgmodel,
                                    }

    with open(path_config, 'w', encoding = 'utf-8') as file:
        yaml.dump(dict_bbcompsep, file, default_flow_style=True)

    return dict_bbcompsep
