'''
params.py

reads-in parameters from config file
'''

import yaml

class GlobalConfig:

    '''
    2nd level class of config file
    contains info on global parameters
    '''

    def __init__(self, param):
        self.r     = param['r']
        self.nside = param['nside']
        self.Lmax  = param['Lmax']
        self.coord_out = param['coord_out']
        self.out_path  = param['output_path']

class DelensConfig:

    '''
    2nd level class of config file
    contains info on delensing parameters
    '''

    def __init__(self, param):
        self.alpha_phi = param['alpha_phi']
        self.alpha_psi = param['alpha_psi']

class EllConfig:

    '''
    2nd level class of config file
    contains info on bandpower range and binning
    '''

    def __init__(self, param):
        self.lmax_out = param['lmax_out']
        self.lmin_out = param['lmin_out']
        self.dell_out = param['dell_out']
        self.lmax_phi = param['lmax_phi']
        self.lmin_phi = param['lmin_phi']
        self.lmax_psi = param['lmax_psi']
        self.lmin_psi = param['lmin_psi']

class ConstantsConfig:

    '''
    2nd level class of config file
    contains value of physical constants
    '''

    def __init__(self, param):
        self.Tcmb = param['Tcmb']

class ForegroundConfig:

    '''
    2nd level class of config file
    contains foreground especifications
    '''

    def __init__(self, param):
        self.mask = param['mask']
        self.dust = param['dust']



class Config:

    '''
    1st level class for config file
    '''

    def __init__(self, param):
        self.global_param = GlobalConfig(param['global'])
        self.delens_param = DelensConfig(param['delens'])
        self.ell_param    = EllConfig(param['ell'])
        self.constants    = ConstantsConfig(param['constants'])
        self.foreground   = ForegroundConfig(param['foreground'])

with open('./QEfgs/config.yml', 'r', encoding = 'utf-8') as config_file:
    config = Config(yaml.load(config_file, yaml.FullLoader))

global_params    = config.global_param
delens_params    = config.delens_param
ell_params       = config.ell_param
constants        = config.constants
foreground       = config.foreground

R         = global_params.r
NSIDE     = global_params.nside
LMAX      = global_params.Lmax
COORD_OUT = global_params.coord_out
OUTPUT_PATH = global_params.out_path

ALPHA_PHI = delens_params.alpha_phi
ALPHA_PSI = delens_params.alpha_psi

LMAX_OUT = ell_params.lmax_out
LMIN_OUT = ell_params.lmin_out
DELL_OUT = ell_params.dell_out

LMAX_PHI = ell_params.lmax_phi
LMIN_PHI = ell_params.lmin_phi

LMAX_PSI = ell_params.lmax_psi
LMIN_PSI = ell_params.lmin_psi

TCMB = float(constants.Tcmb)

MASK = foreground.mask
DUST = foreground.dust

FOOTPRINT      = MASK['footprint']
FOOTPRINT_PATH = MASK['footprint_path']
APO_DEG        = MASK['apo_deg']
COORD_MASK     = MASK['coord']

DUST_TYPE      = DUST['type']
DUST_FREQ      = DUST['freq']
DUST_PATH_BASE = DUST['path']
COORD_DUST     = DUST['coord']

if DUST_TYPE.split('_')[0] == 'pysm':
    DUST_TYPE, DUST_SUBTYPE =DUST_TYPE.split('_')
else:
    DUST_SUBTYPE = '0'

FOOTPRINT_TYPES = ['so']
DUST_TYPES      = ['DF', 'van', 'pysm']

SO_FOOTPRINT_PATH = '/global/cfs/cdirs/act/data/iabril/SO/mask_apodized_david_nside512.fits'

DF_NAMEFILE = f'DustFilaments_TQU_NS2048_Nfil180p5M_BKpatchNormalization_{int(DUST_FREQ)}p0GHz.fits'
VAN_NAMEFILE = f'vans_d1_SOS4_{int(DUST_FREQ):03}_tophat_map_2048.fits'

## assertions:
assert 3*NSIDE + 1 > LMAX, 'not enough nside resolution to reach LMAX'
assert LMAX > LMAX_OUT, 'global Lmax is less than Lmax_out'
assert LMAX > LMAX_PHI, 'global Lmax is less than Lmax_phi'
assert LMAX > LMAX_PSI, 'global Lmax is less than Lmax_psi'

assert FOOTPRINT in FOOTPRINT_TYPES, 'undefined footprint'
assert DUST_TYPE in DUST_TYPES, 'undefined dust foreground'

if FOOTPRINT == 'so':
    assert FOOTPRINT_PATH == SO_FOOTPRINT_PATH
    assert COORD_MASK == 'C'

if DUST_TYPE == 'DF':
    assert COORD_DUST == 'G'
elif DUST_TYPE == 'van':
    assert COORD_DUST == 'G'
elif DUST_TYPE == 'pysm':
    assert COORD_DUST == 'G'

if DUST_TYPE == 'DF':
    DUST_PATH = DUST_PATH_BASE + 'DF/' + DF_NAMEFILE
elif DUST_TYPE == 'van':
    DUST_PATH = DUST_PATH_BASE + 'Vansyngel/' + VAN_NAMEFILE
elif DUST_TYPE == 'pysm':
    DUST_PATH = 'NA'

NAME_RUN = f'{LMAX}_{FOOTPRINT}_{DUST_TYPE}_{DUST_SUBTYPE}_{int(DUST_FREQ)}'
