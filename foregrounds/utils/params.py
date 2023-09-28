'''
parameter file read-in
'''
import yaml
from utils.bandpowers import get_ell_arrays

def pathnames():

    '''
    Returns path to data depending on machine

    ** Parameters **
    machine: str
            machine where code runs ('cori' or 'perl')

    ** Returns **
    dict: stores path to data
        'planck_path': Planck data
        'camb_cmb_lens_nobb': CMB power spectrum, including lensing, without bb
        'output_path': location of file output storage
        'input_path': location of file input storage
        'BK15_data': BK15 data (bandpasses, noise)
        'bbpipe_path': path to bbpower code. contains SO bandpasses and noise in examples folder
        'so_path': path where SO mask is stored
    '''

    dict_path = {}

    
    bbpower_path = config.path_param.bbpower
    output_path = config.path_param.output
     
    dict_path['bbpower_path'] = bbpower_path
    dict_path['output_path'] = output_path
    

    return dict_path

class PathConfig:

    '''
    2nd level class of config file
    contains path info

    bbpower: str
        path to bbpower code
    '''

    def __init__(self, param):

        self.bbpower  = param['bbpower']
        self.output   = param['output']


class CosmoConfig:

    '''
    2nd level class of config file
    contains info on cosmo_params
        CMB: CMB parametrization
        dust: dust parametrization
        model: model parameters
    '''

    def __init__(self, param):
        self.CMB = param['CMB']
        self.dust = param['dust']
        self.sync = param['sync']
        self.model = param['model']

class BpwConfig:

    '''
    2nd level class info on config file
    contains bandpower parameters

    lmin: int
        min ell
    nbands: int
        number of ell bands
    dell: int
        width of ell bands

    '''
    def __init__(self, param):
        self.lmin = param['lmin']
        self.nbands = param['nbands']
        self.dell = param['dell']

class GlobalConfig:

    '''
    2nd level class info on config file
    contains global config parameters

    nside: 2**N
        resolution of maps

    '''

    def __init__(self, param):
        self.nside      = param['nside']
        self.experiment = param['experiment']
        self.polarization = param['polarization']

class BandConfig:

    '''
    2nd level class of config file
    contains info on band_names
    s4sat: list
        CMB S4 SAT band names
    '''

    def __init__(self, param):
        self.s4sat = param['S4_SAT']
        self.so    = param['SO']
        self.lbrd  = param['lbrd']

class Config:

    '''
    1st level class for config file
    cosmo_param: parameters on CMB and dust P(k)

    '''

    def __init__(self, param):
        self.cosmo_param  = CosmoConfig(param['cosmology'])
        self.path_param   = PathConfig(param['paths'])
        self.bpw_param    = BpwConfig(param['bandpowers'])
        self.global_param = GlobalConfig(param['global'])
        self.band_names   = BandConfig(param['band_names'])

with open('config.yml', 'r', encoding = 'utf-8') as config_file:
    config = Config(yaml.load(config_file, yaml.FullLoader))

cosmo_params    = config.cosmo_param
cmb_params      = cosmo_params.CMB
dust_params     = cosmo_params.dust
sync_params     = cosmo_params.sync
model_params    = cosmo_params.model

amp_d_bb        = dust_params['amp_d_bb']
alpha_d_bb      = dust_params['alpha_d_bb']
l0_d_bb         = dust_params['l0_d_bb']
beta_d          = dust_params['beta_d']
temp_d          = dust_params['temp_d']
nu0_d           = dust_params['nu0_d']

A_lens          = cmb_params['A_lens']
r_tensor        = cmb_params['r_tensor']

epsilon_ds      = model_params['epsilon_ds']

amp_s_bb        = sync_params['amp_s_bb']
alpha_s_bb      = sync_params['alpha_s_bb']
l0_s_bb         = sync_params['l0_s_bb']
beta_s          = sync_params['beta_s']
nu0_s           = sync_params['nu0_s']

decorr_amp_d    = dust_params['decorr_amp_d']
decorr_nu01_d   = dust_params['decorr_nu01_d']
decorr_nu02_d   = dust_params['decorr_nu02_d']

decorr_amp_s    = sync_params['decorr_amp_s']
decorr_nu01_s   = sync_params['decorr_nu01_s']
decorr_nu02_s   = sync_params['decorr_nu02_s']

PATH_DICT = dict(pathnames())

LMIN             = config.bpw_param.lmin
DELL             = config.bpw_param.dell
NBANDS           = config.bpw_param.nbands

LMAX, LARR_ALL, LBANDS, LEFF = get_ell_arrays(LMIN, DELL, NBANDS)

NSIDE            = config.global_param.nside
EXPERIMENT       = config.global_param.experiment
POLARIZATION     = config.global_param.polarization

band_names_config = config.band_names
