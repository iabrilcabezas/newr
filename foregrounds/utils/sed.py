'''
sed
    computes component spectra and SEDs
'''

import numpy as np
from utils.params import EXPERIMENT, PATH_DICT
from utils.params import band_names_config
from utils.params import amp_d_bb, alpha_d_bb
from utils.params import beta_d, temp_d, nu0_d
from utils.params import A_lens, l0_s_bb, l0_d_bb
from utils.params import amp_s_bb, alpha_s_bb
from utils.params import beta_s, nu0_s



def get_band_names():

    '''
    Define name of frequency bands according to experiment

    ** Parameters**
    experiment: str
                's4_sat'

    **Returns**
    list of band names
    '''

    if EXPERIMENT == 'S4_SAT':
        band_names = band_names_config.s4sat
    if EXPERIMENT == 'so':
        band_names = band_names_config.so
    if EXPERIMENT == 'lbrd':
        band_names = band_names_config.lbrd

    return band_names

def s4_frequencies():
    '''
    returns S4 SAT frequency channels in GHz
    from https://arxiv.org/pdf/2008.12619.pdf, section 2.3
    '''

    s4_channels = [30,40, 85, 95, 145, 155, 220, 270]

    return s4_channels

def fcmb(nu):

    '''
    spectral energy density in CMB units
    '''

    x = 0.017608676067552197*nu
    ex = np.exp(x)
    return ex*(x/(np.expm1(x)))**2

#Component power spectra
def dl_plaw(A,alpha,ls, lnorm_PL):
    '''power law of index alpha'''
    return A*((ls+0.001)/lnorm_PL)**alpha

#All spectra
def comp_sed(nu,nu0,beta,temp,typ):

    '''
    SED for CMB and dust components
    '''

    if typ == 'cmb':
        return fcmb(nu)
    if typ == 'dust':
        x_to=0.04799244662211351*nu/temp
        x_from=0.04799244662211351*nu0/temp
        return (nu/nu0)**(1+beta)*np.expm1(x_from)/np.expm1(x_to)*fcmb(nu0)
    elif typ == 'sync':
        return (nu/nu0)**beta*fcmb(nu0)
    return None

def read_camb(fname, lmax):
    ''''reads camb data for CMB P(k)'''
    larr_all = np.arange(lmax+1)
    l,dtt,dee,dbb,dte = np.loadtxt(fname,unpack=True)
    l = l.astype(int)
    msk = l <= lmax
    l = l[msk]
    dltt = np.zeros(len(larr_all))
    dltt[l] = dtt[msk]
    dlee = np.zeros(len(larr_all))
    dlee[l] = dee[msk]
    dlbb = np.zeros(len(larr_all))
    dlbb[l] = dbb[msk]
    dlte = np.zeros(len(larr_all))
    dlte[l] = dte[msk]
    return dltt,dlee,dlbb,dlte



class Bpass_delta:
    '''bandpass class object'''
    def __init__(self, freq):
        self.name = str(freq)
        self.nu = freq
        self.bnu = 1.
        # CMB units
        norm = self.nu**2*fcmb(self.nu)
        self.bnu /= norm

    def convolve_sed(self,f):
        '''convolves SED to get normalization'''
        sed = self.bnu*self.nu**2*f(self.nu)
        return sed

#Bandpasses
class Bpass_band:
    '''bandpass class object'''
    def __init__(self,name,fname):
        self.name = name
        self.nu,self.bnu = np.loadtxt(fname, usecols = (0,1), unpack= True)
        self.dnu = np.zeros_like(self.nu)
        self.dnu[1:] = np.diff(self.nu)
        self.dnu[0] = self.dnu[1]
        # CMB units
        norm = np.sum(self.dnu*self.bnu*self.nu**2*fcmb(self.nu))
        self.bnu /= norm

    def convolve_sed(self,f):
        '''convolves SED to get normalization'''
        sed = np.sum(self.dnu*self.bnu*self.nu**2*f(self.nu))
        return sed

def get_component_spectra(lmax):
    '''gets component spectra of CMB and dust, E and B polarization'''

    larr_all = np.arange(lmax+1)

    dls_dust_bb=dl_plaw(amp_d_bb,alpha_d_bb,larr_all,l0_d_bb)

    dls_sync_bb=dl_plaw(amp_s_bb,alpha_s_bb,larr_all,l0_s_bb)

    _,_,dls_cmb_bb,_=read_camb( PATH_DICT['bbpower_path'] + 'examples/data/camb_lens_nobb.dat', lmax)

    return (dls_dust_bb,
            A_lens*dls_cmb_bb,
            dls_sync_bb)

def get_convolved_seds(names, bpss):
    '''convolves SED with bandpasses'''
    nfreqs = len(names)
    seds = np.zeros([3,nfreqs])
    for ib, n in enumerate(names):
        b = bpss[n]
        seds[0,ib] = b.convolve_sed(lambda nu : comp_sed(nu,None,None,None,'cmb'))
        seds[1,ib] = b.convolve_sed(lambda nu : comp_sed(nu,nu0_d,beta_d,temp_d,'dust'))
        seds[2,ib] = b.convolve_sed(lambda nu : comp_sed(nu,nu0_s,beta_s,None,'sync'))
    return seds
