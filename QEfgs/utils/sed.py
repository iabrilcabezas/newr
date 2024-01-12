'''
sed
    computes component spectra and SEDs for dust
'''

import numpy as np

dust_params = {'A_dust_BB': 28., # from BICEP field 1510.09217,
               'EB_dust': 2., # PL amplitude 1801.04945,
               'alpha_dust_EE': -0.42, # PL exponent 1801.04945
               'alpha_dust_BB': -0.54, # PL exponent 1801.04945
               'beta_dust': 1.59, # modified BB emission 1409.5738
               'temp_dust': 19.6, # modified BB emission 1409.5738
               'nu0_dust': 353.}


def fcmb(nu):

    '''
    spectral energy density in CMB units
    '''

    x = 0.017608676067552197*nu
    ex = np.exp(x)
    return ex*(x/(np.expm1(x)))**2

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

fdust = lambda x: comp_sed(x, dust_params['nu0_dust'], dust_params['beta_dust'], dust_params['temp_dust'], 'dust')