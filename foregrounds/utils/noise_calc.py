'''
noise_calc 

calculates noise 
'''

import numpy as np

def S4_beam_FWHM():

    '''
    returns the S4 SAT beam in arcminutes as a function of frequency
    https://arxiv.org/pdf/2008.12619.pdf Table 1
    N.B.: Table 1 Namikawa & Sherwin 2023 only used 1 frequency (20' beam)
    '''

    beam_SAT_30  = 77.
    beam_SAT_40  = 58.
    beam_SAT_85  = 27.
    beam_SAT_95  = 24.
    beam_SAT_145 = 16.
    beam_SAT_155 = 15.
    beam_SAT_220 = 11.
    beam_SAT_270 = 8.5

    return np.array([beam_SAT_30, beam_SAT_40, beam_SAT_85, beam_SAT_95, beam_SAT_145, beam_SAT_155, beam_SAT_220, beam_SAT_270])


def S4_SAT_beams(ell):
    SA_beams = S4_beam_FWHM() / np.sqrt(8. * np.log(2)) /60. * np.pi/180.
    ## SAT beams as a sigma expressed in radians
    return [np.exp(-0.5*ell*(ell+1)*sig**2.) for sig in SA_beams]

def Simons_Observatory_V3_SA_beam_FWHM():
    ## returns the SAT beams in arcminutes
    beam_SAT_27 = 91.
    beam_SAT_39 = 63.
    beam_SAT_93 = 30.
    beam_SAT_145 = 17.
    beam_SAT_225 = 11.
    beam_SAT_280 = 9.
    return(np.array([beam_SAT_27,beam_SAT_39,beam_SAT_93,beam_SAT_145,beam_SAT_225,beam_SAT_280]))

def Simons_Observatory_V3_SA_beams(ell):
    SA_beams = Simons_Observatory_V3_SA_beam_FWHM() / np.sqrt(8. * np.log(2)) /60. * np.pi/180.
    ## SAT beams as a sigma expressed in radians
    return [np.exp(-0.5*ell*(ell+1)*sig**2.) for sig in SA_beams]