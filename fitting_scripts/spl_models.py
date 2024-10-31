from astropy.modeling.models import Voigt1D
import numpy as np

# Analytical single distributions:
def gauss(x, loc, a0, sd):
    diff = x-loc
    return  a0 * np.exp(-0.5*(((diff) / sd)**2))

def asym_gauss(x, loc, a0, A, w):
    diff = x-loc
    sigma_assym = A*diff + w
    return a0 * np.exp(-0.5*(((diff) / sigma_assym)**2))

def lorentzian(x, loc, a0, gamma):
    diff = x-loc
    return  a0 * (1/(1+((diff) / gamma)**2))

def voigt(x, loc, a0, sd, gamma):
    fwhm_L = 2*gamma
    fwhm_G = 2*np.sqrt(2*np.log(2))*sd
    return Voigt1D(x_0=loc, amplitude_L=a0, fwhm_L=fwhm_L, fwhm_G=fwhm_G)(x)

def continuum_function(x, a, b, c):
    if b == 0:
        a == 0
    return a*(x/b)**(-c)

# Velocity conversion:
def vel_correct(l0, l):
    '''
    Transform wavelength to units of velocity (km/s) given a reference wavelength l0 (in case of a single gaussian-fitted-line analysis, it is the centroid).
    '''
    c = 299792.458
    vel = c*((l-l0)/l0)
    return vel
