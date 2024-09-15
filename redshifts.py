# Import modules
from scipy.signal import find_peaks
import os

import numpy as np
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt

from specutils import Spectrum1D
from specutils.analysis import template_match
from astropy.nddata import StdDevUncertainty

# Emission lines to compare
SNR_emission_lines = {
    'He-II,1': {'wavelength':[3202.15]},
    'He-II,2': {'wavelength':[4685.74]},
    'Ne-V,1':  {'wavelength':[3345.81]},
    'Ne-V,2':  {'wavelength':[3425.81]},
    'O-II,1':  {'wavelength':[3726.03]},
    'O-II,2':  {'wavelength':[3728.73]},
    'Ne-III,1': {'wavelength':[3868.69]},
    'Ne-III,2': {'wavelength':[3967.40]},
    'H-ζ':    {'wavelength':[3889.05]},
    'H-ε':    {'wavelength':[3970.07]},
    'H-δ':    {'wavelength':[4101.73]},
    'H-γ':    {'wavelength':[4340.46]},
    'O-III,0':  {'wavelength':[4363.15]},
    'O-III,1':  {'wavelength':[4958.83]},
    'O-III,2':  {'wavelength':[5006.77]},
    'Ar-IV,1':  {'wavelength':[4711.30]},
    'Ar-IV,2':  {'wavelength':[4740.10]},
    'H-β':     {'wavelength':[4861.32]},
    'N-I,1':    {'wavelength':[5197.90]},
    'N-I,2':    {'wavelength':[5200.39]},
    'He-I':   {'wavelength':[5875.60]},
    'O-I,1':   {'wavelength':[6300.20]},
    'O-I,2':   {'wavelength':[6363.67]},
    'N-II,1':   {'wavelength':[6547.96]},
    'N-II,2':   {'wavelength':[6583.34]},
    'H-α':     {'wavelength':[6562.80]},
    'S-II,1':  {'wavelength': [6716.31]},
    'S-II,2':   {'wavelength':[6730.68]},
    'Ar-III': {'wavelength':[7135.67]},
    }

# Initialise template data
def initialise_templates(dir):
    '''
    This function loads template spectra from a specified directory and converts them to `specutils.Spectrum1D` objects.
    The original data used for this script was taken from https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/agn/

    Args:
        dir (str): The directory containing the template FITS files.
    Returns:
        list: A list of `specutils.Spectrum1D` objects representing the loaded template spectra.
    Raises:
        KeyError: If the expected header keywords ("WAVELENGTH" and "FLUX") are not found in the FITS files.
    '''
    template_spectra = []
    template_filenames = [f for f in os.listdir(dir) if f.endswith('.fits')]

    for filename in template_filenames:
        full_dir = os.path.join(dir, filename)
        template_data = fits.open(full_dir)

        try:
            wavelength = template_data[1].data['WAVELENGTH'] * u.Angstrom
            flux = template_data[1].data['FLUX'] * u.Jy

            spect0 = Spectrum1D(flux=flux, spectral_axis=wavelength)
            template_spectra.append(spect0)

            return template_spectra

        except KeyError:
            print("Please ensure correct keys are used to call headers.")

# Redshift function
def redshift_calc(spectrum, template_spectra, Plot=False):
    """
    This function calculates the redshift of an observed spectrum by comparing it to a set of template spectra.

    Args:
        spectrum (np.ndarray): A 2D NumPy array containing the observed spectrum data.
            - spectrum[:, 0]: Wavelength (Angstroms)
            - spectrum[:, 1]: Flux (Jy)
            - spectrum[:, 2]: Errors in flux (Jy)
        template_spectra (list): A list of `specutils.Spectrum1D` objects representing the template spectra.
        Plot (bool, optional): If True, creates a plot comparing the observed and redshifted spectra. Default is False.

    Returns:
        float: The estimated redshift of the observed spectrum.
    """
    wavelength = spectrum[:, 0]
    fluxes = spectrum[:, 1]
    errors = spectrum[:, 2]

    spec_data = Spectrum1D(flux=fluxes * u.Jy, spectral_axis=wavelength * u.Angstrom,
                               uncertainty=StdDevUncertainty(errors), unit='Jy')
    redshift_range = np.arange(0.05, 2.49, 0.01)

    _, redshiftspec, _, _, _ = template_match(spec_data, template_spectra, redshift=redshift_range)

    if Plot == True:
        plt.plot(wavelength, fluxes, label='observed', color='red', alpha=0.5)
        plt.plot(wavelength / (1 + redshiftspec), fluxes, label='rest', color='blue', alpha=1)

        # Plot vertical lines for emission lines
        for line_name, line_data in SNR_emission_lines.items():
            for rest_wavelength in line_data['wavelength']:
                plt.vlines(rest_wavelength, ymin=min(fluxes), ymax=max(fluxes), color='black', linestyle='-.',
                           alpha=0.8)

        plt.xlabel('Wavelength')
        plt.ylabel('Flux')

        plt.legend()
        plt.show()

    return redshiftspec