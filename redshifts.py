# Import modules
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_widths
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum1D
from specutils.analysis import template_match


# Flatten data
def extract_data(array):

    """
    Flattens multidimensional data arrays if required

    Args:
        array (np.ndarray): The input NumPy array.

    Returns:
        np.ndarray: A flattened 1D NumPy array containing the extracted data.

    Raises:
        ValueError: If the array has more than two dimensions.
    """

    if array.ndim == 1:
        return array
    elif array.ndim == 2:
        return array.flatten()
    else:
        raise ValueError("Unsupported array dimension: {}".format(array.ndim))

#
def extract_astronomical_data(filename, verbose=True, wavelength_name=None, flux_name=None, error_name=None):

    """
    Extracts wavelength, flux, and error (when applicable) from a FITS file, handling variations
    in keyword syntax and layout across surveys.

    Add your own keywords if it doesn't work.

    Args:
        filename (str): The input FITS file.
        verbose (bool, optional): If True, prints out the extracted data.
        wavelength_name (str, optional): The name of the wavelength column in the input FITS file.
        flux_name (str, optional): The name of the flux column in the input FITS file.
        error_name (str, optional): The name of the error column in the input FITS file.

    """

    try:
        with fits.open(filename) as hdul:

            # Search for keywords
            potential_keywords = ['wavelength', 'WAVE', 'lambda']
            for keyword in potential_keywords:
                if keyword in hdul[1].columns.names:
                    wavelength_name = keyword
                    break

            potential_keywords = ['flux', 'FLUX', 'f_lambda', 'flux_lines','total flux']
            for keyword in potential_keywords:
                if keyword in hdul[1].columns.names:
                    flux_name = keyword
                    break

            potential_keywords = ['error', 'ERR', 'flux_error', 'ERR_FLUX','flux error']
            for keyword in potential_keywords:
                if keyword in hdul[1].columns.names:
                    error_name = keyword
                    break

            # Define arrays

            if wavelength_name:                             # wavelength
                wave = hdul[1].data[wavelength_name]
            else:
                if verbose:
                    print("Warning: Keyword for wavelength not found.")
                pass

            if flux_name:                                   # flux
                flux = hdul[1].data[flux_name]
            else:
                if verbose:
                    print("Warning: Keyword for flux not found.")
                pass

            if error_name:                                  # errors
                err = hdul[1].data[error_name]
            else:
                if verbose:
                    print("Warning: Keyword for error not found.")
                pass

            # Create spectra with correct form

            wave = extract_data(wave)
            flux = extract_data(flux)
            err = extract_data(err)

            spec = np.array([wave, flux, err]).T

            return spec

    except Exception as e:
        print(f"Error processing {filename}: {e}")
        return None

# Emission lines to compare

SNR_emission_lines = {
    'He-II,1':  {'wavelength':[3202.15]},
    'He-II,2':  {'wavelength':[4685.74]},
    'Ne-V,1':   {'wavelength':[3345.81]},
    'Ne-V,2':   {'wavelength':[3425.81]},
    'O-II,1':   {'wavelength':[3726.03]},
    'O-II,2':   {'wavelength':[3728.73]},
    'Ne-III,1': {'wavelength':[3868.69]},
    'Ne-III,2': {'wavelength':[3967.40]},
    'H-ζ':      {'wavelength':[3889.05]},
    'H-ε':      {'wavelength':[3970.07]},
    'H-δ':      {'wavelength':[4101.73]},
    'H-γ':      {'wavelength':[4340.46]},
    'O-III,0':  {'wavelength':[4363.15]},
    'O-III,1':  {'wavelength':[4958.83]},
    'O-III,2':  {'wavelength':[5006.77]},
    'Ar-IV,1':  {'wavelength':[4711.30]},
    'Ar-IV,2':  {'wavelength':[4740.10]},
    'H-β':      {'wavelength':[4861.32]},
    'N-I,1':    {'wavelength':[5197.90]},
    'N-I,2':    {'wavelength':[5200.39]},
    'He-I':     {'wavelength':[5875.60]},
    'O-I,1':    {'wavelength':[6300.20]},
    'O-I,2':    {'wavelength':[6363.67]},
    'N-II,1':   {'wavelength':[6547.96]},
    'N-II,2':   {'wavelength':[6583.34]},
    'H-α':      {'wavelength':[6562.80]},
    'S-II,1':   {'wavelength':[6716.31]},
    'S-II,2':   {'wavelength':[6730.68]},
    'Ar-III':   {'wavelength':[7135.67]},
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
            print("Please ensure correct keys are used to call headers (i.e. 'WAVELENGTH' & 'FLUX').")

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
    redshift_range = np.arange(0.00, 2.49, 0.01)

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

def extract_snr(spectrum, redshift, line=None):

  wavelength = spectrum[:, 0]/(1+redshift)
  flux = spectrum[:, 1]
  errors = spectrum[:, 2]

  # Select which line from the emission lines to use
  rest_wavelength = SNR_emission_lines[line]['wavelength'][0] if line is not None else SNR_emission_lines['O-III,0']['wavelength'][0]

  # Identify the correct peaks
  lower_bound = rest_wavelength - 30
  upper_bound = rest_wavelength + 30
  mask = (wavelength >= lower_bound) & (wavelength <= upper_bound)
  window_flux = flux[mask]
  window_wavelengths = wavelength[mask]

  peaks, _ = find_peaks(window_flux, prominence=0.5)
  closest_peak_idx = None
  min_diff = float('inf')

  for peak_idx in peaks:

    observed_wavelength = window_wavelengths[peak_idx]
    diff = abs(observed_wavelength - rest_wavelength)

    if diff < min_diff:
        min_diff = diff
        closest_peak_idx = peak_idx

  if closest_peak_idx is not None:
      peak_flux = window_flux[peak_idx]

  else:
      print("Closest peak index not found for ", line)

  # Define continuum

  continuum_mask = (wavelength >= lower_bound) & (wavelength <= upper_bound)
  continuum_data = flux[continuum_mask]
  continuum_wavelengths = wavelength[continuum_mask]

  continuum_fit = np.polyfit(continuum_wavelengths, continuum_data, deg=1)
  continuum_model = np.polyval(continuum_fit, continuum_wavelengths)

  local_std = np.std(continuum_model)
  local_mean = np.mean(continuum_model)

  # Define SNR
  snr = (peak_flux - local_mean) / local_std if local_std > 0 else print('')

  return snr