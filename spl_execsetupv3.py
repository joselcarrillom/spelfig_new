'''
This is the script containing the functions to set up properly the data for the spelfig initial guess fitting
methods.
'''
# Import modules
import os
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


from scipy.signal import find_peaks
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

# Create spectrum array
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

def analyze_emission_lines(x, y, lines_dict, window=20.):

    observed_wavelengths = x
    flux = y

    matched_lines = []
    synthetic_flux = np.copy(flux)  # This maintains the original 1D flux array structure

    # Ensure the synthetic_flux initialization doesn't start with NaN values
    if np.isnan(synthetic_flux).all():
        synthetic_flux[:] = 0  # Set to zero or some baseline if entirely NaN

    for name, params in lines_dict.items():
        rest_wavelength = params['wavelength'][0]
        if rest_wavelength is not None:
            lower_bound = rest_wavelength - window
            upper_bound = rest_wavelength + window
            mask = (observed_wavelengths >= lower_bound) & (observed_wavelengths <= upper_bound)
            window_flux = flux[mask]
            window_wavelengths = observed_wavelengths[mask]

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
                peak_flux = window_flux[closest_peak_idx]
                observed_wavelength = window_wavelengths[closest_peak_idx]
                matched_lines.append({
                    'line': name,
                    'restframe_wavelength': rest_wavelength,
                    'observed_wavelength': observed_wavelength,
                    'peak_flux': peak_flux,
                    'peak_idx': np.where(observed_wavelengths == observed_wavelength)[0][0]
                })

    return matched_lines

# Signal-to-noise ratio
def extract_redshift_snr(spectrum, redshift, line=None):

    """
    This function calculates the S-N ratio of a single emission line in a spectra

    Args:
        spectrum (np.ndarray): A 2D NumPy array containing the observed spectrum data.
        redshift (float): The calculated redshift of the observed spectrum.
        line (str, optional): A line to extract SNR from, based on the dictionary SNR_emission_lines. Default is O-III.

    Returns:
        float: The S-N ratio of the observed spectrum.

    """

    wavelength = spectrum[:, 0]/(1+redshift)
    flux = (spectrum[:, 1])
    norm_flux = flux/max(flux)
    errors = spectrum[:, 2]

    # Select which line from the emission lines to use
    rest_wavelength = SNR_emission_lines[line]['wavelength'][0] if line is not None else SNR_emission_lines['O-III,0']['wavelength'][0]

    # peak_flux = emission_identify(rest_wavelength, wavelength, flux)

    # Define continuum

    continuum_mask = (wavelength >= 5050) & (wavelength <= 5150)
    continuum_data = norm_flux[continuum_mask]
    continuum_wavelengths = wavelength[continuum_mask]

    continuum_fit = np.polyfit(continuum_wavelengths, continuum_data, deg=1)
    continuum_model = np.polyval(continuum_fit, continuum_wavelengths)

    continuum_snr = np.std(continuum_model)
    continuum_mean = np.mean(continuum_model)

    # Define SNR
    # snr = (peak_flux - continuum_mean) / continuum_snr if continuum_snr > 0 else print('')

    matched_lines = analyze_emission_lines(wavelength, norm_flux, SNR_emission_lines)

    return continuum_snr, matched_lines

# Initialise Templates
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
            print("spec:", spect0)

        except KeyError:
            print("Please ensure correct keys are used to call headers (i.e. 'WAVELENGTH' & 'FLUX').")

    return template_spectra

# Calculate redshift
def redshift_calcu(spectrum, template_spectra, filename, Plot=False, user_template=None):
    """
    This function calculates the redshift of an observed spectrum by comparing it to a set of template spectra.

    Args:
        spectrum (np.ndarray): A 2D NumPy array containing the observed spectrum data.
            - spectrum[:, 0]: Wavelength (Angstroms)
            - spectrum[:, 1]: Flux (Jy)
            - spectrum[:, 2]: Errors in flux (Jy)
        template_spectra (list): A list of `specutils.Spectrum1D` objects representing the template spectra.
        Plot (bool, optional): If True, creates a plot comparing the observed and redshifted spectra. Default is False.
        user_template (str, optional): Name of a specific template to use. Overrides other selection methods.

    Returns:
        float: The estimated redshift of the observed spectrum.
    """

    wavelength = spectrum[:, 0]
    fluxes = spectrum[:, 1]
    errors = spectrum[:, 2]

    spec_data = Spectrum1D(flux=fluxes * u.Jy, spectral_axis=wavelength * u.Angstrom,
                           uncertainty=StdDevUncertainty(errors), unit='Jy')
    redshift_range = np.arange(0.00, 2.49, 0.01)

    print(filename)

    if user_template is None:
        selected_templates = []
        if 'manga' in filename.lower() or 'magpi' in filename.lower():
            selected_templates = [template_spectra[0]]
        elif 'paqs' in filename.lower():
            selected_templates = [template_spectra[1]]
        if not selected_templates:
            selected_templates = [template_spectra[0]]  # Default to first templates

    # Template selection based on user input (overrides filename keywords)
    elif user_template:
        selected_templates = [t for t in template_spectra if t.meta['name'] == user_template]
        if not selected_templates:
            raise ValueError(f"Template '{user_template}' not found in provided spectra")

    _, redshiftspec, _, _, _ = template_match(spec_data, selected_templates, redshift=redshift_range)

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