'''
This is the script containing the functions to set up properly the data for the spelfig fitting
methods.
'''
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths
from scipy.optimize import curve_fit

import spl_models as spm

# Set up the lines initial dataframes:
def analyze_emission_lines(x, y, lines_dict, window=10.):
    '''
    Identify emission lines in a given observed spectrum and return the results.
    :param spectrum:
    :param lines_dict:
    :param window:
    :return:
    '''
    observed_wavelengths = x
    flux = y

    matched_lines = []
    results = []
    synthetic_flux = np.copy(flux)  # This maintains the original 1D flux array structure

    # Ensure the synthetic_flux initialization doesn't start with NaN values
    if np.isnan(synthetic_flux).all():
        synthetic_flux[:] = 0  # Set to zero or some baseline if entirely NaN

    continuum_mask = np.ones(len(flux), dtype=bool)

    for name, params in lines_dict.items():
        rest_wavelength = params['wavelength'][0]
        print('This name', name)
        print('This rest_wavelength', rest_wavelength)
        if rest_wavelength is not None:
            print('This rest_wavelength', rest_wavelength)
            lower_bound = rest_wavelength - window
            upper_bound = rest_wavelength + window
            mask = (observed_wavelengths >= lower_bound) & (observed_wavelengths <= upper_bound)
            window_flux = flux[mask]
            window_wavelengths = observed_wavelengths[mask]

            if not window_wavelengths.size:
                continue

            peaks, _ = find_peaks(window_flux, prominence=0.5)
            closest_peak_idx = None
            min_diff = float('inf')

            for peak_idx in peaks:
                ## This lines were meant to identify the closest peak, but they might screw up
                # the things if the redshift of the spectrum is not quite correct.

                observed_wavelength = window_wavelengths[peak_idx]
                diff = abs(observed_wavelength - rest_wavelength)

                if diff < min_diff:
                    min_diff = diff
                    closest_peak_idx = peak_idx

            if closest_peak_idx is not None:
                peak_flux = window_flux[peak_idx]
                observed_wavelength = window_wavelengths[peak_idx]
                matched_lines.append({
                    'line': name,
                    'restframe_wavelength': rest_wavelength,
                    'observed_wavelength': observed_wavelength,
                    'peak_flux': peak_flux,
                    'peak_idx': np.where(observed_wavelengths == observed_wavelength)[0][0]
                })

    # Update continuum mask and calculate synthetic data for gaps
    for line in matched_lines:
        peak_idx = line['peak_idx']
        widths, width_heights, left_ips, right_ips = peak_widths(flux, [peak_idx], rel_height=0.5)
        fwhm = np.interp(left_ips + widths, np.arange(len(flux)), observed_wavelengths) - \
               np.interp(left_ips, np.arange(len(flux)), observed_wavelengths)

        sigma = fwhm[0] / 2.355
        line_start = line['observed_wavelength'] - fwhm[0] / 2 - 3 * sigma
        line_end = line['observed_wavelength'] + fwhm[0] / 2 + 3 * sigma

        results.append({
            'line': line['line'],
            'restframe_wavelength': line['restframe_wavelength'],
            'observed_wavelength': line['observed_wavelength'],
            'fwhm': fwhm[0],
            'peak_flux': line['peak_flux']
        })

        mask = ((observed_wavelengths >= line_start) & (observed_wavelengths <= line_end))
        continuum_mask &= ~mask
        if np.any(mask):
            valid_flux = flux[~mask]
            if valid_flux.size > 0 and not np.isnan(valid_flux).all():
                mean_flux = np.nanmean(valid_flux)
                std_flux = np.nanstd(valid_flux)
                synthetic_flux[mask] = np.random.normal(0.0, 0.1 * std_flux, np.sum(mask))
            else:
                synthetic_flux[mask] = 0  # Fallback if no valid data is available

    # Apply the mask to keep only continuum points
    continuum_spectrum = np.column_stack((observed_wavelengths, synthetic_flux))
    std_cont = np.nanstd(continuum_spectrum[:,1])
    return results, std_cont, continuum_spectrum


def filter_and_prepare_linelist(line_results, continuum_spec0, wavelength_range, snr_ext,
                                window_width=10.):
    min_wavelength, max_wavelength = wavelength_range
    filtered_linelist = []
    # If noise standard deviation is not provided, assume a default or calculate externally

    three_sigma = 3 * snr_ext  # 3 sigma threshold for noise
    print(f'This is 3 sigma: {three_sigma}')

    for line in line_results:
        name = line['line']
        line_wavelength = line['observed_wavelength']
        line_fwhm = line['fwhm']
        line_flux = line['peak_flux']
        sigma_max = line_fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))  # Convert FWHM to sigma

        # Filter based on the wavelength range and FWHM constraints
        if min_wavelength <= line_wavelength <= max_wavelength:
            lineloc_min = line_wavelength - 2 * sigma_max
            lineloc_max = line_wavelength + 2 * sigma_max

            window_left = continuum_spec0[(continuum_spec0[:, 0] <= lineloc_min) &
                                          (continuum_spec0[:, 0] >= lineloc_min - window_width)]
            window_right = continuum_spec0[(continuum_spec0[:, 0] <= lineloc_max + window_width) &
                                           (continuum_spec0[:, 0] >= lineloc_max)]

            local_std = np.mean([np.std(window_left[:, 1]), np.std(window_right[:, 1])])
            local_mean = np.mean([np.mean(window_left[:, 1]), np.mean(window_right[:, 1])])

            # Evaluate inf snr in the line region is good regardless of entire spectrum
            snr = (line_flux - local_mean) / local_std if local_std > 0 else snr_ext

            # Check if the SNR is above the threshold and the peak flux is significant above 3 sigma
            if snr >= snr_ext or line_flux > three_sigma:
                line_details = {
                    'name': name,
                    'wavelength': line_wavelength,
                    'sigma': sigma_max,
                    'min_loc': lineloc_min,
                    'max_loc': lineloc_max,
                    'min_sd': 2.0,
                    'max_sd': 1.5*sigma_max,
                    'max_flux': line_flux,
                    'SNR': snr,  # Including SNR in the output for reference
                }
                filtered_linelist.append(line_details)
    return filtered_linelist

def continuum_init(continuum_spec, g_init):
    '''
    Initial guess for the parameters of the continuum. It uses a first approach fit with
    :param continuum_spec:
    :return:
    '''
    x_continuum = continuum_spec[:, 0]
    y_continuum = continuum_spec[:, 1]
    a_init = np.mean(y_continuum)
    loc0_init = np.min(x_continuum)
    p0 = [a_init, loc0_init, g_init]
    params, params_covariance = curve_fit(spm.continuum_function, x_continuum, y_continuum, p0=p0)

    return params


def initial_dataframe(emlines_dict, filtered_linelist, continuum_pars=None):
    '''
    This function creates an initial parameters dataframe given the emission lines dictionary.
    :param emlines_dict: The dictionary of emission lines
    :return:
    '''
    # Initialize lists to hold data for DataFrame
    line_names = []
    models = []
    ncomp = []
    parameters = []
    min_limits = []
    max_limits = []

    emlines = emlines_dict.keys()
    for line in filtered_linelist:
        line_name = line['name']
        components = emlines_dict[line_name]['components']
        Ncomp = len(components)
        for j, component in enumerate(components):
            if component == 'Voigt':
                params_i = [line['wavelength'], line['max_flux']/Ncomp, line['sigma'],
                            1.11*line['sigma']]
                max_i = [line['max_loc'], line['max_flux'], line['max_sd'], 1.11*line['max_sd']]
                min_i = [line['min_loc'], 0., line['min_sd'], 1.11*line['min_sd']]
            elif component == 'Gaussian':
                params_i = [line['wavelength'], line['max_flux']/Ncomp, line['sigma']]
                max_i = [line['max_loc'], line['max_flux'], line['max_sd']]
                min_i = [line['min_loc'], 0., line['min_sd']]
            elif component == 'Lorentzian':
                params_i = [line['wavelength'], line['max_flux']/Ncomp, 1.11*line['sigma']]
                max_i = [line['max_loc'], line['max_flux'], 1.11*line['max_sd']]
                min_i = [line['min_loc'], 0., 1.11*line['min_sd']]

            line_names.append(line_name)
            models.append(component)  # Assuming wavelength as centroid
            ncomp.append(j+1)  # Using max_flux as initial guess for amplitude
            parameters.append(params_i)
            min_limits.append(min_i)
            max_limits.append(max_i)

    if continuum_pars is not None:
        print('Appending continuum')
        print(continuum_pars)
        parameters.append(continuum_pars)
        line_names.append('Continuum')
        models.append('Continuum')
        ncomp.append(0)
        min_limits.append([-np.inf, 0, -np.inf])
        max_limits.append([np.inf, np.inf, np.inf])

    dfparams = pd.DataFrame({
        'Line Name': line_names,
        'Model': models,  # Initial components set to 1
        'Component': ncomp,
        'Parameters': parameters,
        'Max Limits': max_limits,  # Using max_flux as initial guess for amplitude': sigmas,
        'Min Limits': min_limits,
    })

    return dfparams


def init_setup(spectrum, emlines_dict, wavelength_range, gamma_init):
    '''
    This function sets up all the objects necessary for the execution of the fit functions
    :param spectrum:
    :param lines_dict:
    :param gamma_init: Initial guess for
    '''

    spectrum = spectrum[spectrum[:,0]>=wavelength_range[0] & spectrum[:,0]<=wavelength_range[1]]
    x = spectrum[:, 0]
    y = spectrum[:, 1]
    dy = spectrum[:, 2]
    
    # First restrict the spectrum to the wavelength range:

    # Find the lines present in the spectrum and estimate first guesses for parameters:
    lines_init, snr_cont, continuum0 = analyze_emission_lines(x, y, emlines_dict)

    # Initial guess for the parameters of the continuum:
    continuum_pars = continuum_init(continuum0, gamma_init)

    # Filter the list of lines present in the spectrum:
    linelist0 = filter_and_prepare_linelist(lines_init, continuum0, wavelength_range, snr_cont)

    # Create an initial parameters dataframe:
    dfparams = initial_dataframe(emlines_dict, linelist0, continuum_pars=continuum_pars)

    return dfparams

















# Plotting function:   --------------

def spl_plot(x, y, dy, dfparams, x_zoom=None, y_zoom=None, goodness_marks=None):
    # Extract data from the DataFrame

    x_fit = np.linspace(min(x) - 10, max(x) + 10, 10000)
    y_fit = np.zeros_like(x_fit)
    y_evaluated = np.zeros_like(x)

    # Create a figure with two subplots (upper and lower panels)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5), sharex=True,
                                   gridspec_kw={'hspace': 0.0, 'height_ratios': [5, 2]})

    # Upper panel: Observed spectrum, model, and individual Gaussian components
    for _, row in dfparams.iterrows():
        line_name = row['Line Name']
        component = int(row['Component'])
        model = row['Model']

        if model == 'Continuum':
            component_y = spm.continuum_function(x_fit, *row['Parameters'])
            component_y_ev = spm.continuum_function(x, *row['Parameters'])
        if model == 'Gaussian':
            component_y = spm.gauss(x_fit, *row['Parameters'])
            component_y_ev = spm.gauss(x, *row['Parameters'])
        elif model == 'Lorentzian':
            component_y = spm.lorentzian(x_fit, *row['Parameters'])
            component_y_ev = spm.lorentzian(x, *row['Parameters'])
        elif model == 'Voigt':
            component_y = spm.voigt(x_fit, *row['Parameters'])
            component_y_ev = spm.voigt(x, *row['Parameters'])
        elif model == 'Asym_Gauss':
            component_y = spm.asym_gauss(x_fit, *row['Parameters'])
            component_y_ev = spm.asym_gauss(x, *row['Parameters'])

        # Plot individually the component:

        color = (random.random(), random.random(), random.random())
        ax1.plot(x_fit, component_y, linestyle='--', linewidth=0.8, color=color)

        y_fit += component_y
        y_evaluated += component_y_ev

    ax1.plot(x_fit, y_fit, color='crimson', linewidth=2.0, label='Total Fitted Spectrum')
    ax1.errorbar(x, y, yerr=dy, color='grey', linestyle='-',
                 marker='.', alpha=0.7, markersize=2, linewidth=0.7, label='Observed Spectrum')

    # Set the y-axis label and legend for the upper panel
    ax1.set_ylabel(r'I [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')
    ax1.legend(loc='upper right', frameon=False)

    # Add Chi-squared and BIC as labels with transparency:
    if goodness_marks:
        chi2 = goodness_marks['Reduced Chi-Squared']
        BIC = goodness_marks['BIC']
        ax1.text(0.05, 0.85, f'Chi-squared: {chi2:.2f}', transform=ax1.transAxes, fontsize=12,
                 color='gray', alpha=0.8)
        ax1.text(0.05, 0.78, f'BIC: {BIC:.2f}', transform=ax1.transAxes, fontsize=12,
                 color='gray', alpha=0.8)

    residuals = (abs(y - y_evaluated) / y) * 100  # Compute percentage residuals
    e_range = [-0.5, 110]
    ax2.plot(x, residuals, marker='.', color='steelblue', linestyle='--', alpha=0.5,
             markersize=2, label='Percentage error')
    ax2.set_ylim(e_range)

    # Set the x-axis and y-axis labels for the lower panel
    ax2.set_xlabel(r'$\lambda$ [$\AA$]')
    ax2.set_ylabel('Residuals')
    ax2.legend(loc='upper right', frameon=False)

    # Set the zoom range if provided
    if x_zoom:
        ax1.set_xlim(x_zoom)
        ax2.set_xlim(x_zoom)
    if y_zoom:
        ax1.set_ylim(y_zoom)

    # Fine-tune the plot layout and remove top and right spines

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.grid(False)

    return fig

































