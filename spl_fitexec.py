import os
import numpy as np
import sps_fitsetup as sps
import spl_config as spc
import spl_fitters as spf

mcfit = spf.mcmc_fit

# Extracting important variables: ------------------------------------------------------------------
inputdir0 = spc.inputdir

# Getting the input directory with all the spectra to be fitted:
inputdir = inputdir0 if inputdir0 is not None else os.getcwd()+'/spl_input/'


# Initializing the process by extracting the spectra and redshift make a redshift correction if
# necessary: ---------------------------------------------------------------------------------------

# Line list:
# This is the line list that will be used in the fitting process:
emlines = spc.EMISSION_LINES
initial_conditions = spc.INITIAL_CONDITIONS

# Component extractor: The following is a method that extracts separate lists of the components.
def extract_component_lists(emission_lines):
    """
    Extract lists of emission lines categorized by their extra components.

    Args:
        emission_lines (dict): Dictionary containing emission lines and their components.

    Returns:
        List[Dict[str, List[str]]]: A list of dictionaries where each dictionary contains emission lines
                                    and their components at a specific level.
    """
    # Initialize dictionary to store lines with their extra components
    base_components = {}
    add_components = {}

    # Populate the add_components dictionary with the lines and their components
    for line, details in emission_lines.items():
        wavelength = details['wavelength']
        components = details['components']
        # Store baseline model:
        base_components[line] = {'wavelength': wavelength, 'components': components[0]}
        # component
        if len(components) > 1:
            # Store only lines with extra components
            add_components[line] = components[1:]  # Exclude 'Gaussian'

    # Initialize dictionaries to store lists of lines for each component level
    component_lists = {}

    # Process each emission line with extra components
    for line, components in add_components.items():
        for idx, comp in enumerate(components):
            if idx not in component_lists:
                component_lists[idx] = {}
            if comp not in component_lists[idx]:
                component_lists[idx][comp] = []
            component_lists[idx][comp].append(line)

    # Create a list of dictionaries for each level of components
    output_lists = []
    for level, comps in component_lists.items():
        level_dict = {}
        for comp, lines in comps.items():
            for line in lines:
                level_dict[line] = [comp]
        output_lists.append(level_dict)

    return base_components, output_lists

# ITERATIVE FUNCTIONS: ============================================================================
# We are going to define the iterative functions to perform the fitting process:
def single_fit(spectrum, emission_lines, redshift, initial_conditions, specrange):

    '''

    :param spectrum:
    :param emission_lines:
    :return: It returns the resulting best fit from a single spectrum. It's purpose is to be used


    Notes: This function performs the fit algorithm for a single spectrum
    '''

    # Store the parts of the spectrum:
    spec0 = spectrum

    # Correct by redshift:
    spec0[:, 0] = spec0[:, 0] / (1 + redshift)

    # Store the parts of the spectrum:
    x = spec0[:, 0]
    y = spec0[:, 1]
    dy = spec0[:, 2]

    # First of all we are going to extract the individual list of components:
    base_components, extra_component_lists = extract_component_lists(emission_lines)

    ## Spectral index for the continuum:

    gamma0 = initial_conditions['powerlaw function']['gamma']
    bic_cut = initial_conditions['goodness of fit']['diff bic cut']
    chi_cut = initial_conditions['goodness of fit']['chi cut']
    niter = initial_conditions['fitting']['niter']

    ## Observed data stored in variables:

    try:
        df0 = sps.init_setup(spec0, base_components, specrange, gamma0)
    except RuntimeWarning:
        print("Error occurred during setup. Skipping...")
        df0 = np.empty(0)

    print('This initial dataframe')
    print(df0)
    ## Run a first plot with the baseline model:
    print('Fitting model with just one model...')
    fit0 = mcfit(df0, x, y, dy, niter=niter)
    print('After first fit')
    
    ## Extract the goodness of fit:
    goodness0 = fit0.goodness
    chi2 = goodness0['reduced chi squared']
    bic = goodness0['BIC']
    diff_bic = np.inf - bic
    df0 = fit0.model_parameters_df


    ## ADDING MORE COMPONENTS:
    i = 0
    fit_final = fit0
    while (diff_bic > bic_cut or chi2 > chi_cut) and i < len(extra_component_lists):
        print('Testing increasing number of components...')
        newdf = sps.update_components(df0, extra_component_lists[i])
        newfit = mcfit(newdf, x, y, dy, niter=niter)
        newbic = newfit.goodness['BIC']
        newchi2 = newfit.goodness['reduced chi squared']
        diff_bic = newbic - bic
        print('Goodness of fit: bic: ', newbic, ' chi2: ', newchi2, ' diff bic: ', diff_bic)
        if diff_bic < bic_cut or newchi2 < chi_cut:
            fit_final = newfit
            print('Converged to best fit!')
            return fit_final
        else:
            print('Increasing number of components...')
            df0 = newfit.model_parameters_df
            bic = newbic
        i += 1
    return fit_final

'''
# TEST TO VERIFY THE SINGLE SPECTRUM FIT:
import numpy as np
from astropy.io import fits

specfile = 'manga-7815-6104_Reff_spec.fits'
dir = os.getcwd()
hdu = fits.open(dir+'/spl_input/'+specfile)
spectra = hdu[1].data
spec = np.array([spectra['wavelength'], spectra['emlines'], spectra['flux_error']]).T

redshift = 0.07971366494894028
specrange = [4800, 6800]


test1 = single_fit(spec, emlines, redshift, initial_conditions, specrange)

print('Test for single fit: DONE')
'''

def multiple_spectra_fitting(superdict, emlines, initial_conditions, specrange):
    # Extract the ID's list:
    IDS = superdict.keys()

    # Perform the fitting:
    for ID in IDS:
        spectrum = superdict[ID]['DATA']
        redshift = superdict[ID]['REDSHIFT']
        fit = single_fit(spectrum, emlines, redshift, initial_conditions, specrange)
        # Store the fit in the dictionary:
        superdict[ID]['BEST FIT'] = fit.model_parameters_df
        superdict[ID]['BEST FIT CHI2'] = fit.goodness['reduced chi squared']
        superdict[ID]['BEST FIT BIC'] = fit.goodness['BIC']


def save_and_plot(superdict, saveflags):
    pass

'''
def lines_df_toprint(resdfparams):

    # Create a dictionary to store the results
    results_dict = {
        'Line Name': [],
        'Model': [],
        'Component': [],
        'Centroid': [],
        'Amplitude': [],
        'Sigma': [],
        'Gamma': [],
        'err_Centroid': [],
        'err_Amplitude': [],
        'err_Sigma': [],
        'err_Gamma': [],
    }

    # Populate the dictionary with the results
    param_start = 0
    for i, line in enumerate(lines):
        Ncomp = components[i]
        for j in range(Ncomp):
            results_dict['Line Name'].append(line)
            results_dict['Component'].append(j)
            models[j] = models[param_start + j]
            results_dict['Model'].append(models[j])
            # Extracting parameters:
            # Centroid:
            centroid = theta_max[param_start + j * 3]
            err_centroid = theta_errors[param_start + j * 3]
            # Amplitude:
            amplitude = theta_max[param_start + j * 3 + 1]
            err_amplitude = theta_errors[param_start + j * 3 + 1]

            # Sigma and Gamma:
            if models[j] == 'Gaussian':
                pn = 3
                sigma = theta_max[param_start + j * 3 + 2]
                err_sigma = theta_errors[param_start + j * 3 + 2]
                gamma = np.nan
                err_gamma = np.nan
            elif models[j] == 'Lorentzian':
                pn = 3
                gamma = theta_max[param_start + j * 3 + 2]
                err_gamma = theta_errors[param_start + j * 3 + 2]
                sigma = np.nan
                err_sigma = np.nan
            elif models[j] == 'Voigt':
                pn = 4
                sigma = theta_max[param_start + j * 3 + 2]
                gamma = theta_max[param_start + j * 3 + 3]
                err_sigma = theta_errors[param_start + j * 3 + 2]
                err_gamma = theta_errors[param_start + j * 3 + 3]

            sigma_kms = spm.vel_correct(sigma, centroid)
            gamma_kms = spm.vel_correct(gamma, centroid)
            err_sigma_kms = spm.vel_correct(err_sigma, centroid)
            err_gamma_kms = spm.vel_correct(err_gamma, centroid)

            # Appending results:
            results_dict['Centroid'].append(centroid)
            results_dict['Amplitude'].append(amplitude)
            results_dict['Sigma'].append(sigma)
            results_dict['Sigma (km/s)'].append(sigma_kms)
            results_dict['Gamma'].append(gamma)
            results_dict['Gamma (km/s)'].append(gamma_kms)
            results_dict['err_Centroid'].append(err_centroid)
            results_dict['err_Amplitude'].append(err_amplitude)
            results_dict['err_Sigma'].append(err_sigma)
            results_dict['err_Sigma (km/s)'].append(err_sigma_kms)
            results_dict['err_Gamma'].append(err_gamma)
            results_dict['err_Gamma (km/s)'].append(err_gamma_kms)

        params_per_line = pn * Ncomp
        param_start += params_per_line

    return pd.DataFrame(results_dict)
 '''



