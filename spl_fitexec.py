import os
import numpy as np
import sps_fitsetup as sps
import spl_config as spc
import spl_fitters as spf

mcfit = spf.mcmc_fit

spectradict = None # IMPORT HERE THE SUPER DICTIONARY FROM AKSITA'S MODULES:

# Extracting important variables: ------------------------------------------------------------------
inputdir0 = spc.inputdir
resdir0 = spc.resdir


# Getting the input directory with all the spectra to be fitted:
inputdir = inputdir0 if inputdir0 is not None else os.getcwd()+'/spl_input/'
resdir = resdir0 if resdir0 is not None else os.getcwd()+'/spl_output/'

# Create the output directory if it doesn't exist:
if not os.path.exists(resdir):
    os.mkdir(resdir)

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


    ## Run a first plot with the baseline model:
    print('Fitting model with just one model...')
    fit0 = mcfit(df0, x, y, dy, niter=niter)
    
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

def multiple_spectra_fitting(superdict, emlines, initial_conditions, specrange):
    # Extract the ID's list:
    IDS = superdict.keys()

    # Perform the fitting:
    for ID in IDS:
        spectrum = superdict[ID]['SPECTRUM']
        redshift = superdict[ID]['REDSHIFT']
        print(f'>>>>>   Executing cycle of fits for object: {ID}   <<<<<')
        fit = single_fit(spectrum, emlines, redshift, initial_conditions, specrange)
        # Store the fit in the dictionary:
        superdict[ID]['BEST FIT'] = fit
    return superdict

def file_saver(superdict, savefitfiles=True, saveplotfiles=True, plotranges=None):
    # Extract the ID's list:
    IDS = superdict.keys()

    for ID in IDS:
        # Create filename:
        fitfilename = resdir + ID + '_parfile.csv'
        plotfilename = resdir + ID + '_fitplot.png'
        # Extract the fit:
        fit = superdict[ID]['BEST FIT']

        # Save the fit:
        if savefitfiles is True:
            sps.spl_savefile(fit, fitfilename)
        # Save the plot:
        if saveplotfiles is True:
            # Extract the data:
            x = superdict[ID]['SPECTRUM'][:, 0]
            y = superdict[ID]['SPECTRUM'][:, 1]
            dy = superdict[ID]['SPECTRUM'][:, 2]
            dfparams = fit.model_parameters_df
            if plotranges:
                x_zoom = plotranges[0]
                y_zoom = plotranges[1]
            else:
                x_zoom = None
                y_zoom = None
            spl_fig = sps.spl_plot(x, y, dy, dfparams, x_zoom=x_zoom, y_zoom=y_zoom,
                            goodness_marks=fit.goodness)
            spl_fig.savefig(plotfilename, format = 'png', dpi=400, bbox_inches='tight')


###  EXAMPLE DICTIONARY TO TEST THE CODE: **************************************

### IN CASE NO DICTIONARY IS PROVIDED, THE CODE WILL EXECUTE THIS PART WHICH IS THE EXAMPLE TEST
if spectradict is None:
    from astropy.io import fits
    import re

    specfiles = ['manga-7815-6104_Reff_spec.fits',
    'manga-7978-12705_Reff_spec.fits',
    'manga-8089-12705_Reff_spec.fits']

    # Extract the correct ID from the file names:
    # Regular expression pattern for MaNGA data:
    pattern = r'(manga-\d+-\d+)'

    # Extract the pattern from each filename
    ids = [re.search(pattern, filename).group(1) for filename in specfiles]

    dir = os.getcwd()
    specrange = [4800, 6800]
    plotrange = ([6500, 6800], None)

    # Build the dictionary:
    superdict_test = {}

    for i, specfile in enumerate(specfiles):
        id = ids[i]
        hdu = fits.open(dir+'/spl_test_example/'+specfile)
        spectra = hdu[1].data
        spec = np.array([spectra['wavelength'], spectra['emlines'], spectra['flux_error']]).T
        redshift = hdu[1].header['Z']
        superdict_test[id] = {'SPECTRUM': spec, 'REDSHIFT': redshift}

    print('>>>>   Executing the fits on the example data set   <<<<<')
    superdict_test = multiple_spectra_fitting(superdict_test, emlines, initial_conditions, specrange)
    print('>>>>  Producing the plots and parameter files for the example data set   <<<<<')
    file_saver(superdict_test, savefitfiles=True, saveplotfiles=True, plotranges=plotrange)
    print('Test example done!')


