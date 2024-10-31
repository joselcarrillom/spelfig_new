import pandas as pd
import numpy as np
import emcee
import spl_models as spm  # Importing the analytical models

def spectral_model_emcee(theta, x, models, continuum):
    # We will go throughout theta in terms of the single distributions for the emissions
    # Initializing the flux:
    flux = np.zeros_like(x)
    param_start = 0
    continuum = list(continuum)

    # Adding the emission components:
    for model in models:
        if model == 'Gaussian':
            param_end = param_start + 3
            try:
                flux += spm.gauss(x, *theta[param_start:param_end])
            except:
                print('Error while calculating model, this complete set of parameters:')
                print('theta', theta)
                print('This subset at error', theta[param_start:param_end])
            param_start = param_end
        elif model == 'Asymmetric Gaussian':
            param_end = param_start + 4
            try:
                flux += spm.asym_gauss(x, *theta[param_start:param_end])
            except:
                print('Error while calculating model, this complete set of parameters:')
                print('theta', theta)
                print('This subset at error', theta[param_start:param_end])
            param_start = param_end
        elif model == 'Lorentzian':
            param_end = param_start + 3
            try:

                flux += spm.lorentzian(x, *theta[param_start:param_end])
            except:
                print('Error while calculating model, this complete set of parameters:')
                print('theta', theta)
                print('This subset at error', theta[param_start:param_end])
            param_start = param_end
        elif model == 'Voigt':
            param_end = param_start + 4
            try:
                flux += spm.voigt(x, *theta[param_start:param_end])
            except:
                print('Error while calculating model, this complete set of parameters:')
                print('theta', theta)
                print('This subset at error', theta[param_start:param_end])
            param_start = param_end

    # Adding the continuum:
    #try:

    #flux += spm.continuum_function(x, *continuum)
    #except:
    #print('Error while calculating continuum, this set of parameters:')
    #print(continuum)

    return flux


# HERE THE EMCEE FUNCTIONS TO EXECUTE THE CHAINS: ------------------------------
# Priors function:
def logpriors(theta, min_values, max_values, components):

    # First evaluate each parameter:
    for i, param in enumerate(theta):
        min_val = min_values[i]
        max_val = max_values[i]
        if not (min_val <= param <= max_val):
            return -np.inf

    # Evaluate the sum of amplitudes for lines with multiple components:
    param_start = 0
    for ncomp in components:
        params_l = ncomp * 3
        total_amplitude = 0
        for j in range(ncomp):
            maxflux = max_values[param_start + j * 3 + 1]
            total_amplitude += theta[param_start + j * 3 + 1]
            if total_amplitude > maxflux:
                return -np.inf
        param_start += params_l

    # Evaluate the ratio between amplitudes:
    # HERE future code to constrain ratio of line doublets.
    return 0.0

# The log_probablity function:
def log_prob(theta, x, y, dy, models, continuum):
    f = spectral_model_emcee(theta, x, models, continuum)
    return -0.5 * np.sum(((f - y) / dy) ** 2)

# The log_posterior function:
def log_posterior(theta, x, y, dy, models, continuum, min_values, max_values, components):
    lp = logpriors(theta, min_values, max_values, components)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_prob(theta, x, y, dy, models, continuum)


# The goodness of fit calculator:
def goodness_of_fit(theta, x, y, dy, models, continuum):
    p = len(theta)
    N = len(x)
    f = spectral_model_emcee(theta, x, models, continuum)
    chi2 = np.sum((y - f) ** 2 / dy ** 2)
    chi2dof = chi2 / (N - p)
    BIC = np.log(N) * p - np.log(np.exp(-0.5 * chi2)) * (N - p)
    return {'chi squared': chi2, 'reduced chi squared': chi2dof, 'BIC': BIC}

def emcee_sampler(theta0, nwalkers, niter, args):
    ndim = len(theta0[0])
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=args)
        # moves=[
        # (emcee.moves.DEMove(), 0.8),
        # (emcee.moves.DESnookerMove(), 0.2),
        # ])

    # Running burn-in
    p0, prob, state = sampler.run_mcmc(theta0, 100)
    sampler.reset()

    pos, prob, state = sampler.run_mcmc(p0, niter, progress=True)

    return sampler, pos, prob, state


def results_df(dfparams_init, theta_max, theta_errors, continuum):
    '''
    Final dataframe of results, in the format of the input dataframe (all the parameters
    condensed in one column, and without NaNs)
    '''

    # Create a dictionary to store the results
    lines_groups = dfparams_init['Line Name'].unique()
    print('This list of lines groups', lines_groups)

    results_dict = {
        'Line Name': [],
        'Model': [],
        'Component': [],
        'Parameters': [],
        'Parameter Errors': [],
    }

    # Populate the dictionary with the results
    param_start = 0

    for line in lines_groups:
        linedf = dfparams_init[dfparams_init['Line Name'] == line]
        linedf = linedf.reset_index(drop=True)
        Ncomp = len(linedf)
        if not line == 'Continuum':
            for j in range(Ncomp):
                results_dict['Line Name'].append(line)
                results_dict['Component'].append(j+1)
                model_j = linedf['Model'][j]
                results_dict['Model'].append(model_j)
                # Extracting parameters for the jth component:
                if (model_j == 'Gaussian') | (model_j == 'Lorentzian'):
                    results_dict['Parameters'].append(theta_max[param_start:param_start + 3])
                    results_dict['Parameter Errors'].append(theta_errors[param_start:param_start + 3])
                    param_start += 3
                elif (model_j == 'Voigt') | (model_j == 'Asymmetric Gaussian'):
                    results_dict['Parameters'].append(theta_max[param_start:param_start + 4])
                    results_dict['Parameter Errors'].append(theta_errors[param_start:param_start + 4])
                    param_start += 4


    results_dict['Line Name'].append('Continuum')
    results_dict['Model'].append('Continuum')
    results_dict['Component'].append(0)
    results_dict['Parameters'].append(continuum)
    results_dict['Parameter Errors'].append(np.zeros_like(continuum))

    return pd.DataFrame(results_dict)

def run_mcmc_chains(dfparams, x, y, dy, save_complete_chains=False, savefile=None, niter=1000):
    # First extract the parameters from the dataframe:
    lines = dfparams['Line Name'].unique()
    components = dfparams.groupby('Line Name')['Component'].max()
    models = dfparams['Model'].values
    continuum = dfparams['Parameters'][dfparams['Line Name'] == 'Continuum']
    p0 = np.concatenate(dfparams['Parameters'].values)
    # print('This p0', p0)

    # Min and max values for the parameters:
    min_values = np.concatenate(dfparams['Min Limits'].values)
    max_values = np.concatenate(dfparams['Max Limits'].values)

    #Preparing the args for the log_posterior function:
    args = [x, y, dy, models, continuum, min_values, max_values, components]

    # Now we are going to initialize the walkers:
    ndim = len(p0)
    nwalkers = 2*ndim

    # Try two wasy to initialize the walkers
    # theta0 = np.random.uniform(min_values, max_values, size=(nwalkers, ndim))
    theta0 = np.array([p0 + 1e-3*np.random.randn(ndim) for _ in range(nwalkers)])
    # print('This theta0 dimensions', np.shape(theta0))

    # Running the chains:
    sampler, pos, prob, state = emcee_sampler(theta0, nwalkers, niter, args)
    samples = sampler.flatchain

    # Best model:
    theta_max = samples[np.argmax(sampler.flatlnprobability)]

    # Errors in the model:
    theta_errors = np.std(samples, axis=0)

    # If specified, saving the ENTIRE samples in a file:
    if save_complete_chains is True:
        CHAINS = pd.DataFrame(samples)
        CHAINS.to_csv('complete_chains.csv', index=False)

    # theta_max = np.array(theta_max)
    continuum = theta_max[ndim - 3:]
    # theta_errors = np.array(theta_errors)

    model_parameters_df = results_df(dfparams, theta_max, theta_errors, continuum)
    print("Successful! Convergence found.")

    return model_parameters_df, theta_max, theta_errors

class mcmc_fit(object):
    '''
    This class is meant to produce the object containing all relevant outcomes from the mcmc fit.
    '''
    def __init__(self, dfparams, x, y, dy, save_complete_chains=False, savefile=None, niter=1000):
        self.dfparams = dfparams
        self.continuum = dfparams['Parameters'].loc[dfparams['Line Name'] == 'Continuum'].values[0]
        self.x = x
        self.y = y
        self.dy = dy
        self.save_complete_chains = save_complete_chains
        self.savefile = savefile
        self.niter = niter
        self.model_parameters_df, self.fit_parameters, self.fit_errors = run_mcmc_chains(
            self.dfparams, self.x, self.y, self.dy, self.save_complete_chains, self.savefile,
            self.niter)

        self.models = self.model_parameters_df['Model'].values
        self.goodness = goodness_of_fit(self.fit_parameters, self.x, self.y, self.dy, self.models,
                                        self.continuum)


