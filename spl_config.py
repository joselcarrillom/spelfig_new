## USER SCRIPT ###
'''
This is meant to be the single spectrum setup script for the user for the SpELFiG
code. It aims to gather the necessary information from the user to set up the configuration for
the fits. its elements are as follows:

The EMISSION_LINES dictionary list. It contains all the emission lines potentially to be present
at the spectra. Not necessary in all of them.
'''

# HERE: SPECIFY THE INPUT DIRECTORY in the case is other different from the default.

inputdir = None


EMISSION_LINES = {
	'He-II,1': {'wavelength':[3202.15], 'components':['Gaussian']},
	'He-II,2': {'wavelength':[4685.74], 'components':['Gaussian']},
	'Ne-V,1':  {'wavelength':[3345.81], 'components':['Gaussian']},
	'Ne-V,2':  {'wavelength':[3425.81], 'components':['Gaussian']},
	'O-II,1':  {'wavelength':[3726.03], 'components':['Gaussian']},
	'O-II,2':  {'wavelength':[3728.73], 'components':['Gaussian']},
	'Ne-III,1': {'wavelength':[3868.69], 'components':['Gaussian']},
	'Ne-III,2': {'wavelength':[3967.40], 'components':['Gaussian']},
	'H-ζ':    {'wavelength':[3889.05], 'components':['Gaussian']},
	'H-ε':    {'wavelength':[3970.07], 'components':['Gaussian']},
	'H-δ':    {'wavelength':[4101.73], 'components':['Gaussian']},
	'H-γ':    {'wavelength':[4340.46], 'components':['Gaussian']},
	'O-III,0':  {'wavelength':[4363.15], 'components':['Gaussian']},
	'O-III,1':  {'wavelength':[4958.83], 'components':['Gaussian']},
	'O-III,2':  {'wavelength':[5006.77], 'components':['Gaussian']},
	'Ar-IV,1':  {'wavelength':[4711.30], 'components':['Gaussian', 'Lorentzian']},
	'Ar-IV,2':  {'wavelength':[4740.10], 'components':['Gaussian', 'Lorentzian']},
	'H-β':     {'wavelength':[4861.32], 'components':['Gaussian']},
	'N-I,1':    {'wavelength':[5197.90], 'components':['Gaussian']},
	'N-I,2':    {'wavelength':[5200.39], 'components':['Gaussian']},
	'He-I':   {'wavelength':[5875.60], 'components':['Gaussian']},
	'O-I,1':   {'wavelength':[6300.20], 'components':['Gaussian']},
	'O-I,2':   {'wavelength':[6363.67], 'components':['Gaussian']},
	'N-II,1':   {'wavelength':[6547.96], 'components':['Gaussian']},
	'N-II,2':   {'wavelength':[6583.34], 'components':['Gaussian']},
	'H-α':     {'wavelength':[6562.80], 'components':['Gaussian']},
	'S-II,1':  {'wavelength': [6716.31], 'components':['Gaussian', 'Voigt', 'Gaussian']},
	'S-II,2':   {'wavelength':[6730.68], 'components':['Gaussian', 'Voigt']},
	'Ar-III': {'wavelength':[7135.67], 'components':['Gaussian']},
	'O VI': {'wavelength': [1033.82], 'components':['Gaussian']},
	'Ly-alpha': {'wavelength': [1215.24], 'components':['Gaussian']},
	'N V': {'wavelength': [1305.53], 'components':['Gaussian']},
	'Si IV + O IV': {'wavelength': [1549.48], 'components':['Gaussian']},
	'C IV': {'wavelength': [1640.40], 'components':['Gaussian']},
	'C III': {'wavelength': [2326.00], 'components':['Gaussian']},
	'Mg II': {'wavelength': [3346.79], 'components':['Gaussian']},
}



INITIAL_CONDITIONS = {
	"powerlaw function": {"gamma": -0.01},
	'goodness of fit': {'diff bic cut': 10., 'chi cut': 2.0},
	'fitting': {'niter': 5000}
}









