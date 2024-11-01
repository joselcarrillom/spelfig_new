<a id="readme-top"></a>

<h3 align="center">SpELFiG: Spectral Emission Line Fitting for Gas Spectra</h3>
    

  <p align="center">
    A prototype pipeline to explore emission line spectra of AGN outflows.
 <br />
    <a href="https://github.com/aksitadeo/spelfig_aksita"><strong>Explore the docs »</strong></a>
    <br />
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About The Project

The SpELFiG pipeline is designed to automatically fit emission line profiles of spectra from Active Galactic Nuclei (AGNs). Built for flexibility, it accepts a wide range of input data formats, resolutions, and file sizes, allowing users to seamlessly process data from various sources without manual adjustments. The scientific objective of this pipeline is to assist in the detection and characterisation of AGN outflows by analysing multi-component emission lines. Through investigating the spatially resolved kinematics of these emission lines, we can better understand how AGN-driven outflows might regulate star formation and influence galaxy growth in host galaxies. 

### Features
_Profile Options:_ Choose from Gaussian, Asymmetric Gaussian, Lorentzian, or Voigt profiles – adding as many components as you like – providing versatility to handle different emission line shapes observed in AGN spectra.

_Iterative Fitting Approach:_ Our pipeline employs an innovative iterative fitting method, enhancing accuracy with each cycle. The approach uses Monte-Carlo Markov-Chain (MCMC), systematically refining model parameters to yield high-quality fits, even for complex or blended lines.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

The code is currently running on Python 3.11.5. To check your version, copy and paste the following into your terminal:
```
python --version
python -V
```
Or into a notebook cell:
```
! python --version
```


For this pipeline, the main library package used is emcee. The recommended way to install the stable version of emcee is using pip:
  ```
python -m pip install -U pip
pip install -U setuptools setuptools_scm pep517
pip install -U emcee
  ```
Or conda:
```
conda update conda
conda install -c conda-forge emcee
```
Other libraries include numpy, scipy, matplotlib, pandas and astropy. 

### Installation

1. Fork and clone the repository onto your local machine
   ```sh
   git clone https://github.com/aksitadeo/spelfig_aksita.git
   ```
2. Navigate to the directory
   ```sh
   cd spelfig_aksita
   ```
3. Install relevant packages
   ```sh
   pip install emcee # for example
   ```
4. Change git remote url to avoid accidental pushes to base project
   ```sh
   git remote set-url origin aksitadeo/spelfig_aksita
   git remote -v # confirm the changes
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

For basic usage, first navigate to your directory containing the spectra you would like to run in .fits format, then import the SpELFiG pipeline. Open up the spl_config.py script, and choose how you would like to set up the configuration for a fit:
```
# Specify input directory in the case it is different from the default.

inputdir = your_dir
resdir = None

# Define emission lines, profiles and components you would like to fit, for example:
EMISSION_LINES = {
	'He-II,1': {'wavelength':[3202.15], 'components':['Gaussian',['Gaussian']},
	'He-II,2': {'wavelength':[4685.74], 'components':['Gaussian']},
	'Ar-IV,1':  {'wavelength':[4711.30], 'components':['Gaussian', 'Lorentzian']},
	'Ar-IV,2':  {'wavelength':[4740.10], 'components':['Gaussian', 'Lorentzian']},
	'H-β':     {'wavelength':[4861.32], 'components':['Gaussian']},
	'S-II,1':  {'wavelength': [6716.31], 'components':['Gaussian', 'Voigt', 'Gaussian']},
	'S-II,2':   {'wavelength':[6730.68], 'components':['Gaussian', 'Voigt']},
	'Ar-III': {'wavelength':[7135.67], 'components':['Voigt']},
}

# OPTIONAL:
# Define any MCMC initial conditions, for example:
INITIAL_CONDITIONS = {
	"powerlaw function": {"gamma": -0.01},
	'goodness of fit': {'diff bic cut': 10., 'chi cut': 2.0},
	'fitting': {'niter': 5000}
}
```

From there, you can simply execute the pipeline using the command:

```
python fittersexec.py
```
The pipeline outputs fitted spectra with corresponding residuals, $\chi^2$, and BIC metrics, along with a summary of fit parameters for each emission line.

### Example Fits

Here are some of the fits you can produce using our pipeline.

Gaussian Fit             |  Lorentzian Fit
:-------------------------:|:-------------------------:
![8318-2gauss](https://github.com/user-attachments/assets/c5cff6e4-c4a8-43e5-bf30-f74ae72b6923)  | ![8318-loren](https://github.com/user-attachments/assets/fbf8776b-a6ad-48c7-837d-05b9ecf47ff7)

_For a more comprehensive look on a prototype notebook, please refer to the [Example Spectra Notebook](https://github.com/aksitadeo/spelfig_aksita/blob/main/fitting_notebooks/spectrafitting_aksita.ipynb)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ROADMAP -->
## Roadmap
- [X] Writable configuration script
- [ ] Upgrade continuum fitting
- [ ] Automatic Emission Line Analysis
    - [ ] Kurtosis, Asymmetry, & FWHM calculation

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create! Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also open an issue with the tag "enhancement".

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Aksita Deo - [@AksitaDeo](https://twitter.com/AksitaDeo) - aksita.deo@students.mq.edu.au

Project Link: [https://github.com/aksitadeo/spelfig_aksita](https://github.com/aksitadeo/spelfig_aksita)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### This work was done in close collaboration with Jose-Luis Carillo Martinez


