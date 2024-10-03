import spl_config as spc
import spl_exsetupv2 as spe
import os

def spectrum_dictionary(directory,templates):

    '''
    Creates a dictionary containing key info such as filename, spectra, redshift and SNR.

    Args:
        directory (str): The directory containing the FITS files to analyse.
        templates (list): A list of `specutils.Spectrum1D` templates for redshift calculations.
    Returns:
        dictionary: A dictionary containing the spectra, redshift and SNR.
    Raises:
        UnboundLocalError: Occurs when no peak has been identified.
    '''

    spec_dict = {}

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".fits"):

            # id
            print("Now analysing ", os.path.join(directory, filename))
            full_dir = os.path.join(directory, filename)

            # spectra
            data = spe.extract_astronomical_data(full_dir)

            # redshift
            redshift = spe.redshift_calc(data, templates, Plot=True)
            # snr
            try:
                snr = spe.extract_snr(data, redshift, line='O-III,0')
            except:
                print('UnboundLocalError')
                snr = 0

            # add to dictionary
            spec_dict[filename] = {'data': data, 'redshift': redshift, 'snr': snr}

            continue

        else:
            continue

    return spec_dict

def table_dictionary(spectrum, redshift, snr, objectid):
    for
    dict = {
         str(objectid) : {'spectrum' : spectrum, 'redshift' : redshift, 'snr' : snr}
    }
    return dict

if spc.table is not None:
    spectrum = spc.table['spectrum']
    redshift = spc.table['redshift']
    snr = spc.table['snr']
    objectid = spc.table['objectid']



else:
