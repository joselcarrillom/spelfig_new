import spl_config as spc
import spl_execsetupv3 as spe
import os
import numpy as np

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

            # snr & redshift
            redshift_range = np.arange(0.00, 2.49, 0.01)
            most_lines_dict = None
            max_num_lines = 0
            best_redshift = None
            realsnr = None

            # try emission line fitting first
            for red in redshift_range:
                try:
                    # try to match emission lines
                    snr, identified_lines = spe.extract_redshift_snr(data, red, line=None)
                    if identified_lines:
                        num_lines = len(identified_lines)
                        # capture list that has the most lines identified
                        if num_lines >= max_num_lines:
                            most_lines_dict = identified_lines  # Store emission lines
                            max_num_lines = num_lines
                            best_redshift = red                 # Store the redshift
                            realsnr = snr                       # Store SNR

                        # # if low confidence, try template matching
                        if num_lines < 5:
                            best_redshift = spe.redshift_calc(data, templates,filename, Plot=True)

                except Exception as e:
                    continue

            # add to dictionary
            spec_dict[filename] = {'data': data, 'redshift': best_redshift, 'snr': realsnr, 'emission lines' : most_lines_dict}

            continue

        else:
            continue

    return spec_dict

# def table_dictionary(spectrum, redshift, snr, objectid):
#     for
#     dict = {
#          str(objectid) : {'spectrum' : spectrum, 'redshift' : redshift, 'snr' : snr}
#     }
#     return dict
#
# if spc.table is not None:
#     spectrum = spc.table['spectrum']
#     redshift = spc.table['redshift']
#     snr = spc.table['snr']
#     objectid = spc.table['objectid']
#
#
# else:
