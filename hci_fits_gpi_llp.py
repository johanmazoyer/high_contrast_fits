import glob
from datetime import date
import numpy as np
from astropy.io import fits

#######################################################################################################
##### J spec KLIP ADI
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/Work/python/python_data/disk_mcmc/resultats_mcmc_2001010/Aurora/160323_J_Spec/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/Work/IDL/IDL_data/hr4796/pyklip_reduc/160323_J_Johan_reduc/pyklipADIonly/hr4796-imageandnoise-pyklipADIonly-160323-J-KL3.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map = fits.getdata(reduc_data_data_snr)[1]
# image_snr = image_reduc / image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'KLIP'  # Reduction algorithm
# REDSTRAT = 'ADI'  #Strategy to build the PSF library
# EXCLANG = 3  # Exclusion angle for ADI-type strategy (deg)
# TKL = 3  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 59  # Number of images in the PSF library
# NOISEMET = 'spatial'  # Method used for the detection limit
###################################################################################################

# #######################################################################################################
# ##### H spec KLIP ADI
# #######################################################################################################
raw_data_dir = '/Users/jmazoyer/Dropbox/Work/python/python_data/disk_mcmc/resultats_mcmc_2001010/Aurora/160318_H_Spec/'
reduc_data_data_snr = '/Users/jmazoyer/Dropbox/Work/IDL/IDL_data/hr4796/pyklip_reduc/160318_H_Johan_reduc/pyklipADIonly/hr4796-imageandnoise-pyklipADIonly-160318-H-KL3.fits'

image_reduc = fits.getdata(reduc_data_data_snr)[0]
image_sensitivity_map= fits.getdata(reduc_data_data_snr)[1]
image_snr = image_reduc/image_sensitivity_map

name_investigator = 'Bruce Macintosh'
BIBREF  = '2020AJ....160...79A'
Polarization = 'I'
REDALGO = 'KLIP' # Reduction algorithm
REDSTRAT = 'ADI' #Strategy to build the PSF library
EXCLANG = 3 # Exclusion angle for ADI-type strategy (deg)
TKL = 3  #Number of eigenimages used for KLIP-type subtraction
NPSFLIBR = 37 # Number of images in the PSF library
NOISEMET = 'spatial'  # Method used for the detection limit
# ###################################################################################################

#######################################################################################################
##### K1 spec KLIP ADI 140325
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/140325_K1_DRP_Elodie/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/Work/IDL/IDL_data/hr4796/pyklip_reduc/140325_K1_Johan_reduc/pyklipADIonly/hr4796-imageandnoise-pyklipADIonly-140325-K1-KL3.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[1]
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Bruce Macintosh'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'KLIP' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 3 # Exclusion angle for ADI-type strategy (deg)
# TKL = 3  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 20 # Number of images in the PSF library
# NOISEMET = 'spatial'  # Method used for the detection limit
###################################################################################################

#######################################################################################################
##### K1 spec KLIP ADI 150402
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/150402_K1_DRP_Elodie/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/Work/IDL/IDL_data/hr4796/pyklip_reduc/150402_K1_Johan_reduc/pyklipADIonly/hr4796-imageandnoise-pyklipADIonly-150402-K1-KL3.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[1]
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'KLIP' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 3 # Exclusion angle for ADI-type strategy (deg)
# TKL = 3  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 47 # Number of images in the PSF library
# NOISEMET = 'spatial'  # Method used for the detection limit
###################################################################################################

#######################################################################################################
##### K1 spec KLIP ADI 150403
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/150403_K1_DRP_Elodie/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/Work/IDL/IDL_data/hr4796/pyklip_reduc/150403_K1_Johan_reduc/pyklipADIonly/hr4796-imageandnoise-pyklipADIonly-150403-K1-KL3.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[1]
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'KLIP' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 3 # Exclusion angle for ADI-type strategy (deg)
# TKL = 3  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 48 # Number of images in the PSF library
# NOISEMET = 'spatial'  # Method used for the detection limit
###################################################################################################

#######################################################################################################
##### K2 spec KLIP ADI 150403
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/150403_K2_Spec_BkgdSub/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/Work/IDL/IDL_data/hr4796/pyklip_reduc/150403_K2_Johan_reduc/pyklipADIonly/hr4796-imageandnoise-pyklipADIonly-150403-K2-KL3.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[1]
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'KLIP' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 3 # Exclusion angle for ADI-type strategy (deg)
# TKL = 3  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 23 # Number of images in the PSF library
# NOISEMET = 'spatial'  # Method used for the detection limit
###################################################################################################

# #######################################################################################################
# ##### H spec KLIP RDI
# #######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/Work/python/python_data/disk_mcmc/resultats_mcmc_2001010/Aurora/160318_H_Spec/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/Deliveries/HR_4796A/Spec/PSF Subtracted/RDI-KLIP/HR4796A_rdi-KLmodes-all-87-OWAKL7_data_mask_noise_snr.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[2]
# image_snr = fits.getdata(reduc_data_data_snr)[3]

# name_investigator = 'Bruce Macintosh'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'KLIP' # Reduction algorithm
# REDSTRAT = 'RDI' #Strategy to build the PSF library
# EXCLANG = -1 # Exclusion angle for ADI-type strategy (deg)
# TKL = 7  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 99 # Number of images in the PSF library
# NOISEMET = 'spatial'  # Method used for the detection limit
# ###################################################################################################

# #######################################################################################################
# ##### H spec KLIP NMF
# #######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/Work/python/python_data/disk_mcmc/resultats_mcmc_2001010/Aurora/160318_H_Spec/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/Work/IDL/IDL_data/hr4796/RDI_NMF_Bin/hr4796-imageandnoise-nmfBin-160318-H-KL150.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[1]
# image_sensitivity_map[np.where(image_sensitivity_map == 0)] = np.nan
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Bruce Macintosh'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'NMF' # Reduction algorithm
# REDSTRAT = 'RDI' #Strategy to build the PSF library
# EXCLANG = -1 # Exclusion angle for ADI-type strategy (deg)
# TKL = 15  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 150 # Number of images in the PSF library
# NOISEMET = 'spatial'  # Method used for the detection limit
# ###################################################################################################

#######################################################################################################
##### J spec mask and interpolate ADI
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/Work/python/python_data/disk_mcmc/resultats_mcmc_2001010/Aurora/160323_J_Spec/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/Deliveries/HR_4796A/Spec/PSF Subtracted/Mask-And-Interpolate/hr4796a_j_160323.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[2]
# image_sensitivity_map[np.where(image_sensitivity_map == 0)] = np.nan
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'mask_and_interpolate' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 0 # Exclusion angle for ADI-type strategy (deg)
# TKL = 0  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 59 # Number of images in the PSF library
# NOISEMET = 'temporal'  # Method used for the detection limit
###################################################################################################

# #######################################################################################################
# ##### H spec mask and interpolate ADI
# #######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/Work/python/python_data/disk_mcmc/resultats_mcmc_2001010/Aurora/160318_H_Spec/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/Deliveries/HR_4796A/Spec/PSF Subtracted/Mask-And-Interpolate/hr4796a_h_160318.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[2]
# image_sensitivity_map[np.where(image_sensitivity_map == 0)] = np.nan
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Bruce Macintosh'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'mask_and_interpolate' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 0 # Exclusion angle for ADI-type strategy (deg)
# TKL = 0  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 37 # Number of images in the PSF library
# NOISEMET = 'temporal'  # Method used for the detection limit
# ###################################################################################################

#######################################################################################################
##### K1 spec KLIP mask and interpolate 140325
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/140325_K1_DRP_Elodie/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/Deliveries/HR_4796A/Spec/PSF Subtracted/Mask-And-Interpolate/hr4796a_k1_140325.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[2]
# image_sensitivity_map[np.where(image_sensitivity_map == 0)] = np.nan
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Bruce Macintosh'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'mask_and_interpolate' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 0 # Exclusion angle for ADI-type strategy (deg)
# TKL = 0  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 20 # Number of images in the PSF library
# NOISEMET = 'temporal'  # Method used for the detection limit
###################################################################################################

#######################################################################################################
##### K1 spec KLIP mask and interpolate 150402
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/150402_K1_DRP_Elodie/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/Deliveries/HR_4796A/Spec/PSF Subtracted/Mask-And-Interpolate/hr4796a_k1_150402.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[2]
# image_sensitivity_map[np.where(image_sensitivity_map == 0)] = np.nan
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'mask_and_interpolate' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 0 # Exclusion angle for ADI-type strategy (deg)
# TKL = 0  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 47 # Number of images in the PSF library
# NOISEMET = 'temporal'  # Method used for the detection limit
###################################################################################################

#######################################################################################################
##### K1 spec KLIP mask and interpolate 150403
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/150403_K1_DRP_Elodie/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/Deliveries/HR_4796A/Spec/PSF Subtracted/Mask-And-Interpolate/hr4796a_k1_150403.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[2]
# image_sensitivity_map[np.where(image_sensitivity_map == 0)] = np.nan
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'mask_and_interpolate' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 0 # Exclusion angle for ADI-type strategy (deg)
# TKL = 0  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 48 # Number of images in the PSF library
# NOISEMET = 'temporal'  # Method used for the detection limit
###################################################################################################

#######################################################################################################
##### K2 spec KLIP mask and interpolate 150403
#######################################################################################################
# raw_data_dir = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/HR_4796/150403_K2_Spec_BkgdSub/'
# reduc_data_data_snr = '/Users/jmazoyer/Dropbox/GPIDATA-LLP/Deliveries/HR_4796A/Spec/PSF Subtracted/Mask-And-Interpolate/hr4796a_k2_150403.fits'

# image_reduc = fits.getdata(reduc_data_data_snr)[0]
# image_sensitivity_map= fits.getdata(reduc_data_data_snr)[2]
# image_sensitivity_map[np.where(image_sensitivity_map == 0)] = np.nan
# image_snr = image_reduc/image_sensitivity_map

# name_investigator = 'Christine Chen'
# BIBREF  = '2020AJ....160...79A'
# Polarization = 'I'
# REDALGO = 'mask_and_interpolate' # Reduction algorithm
# REDSTRAT = 'ADI' #Strategy to build the PSF library
# EXCLANG = 0 # Exclusion angle for ADI-type strategy (deg)
# TKL = 0  #Number of eigenimages used for KLIP-type subtraction
# NPSFLIBR = 23 # Number of images in the PSF library
# NOISEMET = 'temporal'  # Method used for the detection limit
###################################################################################################

filelist_raw = sorted(glob.glob(raw_data_dir + '*_distorcorr.fits'))
parangs = []
integration_time = []
mjd_start = []
ut_time = []
for index_angle in range(0, len(filelist_raw)):

    first_header_of_the_raw_fits = fits.getheader(filelist_raw[index_angle], 0)
    second_header_of_the_raw_fits = fits.getheader(filelist_raw[index_angle],
                                                   1)

    mjd_start.append(first_header_of_the_raw_fits['MJD-OBS'])

    parangs.append(second_header_of_the_raw_fits['AVPARANG'])
    integration_time.append(second_header_of_the_raw_fits['ITIME'])
    ut_time.append(second_header_of_the_raw_fits['EXPSTART'])
    image = fits.getdata(filelist_raw[index_angle])

Combined_Rotation_Angle = np.nanmax(parangs) - np.nanmin(parangs)
Number_of_Exposures = len(filelist_raw)
Exposure_Time = np.nansum(integration_time)
Observation_Start = mjd_start[0]
Observation_End = mjd_start[-1]

UT_Midpoint_Date_of_Observation = first_header_of_the_raw_fits['DATE-OBS']
UT_Midpoint_Time_of_Observation = ut_time[len(filelist_raw) // 2].split('.')[0]
datecompact = UT_Midpoint_Date_of_Observation.replace("-", "")

PROPOSID = first_header_of_the_raw_fits['GEMPRGID']
TELESCOP = first_header_of_the_raw_fits['TELESCOP']
INSTRUME = first_header_of_the_raw_fits['INSTRUME']
obs_mode = first_header_of_the_raw_fits['OBSMODE']
file_type = first_header_of_the_raw_fits['FILETYPE']
OBSSTGY = obs_mode + '_spec' if file_type == 'Spectral Cube' else obs_mode + '_pol'

PR_INV_F = name_investigator.split(' ')[0]
PR_INV_L = name_investigator.split(' ')[1]

PIXSCALE = 0.01414  # GPI pixscale Arcsecond / pix

FILTER = obs_mode.split('_')[0]

if FILTER == 'Y':
    Wavelength = 1.04
    Bandwidth = 0.19
elif FILTER == 'J':
    Wavelength = 1.23
    Bandwidth = 0.23
elif FILTER == 'H':
    Wavelength = 1.64
    Bandwidth = 0.30
elif FILTER == 'K1':
    Wavelength = 2.04
    Bandwidth = 0.29
elif FILTER == 'K2':
    Wavelength = 2.26
    Bandwidth = 0.27

TARGNAME = first_header_of_the_raw_fits['OBJECT']
RA_TARG = first_header_of_the_raw_fits['RA']
DEC_TARG = first_header_of_the_raw_fits['DEC']
EQUINOX = first_header_of_the_raw_fits['EQUINOX']

########################################################################################
# creation of primary header
########################################################################################
Primary_header = fits.Header()
Primary_header.append(('', ''), end=True)
Primary_header.append(('', 'FILE INFORMATION'), end=True)
Primary_header.append(('', ''), end=True)
Primary_header.append(('FILETYPE', 'HLSP', 'type of data found in data file'),
                      end=True)
Primary_header.append(('ORIGIN', 'LLP team', 'FITS file originator'), end=True)
DATE = date.today().strftime("%Y-%m-%d")
Primary_header.append(
    ('DATE', DATE, 'Date this file was written (yyyy-mm-dd)'), end=True)

Primary_header.append(('', ''), end=True)
Primary_header.append(('', 'FITS FILE DATA STRUCTURE'), end=True)
Primary_header.append(('', ''), end=True)
Primary_header.append(('NEXTEND', 4, 'Number of extensions'), end=True)
Primary_header.append(('EXT1NAME', 'DATA_INFORMATION', 'Extension 1 name'),
                      end=True)
Primary_header.append(('EXT1TYPE', 'BINTABLE', 'Extension 1 type'), end=True)

Primary_header.append(('EXT1NAME', 'DATA_INFORMATION', 'Extension 1 name'),
                      end=True)
Primary_header.append(('EXT1TYPE', 'BINTABLE', 'Extension 1 type'), end=True)

Primary_header.append(('EXT2NAME', 'REDUCED_DATA', 'Extension 2 name'),
                      end=True)
Primary_header.append(('EXT2TYPE', 'IMAGE', 'Extension 2 type'), end=True)

Primary_header.append(('EXT3NAME', 'SENSITIVITY_MAP', 'Extension 3 name'),
                      end=True)
Primary_header.append(('EXT3TYPE', 'IMAGE', 'Extension 3 type'), end=True)

Primary_header.append(('EXT4NAME', 'SNR_MAP', 'Extension 4 name'), end=True)
Primary_header.append(('EXT4TYPE', 'IMAGE', 'Extension 4 type'), end=True)

Primary_header.append(('', ''), end=True)
Primary_header.append(('', 'PROGRAM AND INSTRUMENT INFORMATION'), end=True)
Primary_header.append(('', ''), end=True)

Primary_header.append(('TELESCOP', TELESCOP, 'telescope used to acquire data'),
                      end=True)
Primary_header.append(
    ('INSTRUME', INSTRUME, 'instrument used to acquire data'), end=True)
Primary_header.append(('PROPOSID', PROPOSID, 'PEP proposal identifier'),
                      end=True)
Primary_header.append(
    ('PR_INV_L', PR_INV_L, 'last name of principal investigator '), end=True)
Primary_header.append(
    ('PR_INV_F', PR_INV_F, 'first name of principal investigator '), end=True)
Primary_header.append(('OBSSTGY', OBSSTGY, 'observation strategy'), end=True)
Primary_header.append(('FILTER', FILTER, 'filter used during observation'),
                      end=True)
Primary_header.append(('BIBREF', BIBREF, 'ADS bibliographic reference'),
                      end=True)

Primary_header.append(('', ''), end=True)
Primary_header.append(('', 'TARGET INFORMATION'), end=True)
Primary_header.append(('', ''), end=True)

Primary_header.append(('TARGNAME', TARGNAME, 'target name'), end=True)
Primary_header.append(('RA_TARG', RA_TARG, 'RA of target (deg) (J2000)'),
                      end=True)
Primary_header.append(('DEC_TARG', DEC_TARG, 'RA of target (deg) (J2000)'),
                      end=True)
Primary_header.append(
    ('EQUINOX', EQUINOX, 'equinox of celestial coord. system'), end=True)

Primary_header.append(('', ''), end=True)
Primary_header.append(('', 'INFORMATION ON OTHER ASTROPHYSICAL SOURCES'),
                      end=True)
Primary_header.append(('', ''), end=True)
Primary_header.append(
    ('CANDNUM', -1, 'Number of point-source detections (-1: unknown)'),
    end=True)
Primary_header.append(('DISKDET', True, 'Disk detected in this data set'),
                      end=True)

Primary_header.append(('', ''), end=True)
Primary_header.append(('', 'ASTROMETRIC INFORMATION'), end=True)
Primary_header.append(('', ''), end=True)

Primary_header.append(('PIXSCALE', PIXSCALE, 'pixel scale [arcsec/pixel]'),
                      end=True)

empty_primary = fits.PrimaryHDU(header=Primary_header)

########################################################################################
# creation of DATA_INFORMATION bin table
########################################################################################
c1 = fits.Column(name='Image_Number', array=np.array([1]), format='K')
c2 = fits.Column(name='Orientation',
                 array=np.array([0.]),
                 format='D',
                 unit='deg. e of n')
c3 = fits.Column(name='Combined_rotation_angle',
                 array=np.array([Combined_Rotation_Angle]),
                 format='D',
                 unit='degrees')
c4 = fits.Column(name='Number_of_Exposures',
                 array=np.array([Number_of_Exposures]),
                 format='K')
c5 = fits.Column(name='Exposure_Time',
                 array=np.array([Exposure_Time]),
                 format='D',
                 unit='seconds')

c6 = fits.Column(name='Observation_Start',
                 array=np.array([Observation_Start]),
                 format='D',
                 unit='Mod. Julian Date')
c7 = fits.Column(name='Observation_End',
                 array=np.array([Observation_End]),
                 format='D',
                 unit='Mod. Julian Date')

c8 = fits.Column(name='UT_Midpoint_Date_of_Observaton',
                 array=[UT_Midpoint_Date_of_Observation],
                 format='10A',
                 unit='yyyy-mm-dd')
c9 = fits.Column(name='UT_Midpoint_Time_of_Observaton',
                 array=[UT_Midpoint_Time_of_Observation],
                 format='8A',
                 unit='hh:mm:ss')

c10 = fits.Column(name='Wavelength',
                  array=np.array([Wavelength]),
                  format='D',
                  unit='microns')
c11 = fits.Column(name='Bandwidth',
                  array=np.array([Bandwidth]),
                  format='D',
                  unit='microns')

c12 = fits.Column(name='Polarization', array=[Polarization], format='2A')

bintable_DATA_INFORMATION = fits.BinTableHDU.from_columns(
    [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12],
    name='DATA_INFORMATION')

########################################################################################
# creation of REDUCED_IMAGE image
########################################################################################

header_REDUCED_DATA = fits.Header()

header_REDUCED_DATA.append(('EXTNAME', 'REDUCED_DATA', 'Extension name'),
                           end=True)

header_REDUCED_DATA.append(('BUNIT', 'ADU', 'brightness units'), end=True)
header_REDUCED_DATA.append(('REDALGO', REDALGO, 'Reduction algorithm'),
                           end=True)
header_REDUCED_DATA.append(
    ('REDSTRAT', REDSTRAT, 'Strategy to build the PSF library'), end=True)
header_REDUCED_DATA.append(
    ('EXCLANG', EXCLANG, 'Exclusion angle for ADI-type strategy (deg)'),
    end=True)
header_REDUCED_DATA.append(('TKL', TKL, 'Number of eigenimages used (KLIP)'),
                           end=True)
header_REDUCED_DATA.append(
    ('NPSFLIBR', NPSFLIBR, 'Number of images in the PSF library'), end=True)

image_REDUCED_DATA = fits.ImageHDU(image_reduc, header=header_REDUCED_DATA)

########################################################################################
# creation of SNR_MAP image
########################################################################################

header_SNR_MAP = fits.Header()

header_SNR_MAP.append(('EXTNAME', 'SNR_DATA', 'Extension name'), end=True)

header_SNR_MAP.append(('BUNIT', 'UNITLESS', ' brightness units '), end=True)
header_SNR_MAP.append(
    ('NOISEMET', NOISEMET, 'Method used for the detection limit'), end=True)
image_SNR_MAP = fits.ImageHDU(image_snr, header=header_SNR_MAP)

########################################################################################
# creation of SENSITIVITY_MAP image
########################################################################################

header_SENSITIVITY_MAP = fits.Header()

header_SENSITIVITY_MAP.append(('EXTNAME', 'SENSITIVITY_MAP', 'Extension name'),
                              end=True)
header_SENSITIVITY_MAP.append(('BUNIT', 'ADU', 'brightness units'), end=True)
header_SENSITIVITY_MAP.append(
    ('NOISEMET', NOISEMET, 'Method used for the detection limit'), end=True)
header_SENSITIVITY_MAP.append(
    ('NSIGMA', 1, 'Detection limit confidence level (sigma)'), end=True)

image_SENSITIVITY_MAP = fits.ImageHDU(image_sensitivity_map,
                                      header=header_SENSITIVITY_MAP)

########################################################################################
# creation of the HCI FITS
########################################################################################
hci_fits = fits.HDUList([
    empty_primary, bintable_DATA_INFORMATION, image_REDUCED_DATA,
    image_SNR_MAP, image_SENSITIVITY_MAP
])

name_fits_save = 'HLSP_' +TARGNAME.replace(" ", "") + '_' +  INSTRUME + '_' + PROPOSID + '_' + datecompact + '_' + OBSSTGY + '_' + REDALGO + '_' + REDSTRAT + '.fits'

hci_fits.writeto('/Users/jmazoyer/Desktop/' + name_fits_save, overwrite=True)

#check fits values
# print(hci_fits.info())
# DATA_INFORMATION = hci_fits['REDUCED_DATA']
# print(DATA_INFORMATION.header)