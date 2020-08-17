import matplotlib.pyplot as plt
import glob

import numpy as np

from astropy.io import fits
import pyklip.instruments.GPI as GPI

name_idps = '/Users/jmazoyer/Downloads/IDPS-II/1RXS_J034117.1-225228_2006-09-03_NIRC2_Kp_diva_hlsp.fits'

name_dusty = '/Users/jmazoyer/Downloads/YoungDusty/HIP16449_2009-11-27_NaCo_Lp_diva_hlsp.fits'

hdu = fits.open(name_dusty)
# print(hdu.info())
DATA_INFORMATION = hdu['DATA_INFORMATION']

print(DATA_INFORMATION.data)


raw_data_dir = '/Users/jmazoyer/Dropbox/Work/python/python_data/disk_mcmc/resultats_mcmc_2001010/Aurora/160318_H_Spec/'

filelist_raw = sorted(glob.glob(raw_data_dir + '*spdc_distorcorr.fits'))


parangs = []
integration_time = []
mjd_start = []
ut_time = []


for index_angle in range(0, len(filelist_raw)):

    new_fits = fits.HDUList()
    
    first_header_of_the_raw_fits = fits.getheader(filelist_raw[index_angle], 0)  
    second_header_of_the_raw_fits = fits.getheader(filelist_raw[index_angle], 1)    

    mjd_start.append(first_header_of_the_raw_fits['MJD-OBS'])

    parangs.append(second_header_of_the_raw_fits['AVPARANG'])
    integration_time.append(second_header_of_the_raw_fits['ITIME'])
    ut_time.append(second_header_of_the_raw_fits['EXPSTART'])


Combined_Rotation_Angle = np.nanmax(parangs) - np.nanmin(parangs)
Number_of_Exposures = len(filelist_raw)
Exposure_Time = np.nansum(integration_time)
Observation_Start = mjd_start[0]
Observation_End = mjd_start[-1]
print(Observation_Start,Observation_End)

UT_Midpoint_Date_of_Observation = first_header_of_the_raw_fits['DATE-OBS']
UT_Midpoint_Time_of_Observation = ut_time[len(filelist_raw)//2]

PROPOSID = first_header_of_the_raw_fits['GEMPRGID']
TELESCOP = first_header_of_the_raw_fits['TELESCOP']
INSTRUME = first_header_of_the_raw_fits['INSTRUME']
obs_mode = first_header_of_the_raw_fits['OBSMODE']

FILTER = obs_mode.split('_')[0] 

if FILTER == 'Y':
    Wavelength =  1.04
    Bandwidth = 0.19
elif FILTER == 'J':
    Wavelength =  1.23
    Bandwidth = 0.23
elif FILTER == 'H':
    Wavelength =  1.64
    Bandwidth = 0.30
elif FILTER == 'K1':
    Wavelength =  2.04
    Bandwidth = 0.29
elif FILTER == 'K2':
    Wavelength =  2.26
    Bandwidth = 0.27


TARGNAME = first_header_of_the_raw_fits['OBJECT']
RA_TARG = first_header_of_the_raw_fits['RA']
DEC_TARG = first_header_of_the_raw_fits['DEC']
EQUINOX = first_header_of_the_raw_fits['EQUINOX']
