#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:16:06 2023

@author: ceren
"""

import glob
import h5py
import numpy as np
import datetime
from astropy.coordinates import cartesian_to_spherical
import matplotlib.pyplot as plt
import matplotlib.dates as dtformat
from matplotlib.ticker import MultipleLocator

import cavsiopy.orientation_plotter as op
import cavsiopy.ephemeris_importer as ei
import cavsiopy.Direction_Analysis as da
import cavsiopy.use_rotation_matrices as RM
import cavsiopy.misc

path_to_files = '/home/ceren/Documents/USASK/Data/2015-2018/'

# name_of_RRI_file="RRI_20171221_220914_221511_lv1_12.0.0.h5" # slewed
# name_of_RRI_file='RRI_20151127_115523_115750_lv1_12.0.0.h5'
name_of_RRI_file="RRI_20171213_234044_234641_lv1_12.0.0.h5" # ram
# name_of_RRI_file= "RRI_20171220_223344_223941_lv1_12.0.0.h5" # nadir
# name_of_RRI_file ="RRI_20170818_152114_152511_lv1_12.0.0.h5"
# name_of_RRI_file ='RRI_20170512_174914_175311_lv1_12.0.0.h5'
# name_of_RRI_file= "RRI_20160418_222759_223156_lv1_12.0.0.h5"
# name_of_RRI_file="RRI_20160419_220939_221336_lv1_12.0.0.h5"
# name_of_RRI_file="RRI_20160420_215117_215514_lv1_12.0.0.h5"
# name_of_RRI_file= "RRI_20160421_213255_213652_lv1_12.0.0.h5"
# name_of_RRI_file="RRI_20160422_211435_211832_lv1_12.0.0.h5"
# name_of_RRI_file="RRI_20160115_103314_103511_lv1_12.0.0.h5"
# name_of_RRI_file="RRI_20170520_163244_163641_lv1_12.0.0.h5"
# name_of_RRI_file="RRI_20170513_172714_173041_lv1_12.0.0.h5" # ascending, slewed, 6.9285

filedate = str(name_of_RRI_file[4:12])
version = str(name_of_RRI_file[27:37])
# complete file location
file_RRI = path_to_files + name_of_RRI_file
# Import RRI data of selected file
file = h5py.File(file_RRI, 'r')

time_data = \
    file.get(
        'CASSIOPE Ephemeris/Ephemeris MET (seconds since May 24, 1968)')[:]

time_data_rri, datetime_time_array_1sec_rri, time_array_1sec_rri, \
    start_time_rri, end_time_rri, Lat_rri, Lon_rri, Alt_rri, \
        MLat_rri, MLon_rri, MLT_rri, \
        GEOx_rri, GEOy_rri, GEOz_rri, GEOVx_rri, GEOVy_rri, GEOVz_rri,\
        GEIx_rri, GEIy_rri, GEIz_rri, GEIVx_rri, GEIVy_rri, GEIVz_rri, \
    roll_rri, pitch_rri, yaw_rri = ei.import_rri_ephemeris(file_RRI)

date_str = '-'.join([filedate[6:8], filedate[4:6], filedate[0:4]])

# Ottawa coordinates
pLat = 45.4215
pLon = -75.6972
OH = 0  # Ottawa, km

# =============================================================================
# %% Rotation matrices J2000
# =============================================================================
# %% time calculations
start_year = int(name_of_RRI_file[4:8])
start_month = int(name_of_RRI_file[8:10])
start_day = int(name_of_RRI_file[10:12])
start_hour = int(name_of_RRI_file[13:15])
start_min = int(name_of_RRI_file[15:17])
start_sec = int(name_of_RRI_file[17:19])
end_hour = int(name_of_RRI_file[20:22])
end_min = int(name_of_RRI_file[22:24])
end_sec = int(name_of_RRI_file[24:26])

start_date= datetime.datetime(start_year, start_month, start_day,
                               start_hour, start_min, start_sec)
end_date= datetime.datetime(start_year, start_month, start_day,
                             end_hour, end_min, end_sec+1)

GEO_rri = [GEOx_rri, GEOy_rri, GEOz_rri]

RRI_ITRF, j2k_itrf_enc, j2k_itrf_nec = \
    da.inst_direction_in_nec_using_icrf_to_itrf(path_to_files, filedate,
                                      time_array_1sec_rri, start_date, end_date,\
                                          GEIx_rri, GEIy_rri, GEIz_rri, \
                                              GEIVx_rri, GEIVy_rri, GEIVz_rri,\
                                              time_data_rri, GEO_rri, \
                                                  roll_rri, pitch_rri, yaw_rri)

# =============================================================================
# %% plotters
# =============================================================================
import matplotlib.gridspec as gridspec
myFmt = dtformat.DateFormatter("%H:%M:%S")
date = start_date.strftime('%d-%b-%Y')
seconds = dtformat.MicrosecondLocator(interval=10000000)

end_ind, number_of_int = op.indices_and_intervals(
    start_date, time_data, 25)

# Specify altitude, latitude and longitude ranges for plots
# 5 degrees below minimum Latitude
Latmin = np.round(min(Lat_rri)-2, decimals=0)
# 5 degrees above maximum Latitude
Latmax = np.round(max(Lat_rri)+2, decimals=0)
# 5 degrees below minimum Longitude
Lonmin = np.round(min(Lon_rri)-2, decimals=0)
if Lon_rri[0] < pLon:
    # 5 degrees above maximum Longitude
    Lonmax = np.round(pLon+7.5, decimals=0)
else:
    # 5 degrees above maximum Longitude
    Lonmax = np.round(max(Lon_rri)+7.5, decimals=0)
# Altmin = np.round(min(Alt)-10, decimals=0)  # 10 km below minimum Altitude
# Altmax = np.round(max(Alt)+10, decimals=0)  # 10 km above maximum altitude

# set the map coverage:
extent = [Lonmin, Lonmax, Latmin, Latmax]
central_lon, central_lat, left, width, bottom, height = op.coverage(extent)

fig = plt.figure(figsize=(16, 10))
spec = gridspec.GridSpec(ncols=1, nrows=1, figure = fig)
plt.suptitle(date, fontsize = 14, y = 0.99)

title_ITRF = ' '.join(['T4: ORF->ICRF->ITRF->NEC' ])
op.plotter_with_map(fig, spec, extent, Lon_rri, Lat_rri, j2k_itrf_enc, 0, 1,
                      time_array_1sec_rri, pLon, pLat, end_ind, 25, -3,
                  'Geographic Longitude', 'Geographic Latitude (km)',
                  title_ITRF, 'bottom_on','off', 0)
