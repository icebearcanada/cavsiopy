#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 22:50:17 2023

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

import orientation_plotter as op
import ephemeris_importer as ei
import Direction_Analysis as da
import use_rotation_matrices as RM
import misc

path_to_files = '/home/ceren/Documents/USASK/Data/FAI/'

# filedate='20191027'
# time_start= '040300'
# time_end='041300'

# filedate='20180310'
# time_start= '052100'
# time_end='052600'

filedate='20200319'
time_start= '090000'
time_end='090400'


start_time = int(time_start)
end_time = int(time_end)
start_hour=int(time_start[0:2])
start_min=int(time_start[2:4])
start_sec=int(time_start[4:6])
start_year=int(filedate[0:4])
start_month=int(filedate[4:6])
start_day=int(filedate[6:8])
end_hour=int(time_end[0:2])
end_min=int(time_end[2:4])
end_sec=int(time_end[4:6])


start_date= datetime.datetime(start_year, start_month, start_day,
                               start_hour, start_min, start_sec)
end_date= datetime.datetime(start_year, start_month, start_day,
                             end_hour, end_min, end_sec+1)

# %% 

name_of_file=''.join(['CAS_ephemeris_',filedate,'T000000_',filedate,'T235959_1.1.0'])
file_CAS= path_to_files + name_of_file +'.txt'

srow, erow, start, end, time_array, Lat, Lon, Alt, \
    GEIx_fai, GEIy_fai, GEIz_fai, GEIVx_fai, GEIVy_fai, GEIVz_fai, \
        roll, pitch, yaw, accuracy = \
            ei.cas_detailed_ephemeris(file_CAS, time_start, time_end)

# %% import sp3 file for GEO/ITRF position
name_of_ephemeris_file=''.join(['CAS_Orbit_GEO_',filedate,'T000000_',\
                    filedate,'T235959_1.1.0'])
file_SP3 = path_to_files + name_of_ephemeris_file +'.sp3'

time_array_sp3, GEOx_fai, GEOy_fai, GEOz_fai = \
    ei.sp3_ephemeris_short(file_SP3, start_date, end_date)

# %% Ottawa coordinates
pLat = 57
pLon = -107
pAlt = 0.07 # Ottawa, km
P_H = 6371.2 + 0.07;

body_fai = [0, 0, 1]

GEO_fai = [GEOx_fai, GEOy_fai, GEOz_fai]

FAI_ITRF, j2k_itrf_enc, j2k_itrf_nec = \
    da.inst_direction_in_nec_using_icrf_to_itrf(path_to_files, filedate,
                                      time_array_sp3, start_date, end_date,\
                                          GEIx_fai, GEIy_fai, GEIz_fai, \
                                              GEIVx_fai, GEIVy_fai, GEIVz_fai,\
                                              GEO_fai, roll, pitch, yaw, body_fai)
 
j2k_itrf_ned, j2k_itrf_enu = \
    da.inst_direction_in_ned_from_itrf(Lat, Lon, FAI_ITRF)

# =============================================================================
# %%
# =============================================================================
# Specify altitude, latitude and longitude ranges for plots
# 5 degrees below minimum Latitude
Latmin = np.round(min(Lat)-2, decimals=0)
# 5 degrees above maximum Latitude
Latmax = np.round(max(Lat)+2, decimals=0)
# 5 degrees below minimum Longitude
Lonmin = np.round(min(Lon)-2, decimals=0)
if Lon[0] < pLon:
    # 5 degrees above maximum Longitude
    Lonmax = np.round(pLon+7.5, decimals=0)
else:
    # 5 degrees above maximum Longitude
    Lonmax = np.round(max(Lon)+7.5, decimals=0)


sat_elev, sat_azim = da.LA_sat(pLat, pLon, Lat, Lon, P_H, Alt)
start_time = time_array_sp3[0] 
extent= [Lonmin, Lonmax, Latmin, Latmax]
central_lon, central_lat, left, width, bottom, height = op.coverage(extent)
fov=op.fov_plotter(time_array_sp3, extent, Alt, Lat, Lon, 20, central_lon,\
                    central_lat, pLon, pLat, 90, j2k_itrf_enu, 
                    sat_azim, sat_elev, start_date.strftime('%Y-%m-%d'), 
                    start_time, 'FAI', 'ICEBEAR', 10)