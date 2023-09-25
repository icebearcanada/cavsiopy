#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 01:17:42 2023

@author: ceren
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dtformat
from matplotlib.ticker import MultipleLocator

import cavsiopy.attitude_plotter as ap
import cavsiopy.attitude_analysis as aa
import cavsiopy.ephemeris_importer as ei

# Ottawa coordinates
pLat = 45.40717
pLon = -75.55394
pAlt = 0.07 # Ottawa, km
target_name = 'Ottawa'

# ----------------------------------------------------------------------------
# locate RRI data file
path_to_data = 'path/to/data/files/'
name_of_RRI_file= "RRI_20160418_222759_223156_lv1_12.0.0.h5"

file_RRI = path_to_data + name_of_RRI_file

# import RRI ephemeris
dict_rri = ei.rri_ephemeris(file_RRI)

# define body frame unit vectors
a = 1
# define RRI body vector
body_rri = [a, 0 , 0]

# rotate the RRI body vector in orbital frame
rbody = aa.rotate_inst(body_rri, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])

print('finding rri_ned')
# RRI body vector in North-East-Down
rri_ned = aa.find_instrument_attitude(rbody, dict_rri['GEIx'],dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'ned')

# ----------------------------------------------------------------------------
print('finding dipole1_ned edge1')
body_m1 = [0, -a, -a] # dipole1 - edge1
r_m1 = aa.rotate_inst(body_m1, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
m1_ned = aa.find_instrument_attitude(r_m1,dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'ned')

# ----------------------------------------------------------------------------
print('finding dipole2_ned edge1')
body_m3 = [0, a, -a] # dipole2 - edge1
r_m3 = aa.rotate_inst(body_m3, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
m3_ned = aa.find_instrument_attitude(r_m3,dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'ned')

# ----------------------------------------------------------------------------
print('finding spacecraft distance to target')
distance = aa.spacecraft_distance_from_a_point(pLon, pLat, pAlt, 
                                               dict_rri['Lon'], 
                                               dict_rri['Lat'], 
                                               dict_rri['Alt'])

# proximity to target--------------------------------------------------------
# find where range between the satellite and target is minimum
min_index = np.where(distance == np.min(distance))[0][0]
               
# ---------------- -----------------------------------------------------------
print('calculating the reception and complementary reception angles')

ra_rri, sra_rri = aa.calculate_reception_angle(rri_ned, pLat, pLon, pAlt, 
                                      dict_rri['Lat'], dict_rri['Lon'], 
                                      dict_rri['Alt'], inst = 'boresight')

ra_d1, cra_d1 = aa.calculate_reception_angle(m1_ned, pLat, pLon, pAlt, 
                                      dict_rri['Lat'], dict_rri['Lon'], 
                                      dict_rri['Alt'], inst = 'D1')

ra_d2, cra_d2 = aa.calculate_reception_angle(m3_ned, pLat, pLon, pAlt, 
                                      dict_rri['Lat'], dict_rri['Lon'], 
                                      dict_rri['Alt'], inst = 'D2')
# ----------------------------------------------------------------------------
print('calculating the slew conditions')
sa1 = 5; sa2 = 10; sa3 = 15; sa4 = 20 ; sa5 = 25; sa6 = 30
slew_1 = aa.find_slew_rri(ra_rri, ra_d1, ra_d2, sa1)
slew_2 = aa.find_slew_rri(ra_rri, ra_d1, ra_d2, sa2)
slew_3 = aa.find_slew_rri(ra_rri, ra_d1, ra_d2, sa3)
slew_4 = aa.find_slew_rri(ra_rri, ra_d1, ra_d2, sa4)
slew_5 = aa.find_slew_rri(ra_rri, ra_d1, ra_d2, sa5)
slew_6 = aa.find_slew_rri(ra_rri, ra_d1, ra_d2, sa6)

# %%---------------------------------------------------------------------------
print('plot complementary reception and slew conditions')
date = dict_rri['time_array'][0].strftime('%d-%b-%Y')
myFmt = dtformat.DateFormatter("%H:%M:%S")
seconds = dtformat.MicrosecondLocator(interval=10000000)
ylimits = -10, 190
step = 10

fig_slew_angles, (ax_reception, ax_slew) = plt.subplots(figsize=(12.06, 10.79),
                                                        nrows = 2, 
                                                        sharex = True)
plt.suptitle(date, fontsize = 18)

# panel 1- complementary reception angles
ax_reception.plot(dict_rri['time_array'][::step], sra_rri[::step], c='k', 
                  marker = 'o', markersize = '5')
ax_reception.plot(dict_rri['time_array'][::step], cra_d1[::step], c='r', 
                  marker = 'o', markersize = '5')
ax_reception.plot(dict_rri['time_array'][::step], cra_d2[::step], c='#04d8b2', 
                  marker = 'o', markersize = '5')

# Mark closest approach to Ottawa
ax_reception.scatter(dict_rri['time_array'][min_index], ylimits[0]+5, 
                     marker='*',  c='k', edgecolor='k', s=50, zorder= 50)

# axes properties
ax_reception.set_ylim(ylimits)
ax_reception.yaxis.set_major_locator(MultipleLocator(45))
ax_reception.yaxis.set_minor_locator(MultipleLocator(15))
ax_reception.set_ylabel('Complementary\n angles (\N{DEGREE SIGN})',fontsize=14)
ax_reception.grid(True, linestyle='--')
ax_reception.tick_params(axis='both', which='major', labelsize=14)

# panel 2- plot slew conditions
ylim_min, ylim_max = 0, 35
ap.plot_slew_rri(ax_slew, ylim_min, ylim_max, 5, dict_rri['time_array'], 
             slew_1, str(sa1))
ap.plot_slew_rri(ax_slew, ylim_min, ylim_max, 4, dict_rri['time_array'], 
             slew_2, str(sa2))
ap.plot_slew_rri(ax_slew, ylim_min, ylim_max, 3, dict_rri['time_array'], 
             slew_3, str(sa3))
ap.plot_slew_rri(ax_slew, ylim_min, ylim_max, 2, dict_rri['time_array'], 
             slew_4, str(sa4))
ap.plot_slew_rri(ax_slew, ylim_min, ylim_max, 1, dict_rri['time_array'], 
             slew_5, str(sa5))
ap.plot_slew_rri(ax_slew, ylim_min, ylim_max, 0, dict_rri['time_array'], 
             slew_6, str(sa6), cb_axis = 'yes', time='yes')
ax_slew.set_ylabel('Closeness to Slew', fontsize = 18)
ax_slew.get_yaxis().set_label_coords(-0.075,0.5)
ax_slew.yaxis.set_tick_params(labelleft=False)
ax_slew.tick_params(axis='both', which='major', labelsize=14)
