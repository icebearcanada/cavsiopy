#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 11:35:16 2023

@author: ceren
"""

import ephemeris_importer as ei
import attitude_analysis as aa
import attitude_plotter as ap
import use_rotation_matrices as urm

# locate RRI data file
path_to_data = '/home/ceren/Documents/USASK/Data/2015-2018/'
name_of_RRI_file= "RRI_20160418_222759_223156_lv1_12.0.0.h5"

file_RRI = path_to_data + name_of_RRI_file

# point to where the IERS and EOP files are
path_to_sofa_files = '/home/ceren/OneDrive/cavsiopy_1.0.0/cavsiopy/'

# import RRI ephemeris
dict_rri = ei.rri_ephemeris(file_RRI)

# define RRI body vector
body_rri = [1, 0 , 0]

# define initial dipole positions
a = 1
body_m1 = [0, -a, -a]
body_m2 = [0, a, a]  # dipole1
body_m3 = [0, a, -a]
body_m4 = [0, -a, a]  # dipole2

# rotate the RRI body vector in orbital frame
rbody = aa.rotate_rri(body_rri, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
r_m1 = aa.rotate_rri(body_m1, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
r_m2 = aa.rotate_rri(body_m2, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
r_m3 = aa.rotate_rri(body_m3, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
r_m4 = aa.rotate_rri(body_m4, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])

# RRI body vector in North-East-Center
rri_nec = aa.find_instrument_attitude(rbody, dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVx'], dict_rri['GEIVx'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files)

m1_nec = aa.find_instrument_attitude(r_m1, dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVx'], dict_rri['GEIVx'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files)

m2_nec = aa.find_instrument_attitude(r_m2, dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVx'], dict_rri['GEIVx'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files)


m3_nec = aa.find_instrument_attitude(r_m3, dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVx'], dict_rri['GEIVx'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files)

m4_nec = aa.find_instrument_attitude(r_m4, dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVx'], dict_rri['GEIVx'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files)

m1_enu = urm.ned2enu(m1_nec).T ; m2_enu = urm.ned2enu(m2_nec).T
m3_enu = urm.ned2enu(m3_nec).T ; m4_enu = urm.ned2enu(m4_nec).T
rri_enu = urm.ned2enu(rri_nec).T

# Ottawa coordinates
pLat = 45.40717
pLon = -75.55394
pAlt = 0.07 # Ottawa, km
target_name = 'Ottawa'

# extent = [dict_rri['Lon'].min(), dict_rri['Lon'].max,
#           dict_rri['Lat'].min(), dict_rri['Lat'].max]
# %% 
arrow_lengths=[.5, .5, .5, .5, .75]
vec_linestyles=['solid', 'solid', 'solid', 'solid', 'solid']
vec_colors=['red', '#12e193', 'red', '#12e193', 'black']
vec_linewidth = [1.5, 1.5, 1.5, 1.5, 1.5]

arrowhead = [0.01, 0.01, 0.01, 0.01,  0.45]

legend_colors=['red', '#12e193', 'black']

legend_edgecolors=['red', '#12e193', 'black']
legend_markers=['_','_',  r'$\longrightarrow$']

OH = 0

sct_params = {'alpha': 1, 'edgecolor': 'black', 'c': 'lightgrey', 'marker': '*',
              's': 180}

vec_args= m1_enu, m3_enu, m2_enu, m4_enu, rri_enu
vec_legend=['$Dipole_{1}$','$Dipole_{2}$', 'Boresight']


title_vec = dict_rri['time_array'][0].strftime('%d %B %Y')
connected_plot = ap.display_observation_geometry(title_vec, 
                                                 dict_rri['time_array'], 
                                                 len(dict_rri['time_array']), 
                                                 60,
                            pLon, pLat, OH, dict_rri['Lon'], dict_rri['Lat'],
                            dict_rri['Alt']/111,
                            arrow_lengths, 'upper center', target_name, 
                            *vec_args,
                            vector_colors = vec_colors,
                            linestyles = vec_linestyles,
                            linewidth = vec_linewidth,
                            labels = vec_legend,
                            markers = legend_markers,
                            colors = legend_colors,
                            edgecolors = legend_edgecolors,
                            arrowhead = arrowhead,
                            sct_kwargs = sct_params)