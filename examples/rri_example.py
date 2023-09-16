#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example python script for determining and plotting the pointing 
direction of RRI.
"""

import cavsiopy.ephemeris_importer as ei
import cavsiopy.attitude_analysis as aa
import cavsiopy.attitude_plotter as ap
import cavsiopy.use_rotation_matrices as urm


# locate RRI data file
path_to_data = '/home/ceren/Documents/USASK/Data/2015-2018/'
name_of_RRI_file= "RRI_20160418_222759_223156_lv1_12.0.0.h5"

file_RRI = path_to_data + name_of_RRI_file

# point to where the IERS and EOP files are
path_to_sofa_files = '/home/ceren/Documents/GitHub/cavsiopy-pysofa2/cavsiopy/'

# import RRI ephemeris
dict_rri = ei.rri_ephemeris(file_RRI)

# define RRI body vector
body_rri = [1, 0 , 0]

# define initial dipole positions
a = 1

# rotate the RRI body vector in orbital frame
rbody = aa.rotate_inst(body_rri, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])

# RRI body vector in North-East-Center
rri_nec = aa.find_instrument_attitude(rbody, dict_rri['GEIx'],dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files, 
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'ned')
rri_enc = urm.ned2enu(rri_nec).T

# ----------------------------------------------------------------------------
body_m1 = [0, -a, -a]
r_m1 = aa.rotate_inst(body_m1, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
m1_nec = aa.find_instrument_attitude(r_m1,dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files,
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'nec')
m1_enc = urm.ned2enu(m1_nec).T 
# ----------------------------------------------------------------------------
body_m2 = [0, a, a]  # dipole1
r_m2 = aa.rotate_inst(body_m2, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
m2_nec = aa.find_instrument_attitude(r_m2,dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files,
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'nec')
m2_enc = urm.ned2enu(m2_nec).T
# ----------------------------------------------------------------------------
body_m3 = [0, a, -a]
r_m3 = aa.rotate_inst(body_m3, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
m3_nec = aa.find_instrument_attitude(r_m3,dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files,
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'nec')
m3_enc = urm.ned2enu(m3_nec).T 
# ---------------------------------------------------------------------------- 
body_m4 = [0, -a, a]  # dipole2
r_m4 = aa.rotate_inst(body_m4, dict_rri['roll'], dict_rri['pitch'], 
                      dict_rri['yaw'])
m4_nec = aa.find_instrument_attitude(r_m4,dict_rri['GEIx'], dict_rri['GEIy'], 
                                    dict_rri['GEIz'], dict_rri['GEIVx'], 
                                    dict_rri['GEIVy'], dict_rri['GEIVz'], 
                                    dict_rri['GEOx'], dict_rri['GEOy'], 
                                    dict_rri['GEOz'], dict_rri['time_array'], 
                                    dict_rri['start_time'], 
                                    dict_rri['time_array'][-1], 
                                    dict_rri['Lat'], dict_rri['Lon'], 
                                    path_to_sofa_files,
                                    method1='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'nec')
m4_enc = urm.ned2enu(m4_nec).T
# ----------------------------------------------------------------------------
# Ottawa coordinates
pLat = 45.40717
pLon = -75.55394
pAlt = 0.07 # Ottawa, km
target_name = 'Ottawa'

extent = [dict_rri['Lon'].min()-10, dict_rri['Lon'].max()+10,
          dict_rri['Lat'].min()-10, dict_rri['Lat'].max()+10]
# %% 
vec_args= m1_enc, m3_enc, m2_enc, m4_enc, rri_enc
# vec_legend=['$Dipole_{1}$','$Dipole_{2}$', 'Boresight']
# arrow_lengths=[.5, .5, .5, .5, .75]
# vec_linestyles=['solid', 'solid', 'solid', 'solid', 'solid']
# vec_colors=['red', '#12e193', 'red', '#12e193', 'black']
# vec_linewidth = [1.5, 1.5, 1.5, 1.5, 1.5]
# arrowhead = [0.01, 0.01, 0.01, 0.01,  0.45]
# legend_colors=['red', '#12e193', 'black']
# legend_edgecolors=['red', '#12e193', 'black']
# legend_markers=['_','_',  r'$\longrightarrow$']

title_vec = dict_rri['time_array'][0].strftime('%d %B %Y')
# connected_plot = ap.display_observation_geometry(title_vec, 
#                                                   dict_rri['time_array'],
#                                                   pLon, pLat, OH, 
#                                                   dict_rri['Lon'], 
#                                                   dict_rri['Lat'],
#                                                   dict_rri['Alt'] ,
#                                                   target_name,
#                                                   *vec_args)
quiver_plot = ap.attitude_3d_ground_quiver(title_vec, dict_rri['time_array'],
                                        pLon, pLat, pAlt, dict_rri['Lon'], 
                                        dict_rri['Lat'], dict_rri['Alt'] , 
                                        rri_enc, target_name, *vec_args)


# %%
# ap.attitude_2d_on_map(dict_rri['time_array'], extent, dict_rri['Lon'], 
#                       dict_rri['Lat'], dict_rri['Alt'], pLon, pLat, rri_enc, \
#                               target_name, 'RRI', 158, step = 60)
