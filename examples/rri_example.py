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

time_data_rri, time_array_1sec_rri, time_array_1sec_rri, start_time, end_time,\
    Lat_rri, Lon_rri, Alt_rri, MLat_rri, MLon_rri, MLT_rri,\
        GEOx_rri, GEOy_rri, GEOz_rri, GEOVx_rri, GEOVy_rri, GEOVz_rri, \
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

rri_j2k_itrf_nec = \
    da.find_inst_attitude(path_to_files, filedate, time_array_1sec_rri, 
                                                start_date, end_date,
                                                GEIx_rri, GEIy_rri, GEIz_rri, 
                                                GEIVx_rri, GEIVy_rri, GEIVz_rri, 
                                                GEOx_rri, GEOy_rri, GEOz_rri, 
                                                roll, pitch, yaw, 
                                                ‘icrf’, ‘itrf’, ‘nec’)

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

# =============================================================================
# %%
# =============================================================================
# =============================================================================
# %% plot 3D pointing with
# =============================================================================
sct_params = {'alpha': 1, 'edgecolor': 'black', 'c': 'lightgrey', 'marker': '*',
              's': 180}
# =============================================================================
# %% plot trajectory
# =============================================================================
# arrow_lengths=[0.5, 0.5, 0.5, 0.5, 1, 1]
# vec_linestyles=['dashed', 'solid', 'dashed', 'solid', 'solid']
# vec_colors=['black', 'black', 'black' , 'black', 'red']
# legend_colors=['black','black','red','None', 'black']
# legend_edgecolors=['black', 'black','red','black', 'black']
# legend_markers=['_','$--$', r'$\longrightarrow$','H','H']
# OH = 0
# target_name = 'Ottawa'
# if ray_model_input==0:
#     vec_args= M1_enu[1:SizeArr], M3_enu[1:SizeArr], M2_enu[1:SizeArr], \
#         M4_enu[1:SizeArr], RRI_enu[1:SizeArr]

#     vec_legend=['Dipole_1','Dipole_2','Boresight',
#         'North of target','South of target']
# else:
#     vec_args= M1_enu[1:SizeArr], M3_enu[1:SizeArr], M2_enu[1:SizeArr], \
#         M4_enu[1:SizeArr], RRI_enu[1:SizeArr]
#     vec_legend=['Dipole_1','Dipole_2', 'Boresight',
#                 'North of target','South of target']

# op.trajectory_3D(target_name, start_time,time_array_1sec, end_ind, 60,
#                             pLon, pLat, OH, Lon, Lat, Alt/100, RRI_enu,
#                             arrow_lengths, 'upper center',\
#                             *vec_args, vector_colors=vec_colors,
#                             labels=vec_legend, linestyles=vec_linestyles,
#                             markers=legend_markers, colors=legend_colors,
#                             edgecolors=legend_edgecolors)

# =============================================================================
# %% plot with connected points-with ray
# =============================================================================

# arrow_lengths=[0.5, 0.5, 0.5, 0.5, 1, 1]
# vec_linestyles=['solid', 'solid', 'solid', 'solid', 'solid',\
#                 'solid']
# vec_colors=['red', '#12e193', 'red', '#12e193', 'black', \
#               'blue']
# vec_linewidth = [2.5, 2.5, 2.5, 2.5, 2, 2]

# arrowhead = [0.01, 0.01, 0.01, 0.01,  0.25, 0.25]

# legend_colors=['red', '#12e193', 'black','blue',\
#                 'None', 'black']

# legend_edgecolors=['red', '#12e193', 'black','blue',\
#                     'black', 'black']
# legend_markers=['_','_',  \
#                 r'$\longrightarrow$', r'$\longrightarrow$', \
#                 'H','H']

# OH = 0
# if ray_model_input==0:
#     vec_args= M1_enu[1:SizeArr], M3_enu[1:SizeArr], M2_enu[1:SizeArr], \
#         M4_enu[1:SizeArr], RRI_enu[1:SizeArr],  los_enu_arr

#     vec_legend=['Dipole_1','Dipole_2','Boresight', 'Ray',
#         'N','S']
# else:
#     vec_args= M1_enu[1:SizeArr], M3_enu[1:SizeArr], M2_enu[1:SizeArr], \
#         M4_enu[1:SizeArr], RRI_enu[1:SizeArr], los_enu_arr
#     vec_legend=['Dipole_1','Dipole_2', 'Boresight','Ray',
#                 'N','S']


# title_vec=time_array_1sec[0].strftime('%d %B %Y')
# connected_plot = \
#     op.vector_direction_plotter_connected_ground_trajectory_3D(title_vec, \
#                                                     time_array_1sec, SizeArr, 60,
#                             pLon, pLat, OH, Lon, Lat, Alt/111,
#                             arrow_lengths, 'upper center', target_name, *vec_args,
#                             vector_colors = vec_colors,
#                             linestyles = vec_linestyles,
#                             linewidth = vec_linewidth,
#                             labels = vec_legend,
#                             markers = legend_markers,
#                             colors = legend_colors,
#                             edgecolors = legend_edgecolors,
#                             arrowhead = arrowhead,
#                             sct_kwargs = sct_params)

# connected_plot.savefig('example_RRItrajectory_experiment2.png', format='png', dpi = 600)
# =============================================================================
# %% plot with connected points-without ray
# =============================================================================

# arrow_lengths=[0.5, 0.5, 0.5, 0.5, 1.5]
# vec_linestyles=['solid', 'solid', 'solid', 'solid', 'solid']
# vec_colors=['red', '#12e193', 'red', '#12e193', 'black']
# vec_linewidth = [2.5, 2.5, 2.5, 2.5, 2]

# arrowhead = [0.01, 0.01, 0.01, 0.01,  0.45]

# legend_colors=['red', '#12e193', 'black','None', 'black']

# legend_edgecolors=['red', '#12e193', 'black', 'black', 'black']
# legend_markers=['_','_',  r'$\longrightarrow$','H','H']

# OH = 0
# if ray_model_input==0:
#     vec_args= M1_enu[1:SizeArr], M3_enu[1:SizeArr], M2_enu[1:SizeArr], \
#         M4_enu[1:SizeArr], RRI_enu[1:SizeArr]

#     vec_legend=['Dipole_1','Dipole_2','Boresight', 'N','S']
# else:
#     vec_args= M1_enu[1:SizeArr], M3_enu[1:SizeArr], M2_enu[1:SizeArr], \
#         M4_enu[1:SizeArr], RRI_enu[1:SizeArr]
#     vec_legend=['Dipole_1','Dipole_2', 'Boresight', 'N','S']


# title_vec=time_array_1sec[0].strftime('%d %B %Y')
# connected_plot = op.vector_direction_plotter_connect_to_point(title_vec, time_array_1sec, SizeArr, 30,
#                             pLon, pLat, OH, Lon, Lat, Alt/111,
#                             arrow_lengths, 'upper center', target_name, *vec_args,
#                             vector_colors = vec_colors,
#                             linestyles = vec_linestyles,
#                             linewidth = vec_linewidth,
#                             labels = vec_legend,
#                             markers = legend_markers,
#                             colors = legend_colors,
#                             edgecolors = legend_edgecolors,
#                             arrowhead = arrowhead,
#                             sct_kwargs = sct_params)

# plot_name = '_'.join([start_time.strftime('%d-%b-%Y') , '3Dview.png'])

# connected_plot.savefig(plot_name, format='png', dpi = 400)

# %% ========================================================================
#  plot one point at a specific time
arrow_lengths=[0.5, 0.5, 0.5, 0.5, 1, 1]
vec_linestyles=['solid', 'solid', 'solid', 'solid', 'solid', 'solid']
vec_colors=['red', '#12e193', 'red', '#12e193', 'black', 'blue']
vec_linewidth = [3.5, 3.5, 3.5, 3.5, 3, 3]
arrowhead = [0.015, 0.015, 0.015, 0.015,  0.3, 0.3]
legend_colors=['red', '#12e193', 'black','blue', 'None', 'black']
legend_edgecolors=['red', '#12e193', 'black','blue', 'black', 'black']
legend_markers=['_','_',  r'$\longrightarrow$', r'$\longrightarrow$', 'H','H']
# s1 = 176; s2 = 185
# s1 = 105; s2 = 114
s1 = 130; s2 = 140
OH = 0
if ray_model_input==0:
    vec_args= M1_enu[s1:s2], M3_enu[s1:s2], M2_enu[s1:s2], \
        M4_enu[s1:s2], RRI_enu[s1:s2],  los_enu_arr[s1:s2]

    vec_legend=['Dipole_1','Dipole_2','Boresight', 'Ray',
        'N','S']
else:
    vec_args= M1_enu[s1:s2], M3_enu[s1:s2], M2_enu[s1:s2], \
        M4_enu[s1:s2], RRI_enu[s1:s2], los_enu_arr[s1:s2]
    vec_legend=['Dipole_1','Dipole_2', 'Boresight','Ray',
                'N','S']


title_vec=time_array_1sec[0].strftime('%d %B %Y')
connected_plot = \
    op.vector_direction_plotter_connect_to_target(title_vec, 
                            time_array_1sec[s1:s2], SizeArr, 20, pLon, pLat, OH, 
                            Lon[s1:s2], Lat[s1:s2], Alt[s1:s2]/111,
                            arrow_lengths, 'upper center', target_name, *vec_args,
                            vector_colors = vec_colors,
                            linestyles = vec_linestyles,
                            linewidth = vec_linewidth,
                            labels = vec_legend,
                            markers = legend_markers,
                            colors = legend_colors,
                            edgecolors = legend_edgecolors,
                            arrowhead = arrowhead,
                            sct_kwargs = sct_params)
# plot_name = '_'.join([start_time.strftime('%d-%b-%Y') , 'onepoint.png'])
# connected_plot.savefig(plot_name, format='png', dpi = 400)