#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example Python script to find the look direction of an imager.
 - Uses spacecraft ephemeris.
 - Plots pointing in 2D to check if the instrument is pointing 
 in Nadir direction.
 - Plots the ground coverage for the imager.
The coverage plotter (fov_plotter) works strictly for Nadir
pointing imagers.
"""

import cavsiopy.attitude_plotter as ap
import cavsiopy.ephemeris_importer as ei
import datetime

path_to_files = '/home/ceren/Documents/USASK/Data/FAI/'

filedate='20201020'
time_start= '043354'
time_end='043846'

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
                             
# import ephemeris from Swarm-E generic .txt files
dict_cas = ei.cas_ephemeris(file_CAS, start_date, end_date)

# import GEO/ITRF position from the sp3 file
dict_sp3 = ei.sp3_ephemeris(file_SP3, start_date, end_date)

# %% coordinates of the center of ICEBEAR radar field-of-view 
pLat = 58
pLon = -106
pAlt = 0.07 # km

# =============================================================================
# %% check if FAI is pointing in Nadir direction. If not, fov_plotter is not
# capable of plotting correct coverage area for off-nadir pointing geometry,yet.
# =============================================================================

# rotate the RRI body vector in orbital frame
body_fai = [0, 0, 1]

# rotate the RRI body vector in orbital frame
rbody = aa.rotate_inst(body_fai, dict_cas['roll'], dict_cas['pitch'], 
                      dict_cas['yaw'])

# RRI body vector in North-East-Center
fai_nec = aa.find_instrument_attitude(rbody, dict_cas['GEIx'], dict_cas['GEIy'], 
                                    dict_cas['GEIz'], dict_cas['GEIVx'], 
                                    dict_cas['GEIVy'], dict_cas['GEIVz'], 
                                    dict_sp3['ITRFx'], dict_sp3['ITRFy'], 
                                    dict_sp3['ITRFz'], dict_sp3['time_experiment'], 
                                    dict_sp3['time_experiment'][0], 
                                    dict_sp3['time_experiment'][-1], 
                                    dict_cas['Lat'], dict_cas['Lon'], 
                                    path_to_sofa_files, method1 ='ephemeris', 
                                    frame2 = 'itrf', frame3 = 'nec')

fai_enc = np.column_stack((fai_nec[:,1], fai_nec[:,0], fai_nec[:,2]))
                                    
# =============================================================================
# %% plotting
# =============================================================================
# Specify altitude, latitude and longitude ranges for plots
# 5 degrees below minimum Latitude
Latmin = np.round(min(dict_cas['Lat'])-2, decimals=0)
# 5 degrees above maximum dict_cas['Lat']itude
Latmax = np.round(max(dict_cas['Lat'])+2, decimals=0)
# 5 degrees below minimum dict_cas['Lon']gitude
Lonmin = np.round(min(dict_cas['Lon'])-2, decimals=0)
if dict_cas['Lon'][0] < pLon:
    # 5 degrees above maximum longitude
    Lonmax = np.round(pLon+7.5, decimals=0)
else:
    # 5 degrees above maximum longitude
    Lonmax = np.round(max(dict_cas['Lon'])+7.5, decimals=0)

extent= [Lonmin, Lonmax, Latmin, Latmax]

# check altitude plots (2D) to see if the spacecraft is pointing along Nadir
ap.attitude_2d_altitude(dict_sp3['time_experiment'], extent, dict_cas['Lon'], 
                        dict_cas['Lat'], dict_cas['Alt'], pLon, fai_enc, 
                        'FAI', 'ICEBEAR', x_axis = 'lon', step = 60)

# plot the FAI field-of-view for Nadir-pointing                      
fig_fov, ax_fov = op.fov_plotter(extent, time_array, Lon, Lat, Alt, fov_deg = 26, 
                pLon, pLat, step = 90, inst_name = 'FAI', target_name = 'ICEBEAR')
