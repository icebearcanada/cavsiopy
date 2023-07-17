#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 01:17:42 2023

@author: ceren
"""


import matplotlib.pyplot as plt
import attitude_plotter as ap
import misc
import attitude_analysis as aa


              
        
# reception angles ----------------------------------------------
offset_D1 = 90-obs_m1 # complementary reception angle 1
offset_D2 = 90-obs_m3 # complementary reception angle 2
offset_los = 180-obs_los # reception angle for RRI


# find slew mode---------------- ----------------------------------------------
print('calculating the slew parameter')
sa1 = 5; sa2 = 10; sa3 = 15; sa4 = 20 ; sa5 = 25; sa6 = 30
slew_1 = aa.find_slew(offset_los, offset_D1, offset_D2, sa1)
slew_2 = aa.find_slew(offset_los, offset_D1, offset_D2, sa2)
slew_3 = aa.find_slew(offset_los, offset_D1, offset_D2, sa3)
slew_4 = aa.find_slew(offset_los, offset_D1, offset_D2, sa4)
slew_5 = aa.find_slew(offset_los, offset_D1, offset_D2, sa5)
slew_6 = aa.find_slew(offset_los, offset_D1, offset_D2, sa6)

# pass-type-------------------------------------------------------------------
RRIz = abs(np.nanmean(RRI_enu[:,2]))
offset_avg = np.nanmean(offset_los_90)
ang = 25
if offset_avg < ang:
    pass_type = 'Slewed pass'
elif 0.9 < RRIz < 1 and offset_avg > ang:
    pass_type = 'RRI-Nadir pass'
else:
    pass_type = 'RRI-Ram pass'

# proximity to Ottawa--------------------------------------------------------
# find where range between the satellite and Ottawa is minimum
min_dist_O=np.min(distance[:-2])
min_index_O = np.where(distance==min_dist_O)

# find where range between the satellite and Ottawa latitudes is minimum
min_index_Olat=misc.find_index(Lat, pLat, 0.1)[0][0]

# set common properties
legend_location = 'upper right'
left, width = .01, .5
bottom, height = .01, .5
width2, height2 = 0.99, 0.99
right = left + width2
top = bottom + height2
date = start_time.strftime('%d-%b-%Y %H:%M:%S')
myFmt = dtformat.DateFormatter("%H:%M:%S")
seconds = dtformat.MicrosecondLocator(interval=10000000)
step = 1
end_ind, ind1 = op.indices_and_intervals(start_time, time_data, step)
LT = start_time-datetime.timedelta(hours=5)
inchesx, inchesy = 19.2, 10.13

OH = 0; H_avg = str("{0:.0f}".format(np.mean(Alt)))
if offset_avg < ang:
    title_slew = 'Slewed pass'
elif 0.9 < RRIz < 1 and offset_avg > ang:
    title_slew = 'RRI-Nadir pass'
else:
    title_slew = 'RRI-Ram pass'
    
if pLon - .25 > Lon[min_index_Olat] :
    sat_pos = 'West'
elif pLon + .25 < Lon[min_index_Olat] :
    sat_pos = 'East'
elif pLon - 0.25 < Lon[min_index_Olat] < pLon + 0.25:
    sat_pos = 'Overhead'
    
if len(F1)==1:
    title_figs = date + ', ' + H_avg + ' km, ' + sat_pos + ' ' +\
    str("{0:.2f}".format(Freq1[0])) + ' MHz, ' + title_slew 
else:
    title_figs = date + ', ' + H_avg + ' km, ' + sat_pos + ' ' +\
    str("{0:.0f}".format(float(Freq1))) + ' MHz and ' + \
    str("{0:.0f}".format(float(Freq2))) + ' MHz, ' + title_slew 
    
# %% plot slew and reception angles 
fig_slew_angles, (ax_reception, ax_slew) = plt.subplots(figsize=(12.06, 10.79),
                                                        nrows = 2, 
                                                        sharex = True)
plt.suptitle(title_figs, fontsize = 18)

offset_angles = offset_los, offset_D1, offset_D2
ylim_min, ylim_max = misc.find_ylimits(offset_los, offset_D1, offset_D2)

offset_leg = [r'$\theta_{SR}$', r'$\phi_{CD1}$',r'$\phi_{CD2}$']
offset_colors = ['black', 'red', '#04d8b2']
offset_title = 'Offset Angles for %s'
figname = 'offset_Angles'
offset_markers = ['.', '.', '.']
offset_ln = ['solid', 'solid', 'solid']
# xlabel, ylabel='Time', 'Observation\n angle (\N{DEGREE SIGN})'
xlabel, ylabel = 'Time', 'Complementary\nReception(\N{DEGREE SIGN})'
ylimits = [ylim_min, ylim_max]

op.line_subplots(ax_reception, time_array_1sec, int_vec,
                  seconds, date, figname, left, bottom, xlabel, ylabel, ylimits,
                  'upper right', 5, 3, *offset_angles, fontsize = 18, 
                  legends=offset_leg,colors=offset_colors, 
                  markers=offset_markers,
                  linestyles=offset_ln)
ax_reception.set_ylim(ylim_min, ylim_max)
misc.arrange_axes(ax_reception, ylim_min, ylim_max)

# Mark closest approach to Ottawa
ax_reception.scatter(time_array_1sec[min_index_O],ylimits[0]+5, marker='*', 
                  c='k', edgecolor='k', s=50, zorder= 50)


ylim_min, ylim_max = 0, 35
op.plot_slew(ax_slew, ylim_min, ylim_max, 5, time_array_1sec, slew_1, str(sa1))
op.plot_slew(ax_slew, ylim_min, ylim_max, 4, time_array_1sec, slew_2, str(sa2))
op.plot_slew(ax_slew, ylim_min, ylim_max, 3, time_array_1sec, slew_3, str(sa3))
op.plot_slew(ax_slew, ylim_min, ylim_max, 2, time_array_1sec, slew_4, str(sa4))
op.plot_slew(ax_slew, ylim_min, ylim_max, 1, time_array_1sec, slew_5, str(sa5))
op.plot_slew(ax_slew, ylim_min, ylim_max, 0, time_array_1sec, slew_6, str(sa6), 
              cb_axis = 'yes', time='yes')
ax_slew.set_ylabel('Closeness to Slew', fontsize = 18)
ax_slew.get_yaxis().set_label_coords(-0.075,0.5)
ax_slew.yaxis.set_tick_params(labelleft=False)
ax_slew.tick_params(axis='both', which='major', labelsize=18)