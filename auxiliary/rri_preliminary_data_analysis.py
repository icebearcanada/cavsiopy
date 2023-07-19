#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for the preliminary analysis of RRI data.

.. toctree::
    :maxdepth: 2   
    calculate_aspect
    import_quaternions
    nan_counter
    plot_invalid_points
    quick_spectrogram
    voltage_reader
    
@author: ceren
"""

import datetime
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import matplotlib.dates as dtformat

from geopy.distance import geodesic
import pyIGRF

import bisect
from spacepy import pycdf

import cavsiopy.miscellaneous as misc
import cavsiopy.use_rotation_matrices as RM

def voltage_reader(filename):
    """
    function to read induced complex voltages at the RRI

    Parameters
    ----------
    filename : str
        RRI file name including the file path.

    Returns
    -------
    v1 : numpy.ndarray[float]
        voltage1 (mV).
    v2 : numpy.ndarray[float]
        voltage2 (mV).
    v3 : numpy.ndarray[float]
        voltage3 (mV).
    v4 : numpy.ndarray[float]
        voltage4 (mV).

    """
    file = h5py.File(filename,'r')
    # %% import voltages from RRI data set
    v1 = file.get('RRI Data/Radio Data Monopole 1 (mV)')[:,:]
    v2 = file.get('RRI Data/Radio Data Monopole 2 (mV)')[:,:]
    v3 = file.get('RRI Data/Radio Data Monopole 3 (mV)')[:,:]
    v4 = file.get('RRI Data/Radio Data Monopole 4 (mV)')[:,:]
    # %% Calculate the orientation and ellipticity angles from Stokes parameters
    v1 = v1.flatten(); v2 = v2.flatten()
    v3 = v3.flatten(); v4 = v4.flatten()

    return v1, v2, v3, v4

def timeTicks(start_time, x, pos):
    """
    User-defined function for formatting

    Parameters
    ----------
    start_time : datetime
        start time of the experiment.
    x : TYPE
        tick value.
    pos : TYPE
        position.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    d = start_time + datetime.timedelta(seconds=x)
    return d.strftime("%H:%M:%S")

def quick_spectrogram(filename, start_time, \
                          fs = 62500.333, freq =10422000):
    """
    function to plot spectrograms for one frequency experiments of the RRI.
    
    Parameters
    ----------
    filename : str
        filename including the path.
    start_time : datetime
        start of the experiment.
    fs : float, optional
        sampling rate. The default is 62500.333.
    freq : float, optional
        frequency of the experiment in Hz. The default is 10422000.

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module. plot of spectrograms.
    ax : axes._subplots.AxesSubplot
        Axes object.

    """

    f1 = freq/1e6

    v1, v2, v3, v4 = voltage_reader(filename)
    
    vx = np.nan_to_num(v1+1j*v2)
    vy = np.nan_to_num(v3+1j*v4)
    
    date_str = start_time.strftime('%Y-%m-%d')
       
    nFFT = 5208 # 5 msec window
    noverlap = 62 # 15 msec window

    nrows = 2; ncols = 1
    fig, (ax) = plt.subplots(nrows=nrows, ncols = ncols, sharex = True, 
                             sharey = True)
    plt.suptitle(date_str, fontsize = 14)
    spec1, freqx1, ts1, im1 = ax[0].specgram(vx, Fs=fs, 
                                            Fc = freq,
                                            NFFT=nFFT, 
                                          sides='twosided', 
                                          cmap=plt.cm.turbo,
                                           vmax = -60, vmin=-72.5 
                                          )

    ax[0].text(.95, .95, 'ax',  color='w' ,  fontsize = 'large', 
               transform=ax[0].transAxes)    

    ax[0].text(0.03, .95, f'Central Freq = {f1} MHz',  
                 color='w' , fontsize = 'medium', 
               transform=ax[0].transAxes)  
    
    spec2, freqy1, ts2, im2 = ax[1].specgram(vy, Fs=fs, 
                                            Fc = freq,
                                            NFFT=nFFT,
                                           sides='twosided', cmap=plt.cm.turbo,
                                            vmax = -60, vmin=-72.5 
                                           )
    ax[1].text(.95, .95, 'ay',  color='w' ,  fontsize = 'large', 
               transform=ax[1].transAxes) 
    cb2 = fig.colorbar(im2, ax=ax[:])
    cb2.ax.tick_params(labelsize=18) 
      
    # Prettify
    scale = 1e6                     # KHz
    ticks = matplotlib.ticker.FuncFormatter(start_time, lambda x, \
                                            pos: '{0:g}'.format(x/scale))
    
    formatter = matplotlib.ticker.FuncFormatter(timeTicks)
    ax[1].xaxis.set_major_formatter(formatter)             
    ax[1].set_xlabel('Time (UT)', fontsize = 18)
    for j in range(2):
        ax[j].tick_params(axis='both', which='major', labelsize=18)
        ax[j].yaxis.set_major_formatter(ticks)
        ax[j].set_ylabel('Frequency (MHz)', fontsize = 18)
        
    return fig, ax

def nan_counter(data, bin_size = 62500):
    """
    Function to count zeroes or NaNs in RRI voltage arrays.

    Parameters
    ----------
    data : numpy.ndarray[float]
        voltage array.
    bin_size : float, optional
        sampling rate in Hz. The default is 62500.

    Returns
    -------
    percentage_zero : numpy.ndarray[float]
        percentage of data points which are zero (%).
    percentage_nan : numpy.ndarray[float]
        percentage of data points which are Not a Number-NaN (%).
    percent_valid : numpy.ndarray[float]
        percentage of data points which are not zero or NaN (%).
    percent_invalid : numpy.ndarray[float]
        percentage of data points which are zero or NaN (%).

    """


    number_of_data = len(data)
    intervals = int(number_of_data/bin_size)
    remnant = number_of_data%bin_size
    end = number_of_data - remnant
    non_zero = np.empty(intervals);
    number_of_nan = np.empty(intervals);
    number_of_zero = np.empty(intervals);
    total_invalid = np.empty(intervals);
    for i in range(0,end,bin_size):
        j = int(i/bin_size)
        V = data[i:i+bin_size]
        non_zero[j] = np.count_nonzero(V)
        number_of_zero[j] = bin_size-non_zero[j]
        number_of_nan[j] = np.count_nonzero(np.isnan(V))
        total_invalid[j] = number_of_nan[j] + number_of_zero[j]

    # percentage of valid points
    percent_valid = (1-np.true_divide(total_invalid, bin_size))*100
    percent_invalid = (np.true_divide(total_invalid, bin_size))*100
    percentage_zero= (1-np.true_divide(non_zero, bin_size))*100
    percentage_nan = (np.true_divide(number_of_nan, bin_size))*100

    return percentage_zero, percentage_nan, percent_valid, percent_invalid

def plot_invalid_points(ax, time_data, percent_invalid, bin_size_low = 62500):
    """
    Function to insert a bar at the bottom of a figure to show 
    the number of invalid points.

    Parameters
    ----------
    ax : axes
        Axes object of matplotlib.
    time_data : numpy.ndarray[float]
        time data from RRI data file.
    percent_invalid : numpy.ndarray[float]
        percentage of data points which are zero or NaN (%).
    bin_size_low : float, optional
        sampling rate in Hz. The default is 62500.

    Returns
    -------
    None.

    """
    
    multiple_sec=1
    tmin = datetime.datetime(1968, 5, 24)
    tt = np.array([tmin +
                   datetime.timedelta(seconds=time_data[i*multiple_sec])
                   for i in range(0, len(percent_invalid))])
    
    # specify colors
    colors_for_mapping = ['green', 'green', 'green', 'green', 'green',
                          'green', 'green', 'green','red','red','red']
    #  create color map using the specified colors
    cmap1 = LinearSegmentedColormap.from_list("mycmap", colors_for_mapping)

    # set discrete levels for the color map
    norm = BoundaryNorm([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], cmap1.N)

    # define ticks
    ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    bcolors = misc.determine_background_colors(percent_invalid, ticks,
                                               colors_for_mapping)

    for i in range(0, len(percent_invalid)-1):
        width = tt[i+1]-tt[i]
        rect1 = patches.Rectangle((tt[i], 0), width, 100,
                                 edgecolor='None',
                                 facecolor=str(bcolors[i]),
                                 fill=True, lw=0)
        c1 = ax.add_patch(rect1)

    dummy_plot = ax.scatter(tt, percent_invalid,
                               c=percent_invalid, cmap=cmap1,
                               norm=norm, s=0)
    xax = ax.axes.get_xaxis()
    xax = xax.set_visible(False)

    yax = ax.axes.get_yaxis()
    yax = yax.set_visible(False)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.patch.set_alpha(0.001)

    return

def calculate_aspect(time_array, pLat, pLon, pAlt, lat, lon, alt, RRI_ned):
    """
    Function to calculate the magnetic aspect angle.
    Depends on pyIGRF.

    Parameters
    ----------
    time_array : datetime.datetime
        experiment time interval as datetime array.
    pLat : TYPE
        geodetic latitude of the target.
    pLon : TYPE
        DESCRIPTION.
    pAlt : TYPE
        DESCRIPTION.
    lat : numpy.ndarray[float]
        geodetic latitude of the target.
    lon : numpy.ndarray[float]
        geodetic longitude of the target.
    alt : numpy.ndarray[float]
        altitude of the target (km).
    RRI_ned : numpy.ndarray[float]
        RRI pointing vector in NOrth-East-Down (NED) coordinate system).

    Returns
    -------
    aspect_los : numpy.ndarray[float]
        aspect angle.

    """

    start_time = time_array[0]
    DOY=start_time.timetuple().tm_yday

    # Line-of-sight distance between the point and satellite
    distance=np.empty(alt.shape); distance.fill(np.nan)
    # Ray vector
    los_ned=[]; los_enu=[];

    # Magnetic field components from IGRF
    D = np.empty(alt.shape); I = np.empty(alt.shape);
    H = np.empty(alt.shape);
    BX = np.empty(alt.shape); BY = np.empty(alt.shape);
    BZ = np.empty(alt.shape); BF = np.empty(alt.shape);

    # Magnetic field vectors for coordinate transformation
    MF_ned=[]; MF_enu=[];

    # Unit vectors for magnetic field components
    dBx= np.empty(alt.shape); dBy= np.empty(alt.shape);
    dBz= np.empty(alt.shape);

    # angle between magnetic field and LOS ray path
    aspect_los=np.empty(alt.shape); aspect_los.fill(np.nan)

    # =============================================================================
    #  Calculate day of the year for IGRF function.
    # =============================================================================
    date_year=start_time.year
    if (date_year%4==0):
        doy_igrf=start_time.year+DOY/366
    else:
        doy_igrf=start_time.year+DOY/365.25

    for i in range(0, len(alt)):
    # =============================================================================
    # calculate the ray vector
    # =============================================================================
        # positive down vector (D in NED)
        dz_down= pAlt-alt[i]
        # positive up vector (E in ENU)
        dz_up=alt[i]-pAlt

        # find the changes in latitude by keeping longitude constant
        point1 = (pLat, pLon)
        ePOP = (lat[i], pLon)
        # Distance due to changes in latitude (y axis in this case)
        dy = geodesic(point1, ePOP).kilometers

        # find the changes in latitude by keeping spacecraft and ground point
        # latitude the same
        point2 = (lat[i], pLon)
        ePOP = (lat[i], lon[i])
        # Distance due to changes in longitude (x axis in this case)
        dx = geodesic(point2, ePOP).kilometers
        import math
        distance[i] = math.sqrt(dx**2 + dy**2 + dz_down**2)

        # Y is positive northward, X is positive eastward. In other cases change
        # of sign is needed:

        # if spacecraft is in the south of Ottawa
        if lat[i] < pLat:
            dy = -dy
        # if spacecraft is in the west of Ottawa
        if 360+lon[i] < 360+pLon:
            dx = -dx

        # NED: down is positive
        los_ray_ned = [dy/distance[i], dx/distance[i], dz_down/distance[i]]
        # ENU: up is positive
        los_ray_enu = [dx/distance[i], dy/distance[i], dz_up/distance[i]]
        # append the ray vector values
        los_ned.append(los_ray_ned)
        los_enu.append(los_ray_enu)

    # =============================================================================
    # # calculate the aspect angle
    # =============================================================================
    # D: magnetic declination, I: magnetic inclination (dip), H: horizontal comp.
        [D[i], I[i], H[i], BX[i], BY[i], BZ[i], BF[i]] = pyIGRF.igrf_value(
        lat[i], lon[i], alt[i], doy_igrf)

        # Unit vectors along the alt, lon, lat directions
        dBx[i] = BX[i]/BF[i] # north component
        dBy[i] = BY[i]/BF[i] # east component
        dBz[i] = BZ[i]/BF[i] # down component (IGRF gives in NED coordinates)

        # magnetic field unit vectors
        B_ned= [dBx[i], dBy[i], dBz[i]] # NED
        MF_ned.append(B_ned)

        MF_enu_init = RM.NED2ENU(B_ned)
        MF_enu.append(MF_enu_init)

    # =============================================================================
    # Find the angle between the magnetic field vector and los ray vector
        dot_product_los = np.dot(los_ray_ned, B_ned)
        aspect_los[i] = math.acos(dot_product_los)*180.0/np.pi

    return aspect_los

def import_QUA_file(file_QUA, start_date, end_date):
    """
    imports quaternions and accuracy data from the quaternion file:
    (https://epop.phys.ucalgary.ca/data-handbook/attitude-quaternion-file-cas_attquat/)
    
    Extracts data for the experiment time interval specified by start_date and 
    end_date.
    
    Parameters
    ----------
    file_QUA : str
        QUA filename with file path
    start_date : datetime.datetime
        beginning of the data sampling interval
    end_date : datetime.datetime
        end of the data sampling interval

    Returns
    -------
    dict : dict
        time: datetime.datetime
            experiment time interval
        q: numpy.array[float]
            quaternions: qR, qI, qJ, qK
        roll : np.array[float]
            roll angle; rotation around x axis (degrees)
        pitch : np.array[float]
            pitch angle; rotation around y axis (degrees)
        yaw : np.array[float]
            yaw angle; rotation around z axis (degrees)
        ds : int
            data sources:
            0:	Dropout (added to indicate the start/
                        end of a ‘NaN’ filled period)
            1:	Uncalibrated coarse sun sensor and magnetometer solution 
                (very large uncertainty, > 10 degrees)
            2:	Coarse Sun Sensors and Magnetometers 
                (large uncertainty of ~8 degrees)
            3:	Onboard-fused star sensor (uncertainty of ~0.1 degrees)
            4:	Star sensor A solution (SSA) (uncertainty of ~0.06 degrees)
            5:	Star sensor B solution (SSB) (uncertainty of ~0.06 degrees)
            6:	Star sensor A+B (locally fused solution, 
                                uncertainty of ~0.01 degrees)
    """

    cdf = pycdf.CDF(file_QUA)
    start_ind = bisect.bisect_left(cdf['Timestamp'], start_date)
    stop_ind = bisect.bisect_left(cdf['Timestamp'], end_date) + 1
    accuracy = cdf['Data_Source_Next'][start_ind:stop_ind]
    time = cdf['Timestamp'][start_ind:stop_ind]
    q = cdf['q'][start_ind:stop_ind, :]
    roll = cdf['Roll'][start_ind:stop_ind]
    pitch = cdf['Pitch'][start_ind:stop_ind]
    yaw = cdf['Yaw'][start_ind:stop_ind]
    cdf.close()
     
    return {'time': time, 'q':q, 
            'roll': roll, 'pitch': pitch, 'yaw': yaw, 'accuracy': accuracy 
            }
