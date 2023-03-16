#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 14:46:24 2021

epmeris_importer includes functions to 
1. import RRI ephemeris from hdf5 data file (import_rri_ephemeris)
2. import ephemeris only related to GEI frame from hdf5 data file (import_RRI_ephemeris_GEI)
3. import ephemeris accuracy from RRI hdf5 data file
4. read CelesTrak TLE file parameters (import_tle_parameters)
5. compare Celestrak TLE parameters with the calculated TLE parameters (compare_orbital)

@author: ceren
"""
import numpy as np
import h5py
import datetime

def import_rri_ephemeris(file_RRI):
# Import RRI data of selected file
    file = h5py.File(file_RRI,'r')

# Attitude information from RRI
    roll = np.asarray(file["CASSIOPE Ephemeris/Roll (deg)"])  # X
    pitch = np.asarray(file["CASSIOPE Ephemeris/Pitch (deg)"])  # Y
    yaw = np.asarray(file["CASSIOPE Ephemeris/Yaw (deg)"])  # Z
    GEIx, GEIy, GEIz = np.asarray(file["CASSIOPE Ephemeris/GEI Position (km)"])
    GEIVx, GEIVy, GEIVz =file.get('CASSIOPE Ephemeris/GEI Velocity (km per s)')[:]
    GEOx, GEOy, GEOz =file.get('CASSIOPE Ephemeris/GEO Position (km)')[:]
    GEOVx, GEOVy, GEOVz =file.get('CASSIOPE Ephemeris/GEO Velocity (km per s)')[:]

# Location information from RRI
    Alt = file.get('CASSIOPE Ephemeris/Altitude (km)')[:]
    Lat = file.get('CASSIOPE Ephemeris/Geographic Latitude (deg)')[:]
    Lon = file.get('CASSIOPE Ephemeris/Geographic Longitude (deg)')[:]
    MLat= file.get('CASSIOPE Ephemeris/Magnetic Latitude (deg)')[:]
    MLon= file.get('CASSIOPE Ephemeris/Magnetic Longitude (deg)')[:]
    MLT = file.get('CASSIOPE Ephemeris/MLT (hr)')[:]

#  Time information from RRI
    time_data = \
        file.get('CASSIOPE Ephemeris/Ephemeris MET (seconds since May 24, 1968)')[:]
    start_time = datetime.datetime(1968,5,24) + \
        datetime.timedelta(seconds=time_data[0])
    end_time = datetime.datetime(1968,5,24) + \
        datetime.timedelta(seconds=time_data[-1])
    time_array_1sec = np.array([start_time + \
                                datetime.timedelta(seconds = i*1)\
                                for i in range(0,len(time_data))])
    datetime_time_array_1sec = [datetime.datetime(1968,5,24)+ \
                                   datetime.timedelta(seconds=time_data[i])\
                                       for i in range(0,len(Alt))]

    return time_data, datetime_time_array_1sec, time_array_1sec, start_time, \
        end_time, Lat, Lon, Alt, MLat, MLon, MLT, \
            GEOx, GEOy, GEOz, GEOVx, GEOVy, GEOVz,\
            GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz, roll, pitch, yaw

def import_rri_ephemeris_gei(file_RRI):
# import only GEI position and velocity
    # Import RRI data of selected file
    file = h5py.File(file_RRI,'r')

    # Attitude information from RRI
    roll = np.asarray(file["CASSIOPE Ephemeris/Roll (deg)"])  # X
    pitch = np.asarray(file["CASSIOPE Ephemeris/Pitch (deg)"])  # Y
    yaw = np.asarray(file["CASSIOPE Ephemeris/Yaw (deg)"])  # Z
    GEIx, GEIy, GEIz = np.asarray(file["CASSIOPE Ephemeris/GEI Position (km)"])
    GEIVx, GEIVy, GEIVz =file.get('CASSIOPE Ephemeris/GEI Velocity (km per s)')[:]


    return GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz, roll, pitch, yaw
# =============================================================================
# %% RRI accuracy
# =============================================================================
def import_accuracy_from_RRI(file_RRI):
    # %%  Import RRI data of selected file
    file = h5py.File(file_RRI,'r')
    time_data = \
        file.get('CASSIOPE Ephemeris/Ephemeris MET (seconds since May 24, 1968)')[:]
    start_time = datetime.datetime(1968,5,24) + \
        datetime.timedelta(seconds=time_data[0])

    time_array = np.array([start_time + \
                                datetime.timedelta(seconds = i*1)\
                                for i in range(0,len(time_data))])
    date = start_time.strftime('%d-%b-%Y')

    accuracy = np.asarray(file["CASSIOPE Ephemeris/Attitude Accuracy (0=Dropout,1=Rough,2=Coarse,3=Moderate,4=Fine,9=No data)"])

    return time_array, date, accuracy

def import_tle_parameters(filename, filedate, DOY):
    """
    Created on Tue Jul 13 22:57:40 2021

    Reads Celestrak TLE file obtained from: https://celestrak.com/

    Parameters
    ----------
    filename (str): filename for TLE file
    filedate (str): date of RRI passage
    DOY (int): day of year of RRI passage

    Returns
    -------
    tle_epoch (float): epoch of TLE file
    tle_inc (float): satellite inclination from TLE file
    tle_ap (float): satellite argument of perigee from TLE file
    tle_raan(float): right ascension of ascending node from TLE file
    tle_ecc (float): eccentricity from TLE file
    tle_ma (float): mean anomaly

    @author: ceren
    """
    fid = open(filename, 'r')
    epoch_year=filedate[2:4] # read the last two digits of the year
    sat_date=epoch_year+str(DOY)+'.'

    # initialize lines
    line1=[]
    line2=[]
    i =  0; r = 0; a = 0; j = 0;

    for row in fid:
        if row.startswith('1 '):
            l1 = row
            i += 1

            # line 1 for the specific date is when l1[18:24]==sat_date
            if (l1[18:24]==sat_date):
                r= i
                a += 1
                line1.append(l1)

        elif row.startswith('2 '):
            l2 = row
            j += 1
            # line 2 for the specific date is when the index is equal to
            # sat_date index
            if (j == r):
                line2.append(l2)
        else:
            a=0

    fid.close()

    if a>0:
        tle_epoch = float(line1[0][18:32])
        tle_inc = float(line2[0][8:16])
        tle_raan = float(line2[0][17:25])
        tle_ecc = float(line2[0][26:33]) * (10 ** -7)
        tle_ap = float(line2[0][34:42])
        tle_ma = float(line2[0][43:51])
    else:
        print('no values to check for: missing TLE')
        tle_epoch = np.empty(0)
        tle_inc = np.empty(0)
        tle_raan = np.empty(0)
        tle_ecc = np.empty(0)
        tle_ap = np.empty(0)
        tle_ma = np.empty(0)

    return tle_epoch, tle_inc, tle_ap, tle_raan, tle_ecc, tle_ma

# =============================================================================
# # comparison between the TLE file and calculated TLE values
# =============================================================================
def compare_orbital(sat_inc, sat_raan, sat_ap, tle_inc, tle_ap, tle_raan, date):

    if np.size(tle_inc)>0:
        print(" ".join(['for the satellite epoch', str(date), ':\n']))
        sap_comp = np.empty(3)

        sinc = sat_inc*180/np.pi
        sap_comp[0] = tle_inc - sinc
        if np.isclose(tle_inc, sinc , atol=1)==True:
            print(" ".join(['TLE_inc =',str(tle_inc),', sat_inc =', str(sinc),
                            ': calculations are correct\n']))
        else:
            print('problem with satellite inclination: \n',
                  'check your satellite orbital elements calculations,',
                  'the tle or the ephemeris file\n')

        sarn = sat_raan*180/np.pi
        sap_comp[1] = tle_raan - sarn
        if np.isclose(tle_raan, sarn, atol=1)==True:
            print(" ".join(['TLE_raan =',str(tle_raan),', sat_raan =',
                            str(sarn), ': calculations are correct\n']))
        else:
            print(" ".join(['TLE_raan =',str(tle_raan),', sat_raan =',
                            str(sarn), 'problem with satellite raan:',
                            '\n check your satellite orbital elements calculations,',
                            'the tle or the ephemeris file']))

        sat_ap_deg = sat_ap*180/np.pi
        sap_comp[2] = tle_ap - sat_ap_deg
        if np.isclose(tle_ap, sat_ap_deg, atol=1)==True:
            print(" ".join(['TLE_ap =',str(tle_ap),', sat_ap =',
                            str(sat_ap_deg), ': calculations are correct\n']))
        else:
            print(" ".join(['TLE_ap =',str(tle_ap),', sat_ap =',
                            str(sat_ap_deg),'problem with satellite ap:',
                            '\n check your satellite orbital elements calculations,',
                            'the tle or the ephemeris file\n']))
    else:
        print('No TLE to compare with the calculated satellite orbital elements')

    return sap_comp
