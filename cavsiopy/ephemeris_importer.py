#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ephemeris_importer module imports several types of data files to obtain 
spacecraft ephemeris. Values are returned as Python dictionaries.

Swarm-E specific data files are: CAS_ephemeris, Cas_AttQUAT, and RRI

Files standard for all spacecraft: sp3, TLE 

ephemeris_importer also includes functions to compute satellite orbital 
elements using GEIJ2K position and velocity (calculate_orbital_elements), 
and to compare computed values with data in TLE file (compare_orbital).

.. toctree::
  :maxdepth: 2 
  cas_ephemeris
  rri_ephemeris
  sp3_ephemeris
  import_tle
  calculate_orbital_elements
  compare_orbital

@author: ceren
"""
import datetime
import numpy as np
import h5py

def cas_ephemeris(file_cas, time_start, time_end):
    '''
    Imports CAS_ephemeris files of Swarm-E:
    Ref:(https://epop.phys.ucalgary.ca/data-handbook/cas-ephemeris-text-files/)  
    
    Extracts satellite ephemeris info for the experiment time interval
    specified by time_start and time_end.
    
    Parameters
    ----------
    file_cas : str
        file name including path.
    time_start : datetime.datetime
        start of the experiment.
    time_end : datetime.datetime
        end of the experiment.

    Returns
    -------
    dict : dict
        keys and properties listed below
        srow : int
            row of the start time in the cas_ephemeris file.
        erow : int
            row of the end time in the cas_ephemeris file.
        start : datetime.datetime
            start time in CAS_ephemeris file.
        end : datetime.datetime
            end time in CAS_ephemeris file.
        time_array : datetime.datetime
            time of the observations.
        Lat : numpy.ndarray[float]
            Geodetic latitude.
        Lon : numpy.ndarray[float]
            Geodetic longitude.
        Alt : numpy.ndarray[float]
            geodetic altitude (km)  .      
        MLat: numpy.ndarray[float]
            Magnetic latitude (deg).
        MLon: numpy.ndarray[float]
            Magnetic longitude (deg).
        MLT: numpy.ndarray[float]
            Magnetic local time (h).
        GEIx : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (X-km).
        GEIy : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (Y-km).
        GEIz : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (Z-km).
        GEIVx : numpy.ndarray[float]
            X component of the spacecraft velocity in GEI coordinates (Vx-km/s).
        GEIVy : numpy.ndarray[float]
            Y component of the spacecraft velocity in GEI coordinates (Vy-km/s).
        GEIVz : numpy.ndarray[float]
            Z component of the spacecraft velocity in GEI coordinates (Vz-km/s).
        GSMx : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (X-km).
        GSMy : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (Y-km).
        GSMz : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (Z-km).
        roll : numpy.ndarray[float]
            roll angle; rotation around x axis (degrees).
        pitch : numpy.ndarray[float]
            pitch angle; rotation around y axis (degrees).
        yaw : numpy.ndarray[float]
            yaw angle; rotation around z axis (degrees).
        acc : int
            accuracy of the attitude solution.
            (0 = Dropout, 1 = Rough, 2 = Coarse, 3 = Moderate,  4 = Fine)
            
    Examples
    --------
    dict_cas = ei.cas_ephemeris(file_CAS, time_start, time_end)
    
    GEIx = dict_cas['GEIx']  
    '''
    
    cas = np.loadtxt(fname=file_cas, skiprows = 2, usecols= range(0,21)) 
    t0 = time_start.time().strftime('%H%M%S')
    t1 = time_end.time().strftime('%H%M%S')
    start_time = int(t0)
    end_time = int(t1)
       
    # find start row in the file
    srow = int(np.where(cas == start_time)[0])
    # find end row in the file
    erow = int(np.where(cas == end_time)[0]) + 1

    GEIx = cas[srow:erow,2]
    GEIy = cas[srow:erow,3]
    GEIz = cas[srow:erow,4]

    GEIVx = cas[srow:erow,5]
    GEIVy = cas[srow:erow,6]
    GEIVz = cas[srow:erow,7]
    
    GSMx = cas[srow:erow,8]
    GSMy = cas[srow:erow,9]
    GSMz = cas[srow:erow,10]
    
    Lat = cas[srow:erow,11]    
    Lon = cas[srow:erow,12]
    Alt = cas[srow:erow,13]
    
    MLat = cas[srow:erow,14] 
    MLon = cas[srow:erow,15] 
    MLT = cas[srow:erow,16] 
    
    yaw = cas[srow:erow,17]    
    pitch = cas[srow:erow,18]
    roll = cas[srow:erow,19]
    acc = cas[srow:erow,20] # attitude accuracy
    
    time_array = np.array([time_start + datetime.timedelta(seconds = i*1) \
                           for i in range(0,len(Alt))])
    
    return { 'srow':srow, 'erow':erow, 'time_array': time_array, 
            'Lat': Lat, 'Lon': Lon, 'Alt': Alt,
            'MLat': MLat, 'Mlon': MLon, 'MLT': MLT, 
            'GEIx': GEIx, 'GEIy': GEIy, 'GEIz': GEIz, 
            'GEIVx': GEIVx, 'GEIVy': GEIVy, 'GEIVz': GEIVz, 
            'GSMx': GSMx, 'GSMy': GSMy, 'GSMz': GSMz, 
            'roll': roll, 'pitch': pitch, 'yaw': yaw, 'accuracy': acc
        } 

# =============================================================================
# RRI
# =============================================================================
def rri_ephemeris(file_rri):
    
    """
    Imports rri_ephemeris and returns the ephemeris information within a 
    dictionary.
       
    Dictionary keys are listed below:
    
    
    Parameters
    ----------
    file_rri : TYPE
        RRI file name including path

    Returns
    -------
    dict:
        keys and properties of the keys are below:
        time_data : numpy.ndarray[float]
            seconds since May 24, 1968.
        time_array : datetime.datetime
            time of the observations.
        start_time : datetime.datetime
            start_time of the experiment.
        end_time : datetime.datetime
            end time of the experiment.
        Lat : numpy.ndarray[float]
            Geodetic latitude (degrees).
        Lon : numpy.ndarray[float]
            Geodetic longitude (degrees).
        Alt : numpy.ndarray[float]
            geodetic altitude (km).        
        MLat: numpy.ndarray[float]
            Magnetic latitude (deg).
        MLon: numpy.ndarray[float]
            Magnetic longitude (degrees).
        MLT: numpy.ndarray[float]
            Magnetic local time (h).
        GEOx : numpy.ndarray[float]
            Spacecraft position in GEO coordinates (X-km).
        GEOy : numpy.ndarray[float]
            Spacecraft position in GEO coordinates (Y-km).
        GEOz : numpy.ndarray[float]
            Spacecraft position in GEO coordinates (Z-km).
        GEOVx : numpy.ndarray[float]
            X component of the spacecraft velocity in GEO coordinates (Vx-km/s).
        GEOVy : numpy.ndarray[float]
            Y component of the spacecraft velocity in GEO coordinates (Vy-km/s).
        GEOVz : numpy.ndarray[float]
            Z component of the spacecraft velocity in GEO coordinates (Vz-km/s).
        GEIx : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (X-km).
        GEIy : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (Y-km).
        GEIz : numpy.ndarray[float]
            Spacecraft position in GEI coordinates (Z-km).
        GEIVx : numpy.ndarray[float]
            X component of the spacecraft velocity in GEI coordinates (Vx-km/s).
        GEIVy : numpy.ndarray[float]
            Y component of the spacecraft velocity in GEI coordinates (Vy-km/s).
        GEIVz : numpy.ndarray[float]
            Z component of the spacecraft velocity in GEI coordinates (Vz-km/s).
        roll : numpy.ndarray[float]
            roll angle; rotation around x axis (degrees).
        pitch : numpy.ndarray[float]
            pitch angle; rotation around y axis (degrees).
        yaw : numpy.ndarray[float]
            yaw angle; rotation around z axis (degrees).
        accuracy : int
            accuracy of the attitude solution.
            (0 = Dropout, 1 = Rough, 2 = Coarse, 3 = Moderate,  4 = Fine)
            
        
    Examples
    --------
    dict_rri = ei.rri_ephemeris(file_RRI)
    
    Lat = dict_rri['Lat']    

    """
    
    # Import RRI data of selected file
    file = h5py.File(file_rri,'r')

    # Attitude information from RRI
    roll = np.asarray(file['CASSIOPE Ephemeris/Roll (deg)'])  # X
    pitch = np.asarray(file['CASSIOPE Ephemeris/Pitch (deg)'])  # Y
    yaw = np.asarray(file['CASSIOPE Ephemeris/Yaw (deg)'])  # Z
    GEIx, GEIy, GEIz = np.asarray(file['CASSIOPE Ephemeris/' + \
                                       'GEI Position (km)'])
    GEIVx, GEIVy, GEIVz =file.get('CASSIOPE Ephemeris/' + \
                                       'GEI Velocity (km per s)')[:]
    GEOx, GEOy, GEOz =file.get('CASSIOPE Ephemeris/' + \
                                       'GEO Position (km)')[:]   
    GEOVx, GEOVy, GEOVz =file.get('CASSIOPE Ephemeris/' + \
                                       'GEO Velocity (km per s)')[:]

    # Location information from RRI
    Alt = file.get('CASSIOPE Ephemeris/Altitude (km)')[:]
    Lat = file.get('CASSIOPE Ephemeris/Geographic Latitude (deg)')[:]
    Lon = file.get('CASSIOPE Ephemeris/Geographic Longitude (deg)')[:]
    MLat= file.get('CASSIOPE Ephemeris/Magnetic Latitude (deg)')[:]
    MLon= file.get('CASSIOPE Ephemeris/Magnetic Longitude (deg)')[:]
    MLT = file.get('CASSIOPE Ephemeris/MLT (hr)')[:] 
    
    # accuracy
    accuracy = np.asarray(file['CASSIOPE Ephemeris/Attitude Accuracy' + \
                               ' (0=Dropout,1=Rough,2=Coarse,' + \
                                   '3=Moderate,4=Fine,9=No data)'])  

    #  Time information from RRI    
    time_data = \
        file.get('CASSIOPE Ephemeris/Ephemeris MET ' + \
                 '(seconds since May 24, 1968)')[:]
    
    start_time = datetime.datetime(1968,5,24) + \
        datetime.timedelta(seconds=time_data[0])
    
    end_time = datetime.datetime(1968,5,24) + \
        datetime.timedelta(seconds=time_data[-1])
    
    time_array = np.array([start_time + datetime.timedelta(seconds = i*1)\
                                for i in range(0,len(time_data))])
    return {
        'time_data':time_data, 'time_array': time_array, 
        'start_time': start_time, 'end_time': end_time,
        'Lat': Lat, 'Lon': Lon, 'Alt': Alt,
        'MLat': MLat, 'MLon': MLon, 'MLT': MLT, 
        'GEOx':GEOx, 'GEOy': GEOy, 'GEOz': GEOz, 
        'GEOVx': GEOVx, 'GEOVy': GEOVy, 'GEOVz': GEOVz,
        'GEIx': GEIx, 'GEIy': GEIy, 'GEIz': GEIz, 
        'GEIVx': GEIVx, 'GEIVy': GEIVy, 'GEIVz': GEIVz, 
        'roll': roll, 'pitch': pitch, 'yaw': yaw, 'accuracy': accuracy
    } 
                      
            
def sp3_ephemeris(file_SP3, start_date, end_date):
    """
    imports position and velocity data in the International Terrestrial 
    Reference Frame (ITRF) from the sp3 data file:
    
    Ref: (https://epop.phys.ucalgary.ca/data-handbook/orbit-geo-sp3-file/)
    
    Extracts satellite ephemeris info for the experiment time interval
    specified by time_start and time_end.

    Parameters
    ----------
    file_SP3 : str
        SP3 filename with file path
    start_date : datetime.datetime
        beginning of the data sampling interval
    end_date : datetime.datetime
        end of the data sampling interval

    Returns
    -------
    dict : dict 
        keys and properties listed below:
        srow : int
            row of the experiment start time in sp3 file
        erow: int
            row of the experiment end time in sp3 file
        time_array_gps : datetime.datetime
            time array for the whole day: note that sp3 files for Swarm-E may 
            have some offset as they use GPS time rather than the UT time.
            to check the beginning and end times of the data file simply do:
                dict_sp3['time_array'][0]
                dict_sp3['time_array'][-1]
        time_array_ut: datetime.datetime
            time array in ut. corrected from gps time.
        time_experiment: datetime.datetume
            time array for the experiment interval in ut.
        ITRFx : numpy.ndarray[float]
            Spacecraft position in ITRF coordinates (X-km).
        ITRFy : numpy.ndarray[float]
            Spacecraft position in ITRF coordinates (Y-km).
        ITRFz : numpy.ndarray[float]
            Spacecraft position in ITRF coordinates (Z-km).
        ITRFVx : numpy.ndarray[float]
            X component of the spacecraft velocity in ITRF coordinates(Vx-dm/s).
        ITRFVy : numpy.ndarray[float]
            Y component of the spacecraft velocity in ITRF coordinates(Vy-dm/s).
        ITRFVz : numpy.ndarray[float]
            Z component of the spacecraft velocity in ITRF coordinates(Vz-dm/s).

    Examples
    --------
    dict_sp3 = ei.sp3_ephemeris(file_SP3, start_date, end_date)
    
    ITRFx = dict_sp3['ITRFx']    

    """
        
    fid = open(file_SP3, 'r')   
    
    # initialize lines: TL: time line, VL: velocity line, PL = position line
    time_sat = []
    vel_sat = []
    pos_sat = []
    
    k =  0; a = 0; j = 0;
    
    for row in fid:
        # time line
        if row.startswith('* '):
            TL = row
            time_sat.append(TL) 
            k += 1
        
        # position line
        elif row.startswith('PL63 '):
            PL = row
            pos_sat.append(PL) 
            j += 1
        
        # velocity line
        elif row.startswith('VL63 '):
            VL = row
            vel_sat.append(VL) 
        else:
            a=0
    fid.close()  
    
    year_init= int(time_sat[0][3:7])
    month_init= int(time_sat[0][8:10])
    day_init=  int(time_sat[0][11:13])
    hour_init = int(time_sat[0][14:16])
    min_init = int(time_sat[0][17:19])
    sec_init= int(float(time_sat[0][20:31]))
    # find when the sp3 file time starts
    start_sp3 = datetime.datetime(year_init, month_init, day_init, 
                             hour_init, min_init, sec_init)  

    # initialize varibles
    Vx = np.empty(k); Vy = np.empty(k); Vz = np.empty(k)
    x = np.empty(k); y = np.empty(k); z = np.empty(k)
    year = np.empty(k); month = np.empty(k); day = np.empty(k)
    hour = np.empty(k); minute = np.empty(k); sec = np.empty(k)
    
    for i in range (0,k):
        year[i] = int(time_sat[i][3:7])
        month[i] = int(time_sat[i][8:10])
        day[i] =  int(time_sat[i][11:13])
        hour[i] = int(time_sat[i][14:16])
        minute[i] = int(time_sat[i][17:19])
        sec[i] = int(float(time_sat[i][20:31]))
                
        Vx[i] = float(vel_sat[i][4:18]) * (10**(-4))
        Vy[i] = float(vel_sat[i][18:32]) * (10**(-4))
        Vz[i] = float(vel_sat[i][32:46]) * (10**(-4))
        
        x[i] = float(pos_sat[i][4:18])
        y[i] = float(pos_sat[i][18:32])
        z[i] = float(pos_sat[i][32:46])
    
    time_array_gps = np.array([start_sp3 + \
                    datetime.timedelta(seconds = i*1) for i in range(0,k)])
    
    time_array_ut = np.array([start_sp3 + \
                    datetime.timedelta(seconds = i*1) - \
                       datetime.timedelta(seconds = start_sp3.second) \
                           for i in range(0,k)])
      
    # start row
    srow = int(np.where(time_array_ut == start_date)[0])
    # end row    
    erow = int(np.where(time_array_ut == end_date)[0]) + 1
    
    time_experiment = time_array_ut[srow:erow]
             
    length = erow-srow
    
    ITRFVx = np.empty(length); ITRFVy = np.empty(length);ITRFVz = np.empty(length)    
    ITRFx = np.empty(length); ITRFy = np.empty(length); ITRFz = np.empty(length)    
 
    ITRFx = x[srow:erow]
    ITRFy = y[srow:erow]
    ITRFz = z[srow:erow]
    
    ITRFVx = Vx[srow:erow]
    ITRFVy = Vy[srow:erow]
    ITRFVz = Vz[srow:erow]
    
    return { 'srow':srow, 'erow':erow, 'time_array_gps': time_array_gps,
            'time_array_ut': time_array_ut, 'time_experiment': time_experiment,
            'ITRFx': ITRFx, 'ITRFy': ITRFy, 'ITRFz': ITRFz, 
            'ITRFVx': ITRFVx, 'ITRFVy': ITRFVy, 'ITRFVz': ITRFVz, 
        } 
            
def calculate_orbital_elements(pX, pY, pZ, Vx, Vy, Vz):
    """
    Calculates satellite orbital parameters using X, Y, Z, Vx, Vy, Vz in GEIJ2K
    
    Reference: Curtis, H. D. (2014). Orbits in Three Dimensions in 
    Orbital mechanics for engineering students. Butterworth-Heinemann.

    Parameters
    ----------
    pX : numpy.ndarray[float]
        X position in GEIJ2K (km).
    pY : numpy.ndarray[float]
        Y position in GEIJ2K (km).
    pZ : numpy.ndarray[float]
        Z position in GEIJ2K (km).
    Vx : numpy.ndarray[float]
        X component of velocity in GEIJ2K (km/s).
    Vy : numpy.ndarray[float]
        Y component of velocity in GEIJ2K (km/s).
    Vz : numpy.ndarray[float]
        Z component of velocity in GEIJ2K (km/s).

    Returns
    -------
    inc : float
        satellite inclination (radian).
    raan : float
        Right ascension of ascending node (radian).
    ap : float
        Argument of periapsis (radian).
    e : float
        eccentricity (radian).
    TA : float
        True anomaly (radian).
        
    Examples
    --------
    sat_inc, sat_ap, sat_raan, sat_ecc, sat_TA = \
        
        ei.calculate_orbital_elements(dict_rri['GEIx'], dict_rri['GEIy'], \
                                      
                                      dict_rri['GEIz'], dict_rri['GEIVx'], \
                                      
                                      dict_rri['GEIVy'], dict_rri['GEIVz'])
    """
          
    inc = np.empty(pX.shape); inc.fill(np.nan)
    raan = np.empty(pX.shape); raan.fill(np.nan)
    ap = np.empty(pX.shape); ap.fill(np.nan)           
    TA = np.empty(pX.shape); TA.fill(np.nan)
    E = np.empty(pX.shape); E.fill(np.nan)
    Praan = np.empty(pX.shape); Praan.fill(np.nan)
    
    for i in range (0,len(pX)):
        # distance
        d = np.sqrt(pX[i]**2 + pY[i]**2 + pZ[i]**2)
        d_vec = [pX[i], pY[i], pZ[i]]
        # speed
        V = np.sqrt(Vx[i]**2 + Vy[i]**2 + Vz[i]**2)
        V_vec = [Vx[i], Vy[i], Vz[i]]
        
        # Vrad < 0: satellite is flying towards perigee
        # Vrad > 0: satellite is flying away from perigee
        Vrad = (pX[i]*Vx[i] + pY[i]* Vy[i] + pZ[i]* Vz[i]) / d
        
        # components of angular momentum vector
        hx = pY[i] * Vz[i] - Vy[i] * pZ[i]
        hy = Vx[i] * pZ[i] - pX[i] * Vz[i]
        hz = pX[i] * Vy[i] - Vx[i] * pY[i]
        h_vec = [hx, hy, hz]
        
        h = np.linalg.norm(h_vec)
        
        # satellite inclination is between 0 and 180. 
        # if 90 < inc < = 180: retrograde orbit
        inc[i] = np.arccos ( hz / h )
        
        # node line
        N_vec = np.cross([0, 0, 1], h_vec)
        
        N = np.linalg.norm(N_vec)
        
        # nx/N > 0: raan is in first or fourth quadrant
        # nx/N < 0: raan is in the second or third quadrant
        # if ny > 0: ascending node lies on the positive side of the vertical XZ
        #            0- pi
        # if ny < 0: ascending node lies on the negative side of the vertical XZ
        #            pi- 2* pi
        Praan[i] = np.arccos(N_vec[0] / N)    
        if (Praan[i] > 0 and N_vec[1] < 0):
            raan[i] = 2*np.pi - Praan[i]
        elif (Praan[i] <= 0 and N_vec[1] < 0):
            raan[i] = 2*np.pi - Praan[i]
        elif (Praan[i] <= 0 and N_vec[1] > 0):
            raan[i] = Praan[i]            
        else:
            raan[i] = Praan[i]
        
        # gravitational parameter mu: G * (m1 + m2) 
        # G: 6.67408 × 10-11 m3 kg-1 s-2
        # m1: 5.972 × 10^24 kg
        # m2:  500 kg 
        mu = 6.67408 * 10**(-11) * ((5.972 * 10**(24) ) + 500) 
        # convert to km^3: because velocity's unit is [km per s]
        mu = mu * 10 **(-9)
        
        # eccentricity
        e_d = (1/mu) * (V**2 - mu/d)  # vector component along r
        e_v = (-1/mu) * (d * Vrad)    # vector component along V
        
        e_vec_x = e_d * pX[i] + e_v * Vx[i]
        e_vec_y = e_d * pY[i] + e_v * Vy[i]
        e_vec_z = e_d * pZ[i] + e_v * Vz[i]
        e_vec = [e_vec_x, e_vec_y, e_vec_z]
        
        e = np.sqrt(e_vec_x**2 + e_vec_y**2 + e_vec_z**2)
        # calculate the True anomaly
        P1 = np.dot(e_vec,d_vec)/(e*d)
        if np.dot(d_vec, V_vec)< 0:
            TA[i] = 2*np.pi - np.arccos(P1)
        else:
            TA[i] = np.arccos(P1)
            
        # eccentric anomaly
        E[i] = 2*np.arctan2( ( 1-e)**(1/2) , (1+e)**(1/2) )*np.tan(TA[i]/2)
        # else:
        #     cp = np.cross(N_vec,d_vec)
        #     P2 = np.dot(N_vec, d_vec)/(N*d)
        #     if cp[3] >= 0:          
        #         TA = np.arccos(P2)
        #     else:
        #         TA = 2*np.pi - np.arccos(P2)
           
        # # another way to calculate |e|: yields the same results as above
        # Part_1 = (2 * mu - (d * ( V**2) ))* d * (Vrad **2)
        # Part_2 = (mu - d* (V**2))**2
        # e2 = (1/mu) * np.sqrt(Part_1 + Part_2)
           
        # argument of perigee 
        Ne_dot = np.dot(N_vec, e_vec)
        Ne = N*e
        P2 = np.arccos( Ne_dot / Ne)
        # if Ne_dot > 0 ; ap is in the first or fourth quadrant
        # if Ne_dot < 0 ; ap is in the second or third quadrant
        # if e is in positive Z direction, pointing up; ap between 0 and pi
        # if e is in negative Z direction, pointing down; ap between pi and 2pi
        if (e_vec_z < 0):
            ap[i] = (2 * np.pi) - P2
        else:
            ap[i] = P2
  
    return inc, ap, raan, e, TA


def import_tle(filename, filedate, DOY):
    """
    Reads Celestrak TLE file 
    
    File obtained from: https://celestrak.com/   

    Parameters
    ----------
    filename: str
        filename for TLE file.
    filedate: str
        date of RRI passage.
    DOY: int 
        day of year of RRI passage.

    Returns
    -------
    tle_epoch : float
        epoch of TLE file.
    tle_inc : float
       satellite inclination (degrees).
    tle_ap : float
        satellite argument of perigee (degrees).
    tle_raan : float
        right ascension of ascending node (degrees).
    tle_ecc : int
        orbit eccentricity.
    tle_ma : float
        mean anomaly (degrees).
        
    Examples
    --------
    time_start = datetime.datetime(2016, 4, 18, 22, 27, 59) 
    
    DOY = time_start.timetuple().tm_yday
    
    filedate = '20160418'    
    
    tle_epoch, tle_inc, tle_ap, tle_raan, tle_ecc, tle_ma = \
        
        ei.import_tle(file_TLE, filedate, DOY)

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
            if l1[18:24] == sat_date:
                r= i
                a += 1
                line1.append(l1)

        elif row.startswith('2 '):
            l2 = row
            j += 1
            # line 2 for the specific date is when the index is equal to 
            # sat_date index
            if j == r:
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

def compare_orbital(file_TLE, filedate, DOY, pX, pY, pZ, Vx, Vy, Vz):
   
    """   
    Compares values from the TLE file with the calculated average TLE values.
    Prints out the comparison results.
    
    inc: spacecraft inclination angle
    raan: right ascension of ascending node
    ap: argument of periapsis
    
    Parameters
    ----------
    file_TLE: str
        filename for TLE file (including path).
    filedate: str
        date of RRI passage, format: '%H%M%S'.
    DOY: int 
        day of year of RRI passage. 
    pX : numpy.ndarray[float]
        X position in GEIJ2K (km).
    pY : numpy.ndarray[float]
        Y position in GEIJ2K (km).
    pZ : numpy.ndarray[float]
        Z position in GEIJ2K (km).
    Vx : numpy.ndarray[float]
        X component of velocity in GEIJ2K (km/s).
    Vy : Tnumpy.ndarray[float]
        Y component of velocity in GEIJ2K (km/s).
    Vz : numpy.ndarray[float]
        Z component of velocity in GEIJ2K (km/s).

    Returns
    -------
    sap_comp : list[float, float, float]
        
        Element1: the difference between the tle and computed inclination (deg) 
        
        Element2: the difference between the tle and computed raan (deg) 
        
        Element3: the difference between the tle and computed ap (deg) 
    
    Examples
    --------
    time_start = datetime.datetime(2016, 4, 18, 22, 27, 59) 
    
    DOY = time_start.timetuple().tm_yday
    
    filedate = '20160418'
    
    comp = ei.compare_orbital(file_TLE, filedate, DOY, dict_rri['GEIx'], 
                             
                              dict_rri['GEIy'], dict_rri['GEIz'], 
                             
                              dict_rri['GEIVx'], dict_rri['GEIVy'], 
                              
                              dict_rri['GEIVz'])
        
    """
    _, tle_inc, tle_ap, tle_raan, *_ = \
        import_tle(file_TLE, filedate, DOY)
        
    sat_inc, sat_ap, sat_raan, *_ = \
        calculate_orbital_elements(pX, pY, pZ, Vx, Vy, Vz)
        
    if np.size(tle_inc) > 0:
        print(" ".join(['for the satellite epoch', filedate, ':\n']))
        sap_comp = np.empty(3)
        
        sinc = np.rad2deg(np.mean(sat_inc))
        sap_comp[0] = tle_inc - sinc
        if np.isclose(tle_inc, sinc , atol=1) == True:
            print(" ".join(['TLE_inc =',str(tle_inc),', sat_inc =', str(sinc), 
                            ': calculations are correct\n']))
        else: 
            print('problem with satellite inclination: \n',
                  'check computed orbital elements,',
                  'the tle or the ephemeris file\n')
        
        sarn = np.rad2deg(np.mean(sat_raan))
        sap_comp[1] = tle_raan - sarn
        if np.isclose(tle_raan, sarn, atol=1) == True:
            print(" ".join(['TLE_raan =',str(tle_raan),', sat_raan =', 
                            str(sarn), ': calculations are correct\n']))
        else: 
            print(" ".join(['TLE_raan =',str(tle_raan),', sat_raan =', 
                            str(sarn), 'problem with satellite raan:',
                            '\n check computed orbital elements,',
                            'the tle or the ephemeris file']))
        
        sat_ap_deg = np.rad2deg(np.mean(sat_ap))
        sap_comp[2] = tle_ap - sat_ap_deg
        if np.isclose(tle_ap, sat_ap_deg, atol=1) == True:
            print(" ".join(['TLE_ap =',str(tle_ap),', sat_ap =', 
                            str(sat_ap_deg), ': calculations are correct\n']))
        else: 
            print(" ".join(['TLE_ap =',str(tle_ap),', sat_ap =', 
                            str(sat_ap_deg),'problem with satellite ap:',
                            '\n check computed orbital elements,',
                            'the tle or the ephemeris file\n']))
    else:
        print('No TLE to compare with computed orbital elements')
                
    return sap_comp
