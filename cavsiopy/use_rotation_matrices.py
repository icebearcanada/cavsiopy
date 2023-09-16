#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
use_rotation_matrices contains functions to  

1. calculate the GMST time from noon or midnight

2. compute the rotation matrices for the transformations between 
spacecraft orbital frame(ORF) and GEI J2000, J2K and ECEF, J2K and ITRF, 
ECEF and NED, ITRF and NEC, ITRF and NED, and  NED and ENU.
    
Two sets of rotation matrices are provided to build up direction cosine matrix.

1. SET-1:
According to Coordinate Transformations via Euler Angles, Riggs, 2019. Rev. E.
expresses the inertial frame vectors in terms of rotated frame vectors
used for transformation from a rotated frame to inertial frame
gives i, j, k in terms of i', j', k'; abbreviated as r2i in the code
definition of the reference: x (out), y (right), and z (up)

ROLL: positive is from +Y towards +Z
PITCH: positive is from +Z towards +X
YAW: positive is from +X towards +Y

2. SET-2:
According to Curtis, H. D. (2013). Orbital mechanics for engineering students.
expresses the rotated frame vectors in terms of inertial frame vectors
for transformation from an inertial frame to rotated frame
gives i', j', k' in terms of i, j, k; abbreviated as i2r in the code
RRI description: x (out), y (right), and z (down)

ROLL: positive is from +Y towards +Z
PITCH: positive is from +Z towards +X
YAW: positive is from +X towards +Y

.. toctree::
  :maxdepth: 2
  GMST_midnight
  GMST_noon
  RX_i2r
  RX_r2i
  RY_i2r
  RY_r2i
  RZ_i2r
  RZ_r2i
  ecef_to_gei_rm
  gei2ecef
  gei_to_ecef_rm
  gei_to_oe_rm
  icrf2itrf
  icrf_to_itrf_rm
  itrf_to_icrf_rm
  j2k_to_orf_rm
  nec2terrestrial
  nec_to_terrestrial_rm
  ned2enu
  ned_to_terrestrial_rm
  oe_to_gei_rm
  orf_to_j2k_rm
  orf_to_j2k_use_orbital_elements
  orf_to_j2k_use_spacecraft_ephemeris
  terrestrial2nec
  terrestrial2ned
  terrestrial_to_nec_rm
  terrestrial_to_ned_rm
    
@author: ceren, warren
"""
import numpy as np
from astropy.time import Time
import pysofa2 as pysofa
import cavsiopy.ephemeris_importer as ei
import cavsiopy.miscellaneous as misc
import cavsiopy.complement_missing_sofa as cms


def GMST_midnight(utc_dt):
    """
    Calculates GMST in radians.
    
    Includes day fraction corrected with respect to midnight.
    
    Follows: Curtis, H. D. (2013). Orbital mechanics for engineering students.
    Butterworth-Heinemann. This formula gives GMST in degrees.
    
    Parameters
    ----------
    utc_dt: datetime.datetime
        time in UTC.
    
    Returns
    -------
    GMST: float
        GMST in radians.
    
    Examples
    --------
    GMST = calculate_GMST_midnight(utc_dt_array)
    """

    JD = Time( utc_dt ).ut1.jd
    T  = (JD - 2451545)/36525 # Julian Centuries since J2K epoch in TT.
    UT = JD % 1 - 0.5   # Day fraction, correcting to midnight rather than noon.
  
    gmst = 100.460618375 + 36000.770053608336*T + 0.0003879333*T**2 -\
        2.583e-8*T**3
  
    E_deg = 360.98564724 * UT # /24 built-in, UT is just the fractional day.
  
    theta_deg = gmst + E_deg
    GMST = ( theta_deg * np.pi/180 ) % (2*np.pi)
    
    return GMST

def GMST_noon(time_array):
    """
    Follows: Curtis, H. D. (2013). Orbital mechanics for engineering students.
    Butterworth-Heinemann. This formula gives GMST in degrees.

    Includes day fraction corrected with respect to noon.

    Parameters
    ----------
    time_array : datetime.datetime
        time in UTC.

    Returns
    -------
    GMST: float
        GMST in radians.

    Examples
    --------
    GMST = calculate_GMST_noon(time_array)
    """    
   
    # time of epoch, jd: Julian date
    JD=Time(time_array).ut1.jd
    # Time elapsed since JD 2451545 (Jan 1, 2000, 12 UT)
    T = (JD - 2451545)/36525;  

    GMST_deg = 100.460618375 + 36000.770053608336 * T + \
        0.0003879333 * (T**2) - 2.583 * (10**(-8)) * (T**3)

    # Formula above can yield values greater than 360. If that's the case:
    if GMST_deg < 0:
        GMST = GMST_deg + np.fix(GMST_deg/360)*360
    elif GMST_deg >= 360:
        GMST = GMST_deg - (np.fix(GMST_deg/360) - 1)*360
    else:
        GMST = GMST_deg

    # from Earth's rotation
    UT_H = time_array.timetuple().tm_hour +\
        time_array.timetuple().tm_min/60 +\
            ( time_array.timetuple().tm_sec / 3600 )

    E_deg = 360.98564724 *  (UT_H/24)

    theta_deg = GMST + E_deg

    theta_rad = np.deg2rad(theta_deg)

    return theta_rad

def RX_r2i(x):
    """
    Parameters
    ----------
    x : float
        roll angle (radians).

    Returns
    -------
    numpy.ndarray
        Rotation matrix about x axis.
        
    Examples
    --------
    Rx = RX_r2i(roll)
    """
    return np.array([[1,      0,       0],
                      [0, np.cos(x), -np.sin(x)],
                      [0, np.sin(x),  np.cos(x)]])

def RY_r2i(y):
    """
    Parameters
    ----------
    y : float
        pitch angle (radians).

    Returns
    -------
    numpy.ndarray
        Rotation matrix about y axis.
        
    Examples
    --------
    Ry = RY_r2i(pitch)
    """
    return np.array([[np.cos(y), 0, np.sin(y)],
                      [0,      1,        0],
                      [-np.sin(y), 0, np.cos(y)]])

def RZ_r2i(z):
    """
    Parameters
    ----------
    z : float
        yaw angle (radians).

    Returns
    -------
    numpy.ndarray
        Rotation matrix about z axis.
    
    Examples
    --------
    Rz = RZ_r2i(yaw)
    """
    return np.array([[np.cos(z), -np.sin(z), 0],
                      [np.sin(z), np.cos(z), 0],
                      [0,      0,       1]])

def RX_i2r(x):
    """
    Parameters
    ----------
    x : float
        roll angle (radians).

    Returns
    -------
    numpy.ndarray
        Rotation matrix about x axis.

    Examples
    --------
    Rx = RX_i2r(roll)        
        
    """
    return np.array([[1,      0,       0],
                      [0, np.cos(x), np.sin(x)],
                      [0, -np.sin(x),  np.cos(x)]])

def RY_i2r(y):
    """
    Parameters
    ----------
    y : float
        pitch angle (radians).

    Returns
    -------
    numpy.ndarray
        Rotation matrix about y axis.
        
    Examples
    --------
    Ry = RY_i2r(pitch)  
    """
    return np.array([[np.cos(y), 0, -np.sin(y)],
                      [0,      1,        0],
                      [np.sin(y), 0, np.cos(y)]])    
        
def RZ_i2r(z):
    """
    Parameters
    ----------
    z : float
        yaw angle (radians).

    Returns
    -------
    numpy.ndarray
        Rotation matrix about z axis.
        
    Examples
    --------
    Rz = RZ_i2r(yaw)
    """
    return np.array([[np.cos(z), np.sin(z), 0],
                      [-np.sin(z), np.cos(z), 0],
                      [0,      0,       1]])

def oe_to_gei_rm(raan, inc, ap):
    """
    Function to calculate the perifocal frame (OE) to GEI rotation matrix
    
    To transform from OE to J2K using orbital elements:
    RZ_r2i(raan) @ RX_r2i(inc) @ RZ_r2i(ap)
    
    Ref: Frame rotations and quaternions, Pedro A. Capó-Lugo, Peter M. Bainum,
    in Orbital Mechanics and Formation Flying, 2011

    Parameters
    ----------
    inc : float
        satellite inclination (in radians).
    raan : float
        satellite right ascension of ascending node (in radians).
    ap : float
        satellite argument of periapsis (in radians).

    Returns
    -------
    rm_OE2GEI : numpy.ndarray[float]
        ORF to GEI matrix using orbital elements.
        
    Examples
    --------
    R = oe_to_gei(sat_raan, sat_inc, sat_ap)

    """
    rm_OE2GEI = RZ_r2i(raan) @ RX_r2i(inc) @ RZ_r2i(ap)
   
    return rm_OE2GEI


def gei_to_oe_rm(inc, raan, ap):
    """
    Function to calculate the GEI to orbital frame matrix.
    
    GEI2OE matrix is the transpose of OE2GEI matrix.

    Parameters
    ----------
    inc : float
        satellite inclination (in radians).
    raan : float
        satellite right ascension of ascending node (in radians).
    ap : float
        satellite argument of periapsis (in radians).

    Returns
    -------
    rm_GEI2OE : numpy.ndarray[float]
        orbital frame to GEI matrix using satellite orbital 
        elements.
        
    Examples
    --------
    R = gei_to_oe(sat_raan, sat_inc, sat_ap)

    """
    
    OE2GEI = RZ_r2i(raan) @ RX_r2i(inc) @ RZ_r2i(ap)
    rm_GEI2OE = OE2GEI.T
    
    return rm_GEI2OE

def orf_to_j2k_rm(pX, pY, pZ, Vx, Vy, Vz):
    """   
    Given 6 floats of the J2K/ICRF position and velocity data,
    Returns a 3x3 rotation matrix describing the rotation of ORF -> J2K.
    
    Ref: 
    Canuto, E., Novara, C., Carlucci, D., Montenegro, C. P., & Massotti, L. 
    (2018). Orbital Control and Prediction Problems. In Spacecraft Dynamics 
    and Control: The Embedded Model Control Approach. Butterworth-Heinemann.


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
    Vy : Tnumpy.ndarray[float]
        Y component of velocity in GEIJ2K (km/s).
    Vz : numpy.ndarray[float]
        Z component of velocity in GEIJ2K (km/s).

    Returns
    -------
    rNadir_to_J2K : np.ndarray[float]
        rotation matrice for the transformation from ORF to GEIJ2K

    """    
    # Rebuild the Position and Velocity vectors in ICRF/J2K
    Position = np.array([pX, pY, pZ])
    Velocity = np.array([Vx, Vy, Vz])
    
    # Invert vector, as Z is satellite-down rather than center-of-earth-up 
    # as given in the position.
    zVector = -Position  
    yVector = np.cross(zVector, Velocity)
    xVector = np.cross(yVector, zVector)
    
    # Normalize the vectors:
    xVector = xVector / np.linalg.norm(xVector)
    yVector = yVector / np.linalg.norm(yVector)
    zVector = zVector / np.linalg.norm(zVector)
    
    # Build the rotation matrix, which due to the side-by-side construction,
    # is Nadir->J2K
    rNadir_to_J2K = np.dstack([xVector, yVector, zVector])[0]
    
    return rNadir_to_J2K

def j2k_to_orf_rm(pX, pY, pZ, Vx, Vy, Vz):
    """
    Given 6 floats of the J2K/ICRF position and velocity data,
    Returns a 3x3 rotation matrix describing the rotation of J2K-> ORF.
    
    Uses the transpose of rNadir_to_J2K.
    
    Ref: 
    Canuto, E., Novara, C., Carlucci, D., Montenegro, C. P., & Massotti, L. 
    (2018). Orbital Control and Prediction Problems. In Spacecraft Dynamics 
    and Control: The Embedded Model Control Approach. Butterworth-Heinemann.


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
    Vy : Tnumpy.ndarray[float]
        Y component of velocity in GEIJ2K (km/s).
    Vz : numpy.ndarray[float]
        Z component of velocity in GEIJ2K (km/s).

    Returns
    -------
    rJ2K_to_Nadir : np.ndarray[float]
        rotation matrice for the transformation from GEIJ2K to ORF.

    """    
    rNadir_to_J2K = orf_to_j2k_rm(pX, pY, pZ, Vx, Vy, Vz)
    
    # Transpose to get the J2K->Nadir rotation matrix.
    rJ2K_to_Nadir = rNadir_to_J2K.T
    
    return rJ2K_to_Nadir


def gei_to_ecef_rm(theta_rad):
    """ 
    Function to calculate the GEI to ECEF rotation matrix.
    
    One rotation is needed in xy frame (yaw rotation) to transform from 
    GEIJ2K to ECEF. Yaw angle is the azimuth angle between GEIJ2K and 
    ECEF in this case.
    
    Azimuth= Greenwich sidereal time + Earth's rotation speed * UTC.
    
    Ref:
    1. J. Riggs, Coordinate transformations via Euler Angle Rotations,
    Rev. E - 10 April 2019
    
    2. A. Sofyali, Orbital Mechanics notes, 7 February 2019.    

    Parameters
    ----------
    theta_rad : float
        GMST_midnight.

    Returns
    -------
    rm_GEI2ECEF : numpy.ndarray[float]
        Rotation matrix for transformations from GEI to ECEF.

    Examples
    --------
    R = gei_to_ecef_rm(theta_rad)

    """
    
    rm_GEI2ECEF = RZ_i2r(theta_rad)
    
    return rm_GEI2ECEF

def ecef_to_gei_rm(theta_rad): 
    """
    Function to calculate ECEF to GEI rotation matrix.
    
    Uses the transpose of GEIJ2K to ECEF rotation matrix.
    
    Parameters
    ----------
    theta_rad : float
        GMST_midnight.
    
    Returns
    -------
    GEI2ECEF : numpy.ndarray[float]
        Rotation matrix for transformations from GEI to ECEF.
        
    Examples
    --------
    R = ecef_to_gei_rm(theta_rad)
    
    """
    rm_ECEF2GEI = (RZ_i2r(theta_rad)).T
    
    return rm_ECEF2GEI
    
def terrestrial_to_ned_rm(lat, lon):
    """
    Function to calculate the ECEF to NED rotation matrix.
    
    Y-Z (pitch-yaw) sequence with longitude and latitude are needed.
    
    1. Aligning the x-y plane of ECEF and NED using longitude: yaw(longitude)
    
    2. Aligning the x-z plane using latitude: pitch(-(pi/2+latitude))
    
    Ref: 
    1. J. Riggs, Coordinate transformations via Euler Angle Rotations,
    Rev. E - 10 April 2019
    
    2. Cai, G., Chen, B.M., Lee, T.H. (2011). Coordinate Systems and
    Transformations. In: Unmanned Rotorcraft Systems. Advances in Industrial
    Control. Springer, London. https://doi.org/10.1007/978-0-85729-635-1_2

    Parameters
    ----------
    lat : numpy.ndarray[float]
        Geodetic latitude in degrees.
    lon : numpy.ndarray[float]
        Geodetic longitude in degrees.

    Returns
    -------
    rm_ter2ned : numpy.ndarray[float]
        Rotation matrix for transformations from ECEF to NED.

    """
    # Open form ----------------------------------------------------------    
    lamda = np.deg2rad(lon)
    psi = np.deg2rad(lat)
    a = -np.cos(lamda)*np.sin(psi)
    b = -np.sin(lamda)
    c = -np.cos(lamda)*np.cos(psi)
    d = -np.sin(lamda)*np.sin(psi)
    e = np.cos(lamda)
    f = -np.sin(lamda)*np.cos(psi)
    g = np.cos(psi)
    h = 0
    m = -np.sin(psi)

    rm_ter2ned = np.array([[a, d, g], 
                            [b, e, h],
                            [c, f,  m]])
    
    return rm_ter2ned

def ned_to_terrestrial_rm(lat, lon):
    """
    Function to calculate NED to ECEF rotation matrix.
    
    Uses the transpose of ECEF to NED rotation matrix.

    Parameters
    ----------
    lat : numpy.ndarray[float]
        Geodetic latitude in degrees.
    lon : numpy.ndarray[float]
        Geodetic longitude in degrees.

    Returns
    -------
    rm_ned2ter : numpy.ndarray[float]
        Rotation matrix for transformations from NED to ECEF.

    """
    # Open form ----------------------------------------------------------    
    ter2ned_rm = terrestrial_to_ned_rm(lat, lon)
    
    rm_ned2ter = ter2ned_rm.T
    
    return rm_ned2ter

def icrf_to_itrf_rm(path_to_files, input_time):
    """
    Function to calculate the ICRF to ITRF matrix using pysofa routines.

    Parameters
    ----------
    path_to_files : str
        path_to_initialization files (IERS and EOP).
    input_time : datetime.datetime
        time.

    Returns
    -------
    rm_ICRF2ITRF : numpy.ndarray[float]
        ICRF to ITRF rotation matrix.

    """
    # calculate day-of-year = DOY
    start_date = input_time
    DOY = start_date.timetuple().tm_yday
    date_year = start_date.year
    if (date_year%4 == 0):
        doy = DOY / 366
        # calculate year + doy to match with IERS data
        ydoy = date_year + doy
    # year in IERS tabulations for TT-UT1 changes with 0.05
        inc = 0.05 * 366
    else:
        doy = DOY / 365.25
    # calculate year + doy to match with IERS data
        ydoy = date_year + doy
    # year in IERS tabulations for TT-UT1 changes with 0.05
        inc = 0.05 * 365.25

    # find the terrestrial time difference from universal time---------------

    # obtain UT1-TAI from
    # https://datacenter.iers.org/data/latestVersion/38_EOP_C01.1900-NOW_V2013_0138.txt
    filename_UT1_TAI = '38_EOP_C01.1900-NOW_V2013_0138.txt'
    file_UT1_TAI = path_to_files + filename_UT1_TAI
    data_UT1_TAI = np.loadtxt(fname=file_UT1_TAI, skiprows=1, comments='#')
    year_UT1_TAI = data_UT1_TAI[:,0]
    UT1_TAI = data_UT1_TAI[:,5]

    index_year = misc.find_index(year_UT1_TAI, ydoy, 0.025)

    if np.size(index_year)==1:
        UT1_TAI_date = UT1_TAI[index_year[0][0]]
    elif np.size(index_year)>1:
        UT1_TAI_distance_to_lower = ydoy- year_UT1_TAI[index_year[0][0]]
        UT1_TAI_lower = UT1_TAI[index_year[0][0]]

        UT1_TAI_distance_to_upper = year_UT1_TAI[index_year[0][1]]-ydoy
        UT1_TAI_upper = UT1_TAI[index_year[0][1]]
        # calculate the weighted_mean to find deltaT of date
        UT1_TAI_date = UT1_TAI_lower * (UT1_TAI_distance_to_lower/0.05) + \
            UT1_TAI_upper * (UT1_TAI_distance_to_upper/0.05)

    # TT-UT1 = 32.184s - (UT1-TAI)
    # Ref: https://hpiers.obspm.fr/eop-pc/earthor/ut1lod/ut1-tai_pred.html
    dt = 32.184 - (UT1_TAI_date)
    
    JD = Time(input_time).jd  # jd: Julian date
    # Time elapsed since JD 2451545 (Jan 1, 2000, 12 UT)
    dT_JD = (JD - 2451545)/36525;  

    utc1 = np.floor(JD)
    utc2= JD-utc1

    MJD = int(np.floor(Time(input_time).mjd))

    #  IERS tabulations
    #  delta UT (UT1-UTC) can be obtained from the link below
    # link = 
    # "https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt"
    filename_IERS = 'IERS_daily_deltaT.txt'
    file_IERS = path_to_files + filename_IERS
    IER=np.loadtxt(fname=file_IERS, skiprows=14, comments='#')

    # find where MJD equals the MJD in the text file
    index_day = np.where(IER== MJD)[0]
    # delta UT from IERS file
    dut = IER[index_day, 7]

    # find the terrestrial time from coordinated universal time--------
    #  dut = ut1-utc
    uta, utb = pysofa.Utcut1(utc1, utc2, dut)

#----------------------------------------------------------------------
    #  dt = ut-tt
    tta, ttb = pysofa.Ut1tt(uta, utb, dt)

    #  find the x_p, y_p, s: celestial pole coordinates----------------
    x_p, y_p, s = pysofa.Xys00b(tta, ttb)

    #  find the celestial to terrestrial matrix------------------------
    rm_ICRF2ITRF = cms.c2t00b(tta, ttb, uta, utb, x_p, y_p)
    
    return rm_ICRF2ITRF

def itrf_to_icrf_rm(path_to_files, time_array):
    """
    Function to calculate ITRF to ICRF rotation matrix using pysofa routines.
    
    Uses the transpose of ITRF to ICRF rotation matrix.

    Parameters
    ----------
    path_to_files : str
        path_to_initialization files (IERS and EOP).
    time_array : datetime.datetime
        time array.

    Returns
    -------
    rm_ICRF2ITRF : numpy.ndarray[float]
        ICRF to ITRF rotation matrix.

    """
    
    rm_ICRF2ITRF = icrf_to_itrf_rm(path_to_files, time_array)
    
    rm_ITRF2ICRF = rm_ICRF2ITRF.T
    
    return rm_ITRF2ICRF

def terrestrial_to_nec_rm( gX, gY, gZ ):
    """
    Given the Geographic Cartesian Position ephemeris data return a 
    3d rotation matrix describing the rotation from ITRF to NEC
    
    Ref:
    https://earth.esa.int/documents/10174/1514862/Swarm_Level-1b_Processor_Algorithms
    Page 137, "C.3 North East Center (NEC)

    Parameters
    ----------
    gX : numpy.ndarray[float]
        X position in GEO/ITRF (km).
    gY : numpy.ndarray[float]
        Y position in GEO/ITRF (km).
    gZ : numpy.ndarray[float]
        Z position in GEO/ITRF (km).

    Returns
    -------
    R_NEC : numpy.ndarray[float]
        3d rotation matrix describing the rotation from ITRF to NEC.
      
    Examples
    --------
    rter2nec  = terrestrial_to_nec_rm(GEOx, GEOy, GEOz)
    """

    Position_Vector = np.array([ gX,gY,gZ ])
    
    eCenter = -Position_Vector
    # The cross product of geographic north-pole vector and the center vector.
    eEast   = np.cross( eCenter, [0,0,1] ) 
    eNorth  = np.cross( eEast, eCenter )
    
    # And then normalize the vectors:
    eCenter /= np.linalg.norm(eCenter)
    eEast   /= np.linalg.norm(eEast  )
    eNorth  /= np.linalg.norm(eNorth )
    
    ter2nec_rm = np.stack(( eNorth, eEast, eCenter))
    
    return  ter2nec_rm

def nec_to_terrestrial_rm(gX, gY, gZ):
    """
    Function to calculate NEC to ITRF rotation matrix.
    
    Uses the transpose of ITRF to NEC rotation matrix.    

    Parameters
    ----------
    gX : numpy.ndarray[float]
        X position in GEO/ITRF (km)
    gY : numpy.ndarray[float]
        Y position in GEO/ITRF (km)
    gZ : numpy.ndarray[float]
        Z position in GEO/ITRF (km)

    Returns
    -------
    rm_nec2ter : numpy.ndarray[float]
        3d rotation matrix describing the rotation from NEC to ITRF.
        
    Examples
    --------
    r_nec2ter = nec_to_terrestrial_rm(GEOx, GEOy, GEOz)

    """
    
    rm_ter2nec = terrestrial_to_nec_rm( gX, gY, gZ )
    
    rm_nec2ter = rm_ter2nec.T
    
    return rm_nec2ter

def orf_to_j2k_use_orbital_elements(body_vec, pX, pY, pZ, Vx, Vy, Vz,
                                 P= np.array([[0, 0, -1],[1, 0, 0],[0,-1,0]])):
    """
    Employs orbital elements to transform from ORF to GEIJ2K: 
    
    To align the perifocal frame with the orbital frame 
    RZ_r2i(TA) @ P @  body_vec
    
    Ref:Canuto, E., Novara, C., Carlucci, D., Montenegro, C. P., & Massotti, L. 
    (2018). Orbital Control and Prediction Problems. In Spacecraft Dynamics 
    and Control: The Embedded Model Control Approach. Butterworth-Heinemann

    Parameters
    ----------
    body_vec : np.ndarray[float]
        rotated body frame vector.
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
    P  : numpy.ndarray
        optional reordering matrix to align the perifocal frame with the 
    orbital frame.The default is P= np.array([[0, 0, -1],[1, 0, 0],[0,-1,0]])

    Returns
    -------
    inst_GEI : numpy.ndarray[float]
        instrument look direction in GEIJ2K.
        
    Examples
    --------
    inst_GEI = \
        
        orf_to_j2k_use_orbital_elements(body_vec, pX, pY, pZ, Vx, Vy, Vz)

    """

    # Initialize GEI/ECI Vectors
    SC_OE2GEI = np.empty(body_vec.shape)
    SC_OE2GEI2 = np.empty(body_vec.shape)

    #  obtain orbital elements
    inc, ap, raan, e, TA = ei.calculate_orbital_elements(pX, pY, pZ, \
                                                         Vx, Vy, Vz)

    for i in range(0, body_vec.shape[1]):
        # align the orbital frame with the perifocal frame
        inst_OE = RZ_r2i(TA[i]) @ P @  body_vec[:, i]
        # rotation matrix to transform from perifocal to GEIJ2K
        OE2GEI = oe_to_gei_rm(raan[i], inc[i], ap[i])

        SC_OE2GEI[:, i] = OE2GEI @ (inst_OE)
        #  normalize
        denominator = np.sqrt(SC_OE2GEI[0, i]**2 + SC_OE2GEI[1, i]**2 + \
                              SC_OE2GEI[2, i]**2)
            
        SC_OE2GEI2[:, i] = SC_OE2GEI[:, i] / denominator

    inst_GEI = np.array(SC_OE2GEI2)

    return inst_GEI

# =============================================================================
# # %% function for GEI to ECEF transformations
# =============================================================================
def gei2ecef(inst_GEI, time_array):
    """
    function for GEIJ2K to ECEF transformations

    Parameters
    ----------
    inst_GEI : np.ndarray[float]
        instrument look direction in GEIJ2K .  
    time_array : datetime.datetime
        experiment time interval.

    Returns
    -------
    inst_ECEF : np.ndarray[float]
        instrument look direction in GEIJ2K.
        
    Examples
    --------
    inst_ECEF = gei2ecef(inst_GEI, time_array)

    """
    SizeArr = np.size(time_array)

     # Initialize GEI/ECI Vectors
    inst_GEI2ECEF = np.empty(inst_GEI.shape)
    inst_GEI2ECEF2 = np.empty(inst_GEI.shape)
    theta_rad = np.empty(SizeArr)

    for i in range(0, SizeArr):
        theta_rad[i] = GMST_midnight(time_array[i])
        
        inst_GEI2ECEF[:, i] = gei_to_ecef_rm(theta_rad[i]) @ inst_GEI[:, i]

        denominator = np.sqrt(inst_GEI2ECEF[0, i]**2 + \
                              inst_GEI2ECEF[1, i]**2 + \
                              inst_GEI2ECEF[2, i]**2)

        inst_GEI2ECEF2[:, i] = inst_GEI2ECEF[:, i]/denominator


    inst_ECEF = np.array(inst_GEI2ECEF2)

    return inst_ECEF

# =============================================================================
# # %% function for ECEF to NED: (North, East, Down) transformations
# # z is positive downwards
# =============================================================================
def terrestrial2ned(inst_ter, lat, lon):
    """
    function for transformations from ECEF to NED (z is positive downwards)

    Parameters
    ----------
    inst_ter : np.ndarray[float]
        instrument look direction in terrestrial frame.
    lat : numpy.ndarray[float]
        Geodetic latitude in degrees.
    lon : numpy.ndarray[float]
        Geodetic longitude in degrees.

    Returns
    -------
    inst_NED : np.ndarray[float]
        instrument look direction in NED.
        
    Examples
    --------
    inst_NED = terrestrial2ned(inst_ter, lat, lon)

    """
   
    # Initialize NED Vectors
    inst_ter2NED = np.empty(inst_ter.shape)
    inst_NED=np.empty(inst_ter.shape)

    for i in range(0, inst_ter.shape[1]):

        inst_ter2NED[:, i] = terrestrial_to_ned_rm(lat[i], \
                                                    lon[i]) @ inst_ter[:, i]

        denominator = np.sqrt(inst_ter2NED[0, i]**2 + \
                              inst_ter2NED[1, i]**2 + \
                              inst_ter2NED[2, i]**2)

        inst_NED[:, i]= inst_ter2NED[:, i] / denominator

    return inst_NED.T

# =============================================================================
# # %% function for NED to ENU transformations
# =============================================================================
def ned2enu(inst_NED):
    """
    Function for NED to ENU transformations
    
    The transformation matrice from NED to ENU and ENU to NED are the same.
    [0 1 0 ; 1 0 0 , 0 0 -1]
    in ENU, z increases in upward direction.

    Ref: Grewal, M. S., Weill, L. R., & Andrews, A. P. (2008).
    Appendix C: Coordinate Transformations. In Global Positioning Systems,
    Inertial Navigation, and Integration (pp. 456–501). John Wiley & Sons, Inc.
    https://doi.org/10.1002/9780470099728.app3

    Parameters
    ----------
    inst_NED : np.ndarray[float]
        instrument look direction in NED.

    Returns
    -------
    inst_ENU : np.ndarray[float]
        instrument look direction in ENU.

    """
    # # Initialize ECEF Vectors
    SC_NED = np.array(inst_NED).T

    # Initialize ENU Vectors
    inst_ENU = np.empty(SC_NED.shape)

    # rotation matrice from NED to ENU
    RN2E= np.array([[0, 1, 0],
                    [1, 0, 0],
                    [0, 0, -1]])

    inst_ENU = RN2E @ SC_NED

    return inst_ENU


def orf_to_j2k_use_spacecraft_ephemeris(body_vec, pX, pY, pZ, Vx, Vy, Vz):
    """
    Employs spacecraft ephemeris to transform from ORF to GEIJ2K: 
    
    Parameters
    ----------
    body_vec : np.ndarray[float]
        rotated body frame vector.
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
    P  : numpy.ndarray
        optional reordering matrix to align the perifocal frame with the 
    orbital frame. Default is P= np.array([[0, 0, -1],[1, 0, 0],[0,-1,0]])

    Returns
    -------
    inst_GEI : numpy.ndarray[float]
        instrument look direction in GEIJ2K.
        
    Examples
    --------
    inst_GEI = \
        
        orf_to_j2k_use_spacecraft_ephemeris(body_vec, pX, pY, pZ, Vx, Vy, Vz)
    """

    # Initialize GEI/ECI Vectors
    SC_orf_J2K = np.empty(body_vec.shape)
    SC_orf_J2K2 = np.empty(body_vec.shape)

    for i in range(0, body_vec.shape[1]):
        R = orf_to_j2k_rm(pX[i], pY[i], pZ[i], Vx[i], Vy[i], Vz[i])
        
        SC_orf_J2K[:, i] = R @ (body_vec[:, i])
        #  normalize the matrix
        denominator= np.sqrt(SC_orf_J2K[0, i]**2 + SC_orf_J2K[1, i]**2 + \
                             SC_orf_J2K[2, i]**2)

        SC_orf_J2K2[:, i]= SC_orf_J2K[:, i]/denominator

    inst_GEI = np.array(SC_orf_J2K2)

    return inst_GEI

def icrf2itrf(inst_GEI, path_to_files, time_array):
    """
    Function to transform the pointing vector in GEIJ2K to ITRF.
    
    Notes
    -----
    Uses Python wrappers for the Standards of Fundamental Astronomy library
    https://www.iausofa.org/

    Parameters
    ----------
    inst_GEI : numpy.ndarray[float]
        instrument look direction in GEIJ2K.
    path_to_files : str
        IES and EOP file path.
    time_array : datetime.datetime
        time array.

    Returns
    -------
    inst_ITRF : numpy.ndarray[float]
        instrument look direction in ITRF.
    
    Examples
    --------
    inst_ITRF = icrf2itrf(inst_GEI, path_to_files, time_array, 
                          GEOx, GEOy, GEOz)

    """
    row, col = np.shape(inst_GEI)
    inst_ITRF_Vector=np.empty((col,3)); inst_ITRF_Vector.fill(np.nan)

    for i in range(0, col):
        c2t = icrf_to_itrf_rm(path_to_files, time_array[i])
        inst_ITRF_Vector[i, :] = c2t @ inst_GEI[:,i]

    inst_ITRF = inst_ITRF_Vector.T
    
    return inst_ITRF

def terrestrial2nec(inst_ter, gX, gY, gZ):
    """
    Function to transform instrument look direction from terrestrial frame
    to NEC. Uses (Nielsen, 2019) method.
    
    Parameters
    ----------
    inst_ter : numpy.ndarray[float]
        Instrument look direction in terrestrial frame.
    gX : numpy.ndarray[float]
        X position in GEO/ITRF (km).
    gY : numpy.ndarray[float]
        Y position in GEO/ITRF (km).
    gZ : numpy.ndarray[float]
        Z position in GEO/ITRF (km).

    Returns
    -------
    inst_NEC : numpy.ndarray[float]
        instrument look direction in NEC.
        
    Examples
    --------
    inst_NEC = itrf2nec(inst_ITRF, GEOx, GEOy, GEOz)

    """
    row, col = np.shape(inst_ter)
    # transform inst_ITRF to inst_NEC using quaternions
    inst_NEC=np.empty((col,3)); inst_NEC.fill(np.nan)
    
    for i in range(0, col):

        rNEC_from_ITRF = terrestrial_to_nec_rm(gX[i], gY[i], gZ[i])

        inst_NEC[i, :]  = np.dot( rNEC_from_ITRF,  inst_ter[: , i])

    return inst_NEC

def nec2terrestrial(inst_NEC, gX, gY, gZ):
    """
    Function to transform instrument look direction from NEC to terrestrial 
    frame. Uses the transpose of ITRF to NEC rotation matrix.

    Parameters
    ----------
    inst_NEC : numpy.ndarray[float]
        Instrument look direction in NEC.
    gX : numpy.ndarray[float]
        X position in GEO/ITRF (km).
    gY : numpy.ndarray[float]
        Y position in GEO/ITRF (km).
    gZ : numpy.ndarray[float]
        Z position in GEO/ITRF (km).

    Returns
    -------
    inst_ter : numpy.ndarray[float]
        Instrument look direction in ITRF.
        
    Examples
    --------
    inst_ITRF = nec2itrf(inst_NEC, GEOx, GEOy, GEOz)

    """

    row, col = np.shape(inst_NEC)
    # transform inst_ITRF to inst_NEC using quaternions
    inst_ter=np.empty((col,3)); inst_ter.fill(np.nan)
    
    inst_nec_vec = inst_NEC.T

    for i in range(0, col):

        rITRF_from_NEC = nec_to_terrestrial_rm(gX[i], gY[i], gZ[i])

        inst_ter[i, :]  = np.dot( rITRF_from_NEC,  inst_nec_vec[: , i])

    return inst_ter

