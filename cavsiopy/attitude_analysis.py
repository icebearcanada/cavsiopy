#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
attitude_analysis module includes functions for finding the instrument
pointing direction in spacecraft orbital frame(ORF), GEI J2000, ECEF, ITRF,
 NEC, NED and ENU.
Additional functions are available to calculate the look angles of the
spacecraft and the instrument with respect to a ground point, line-of-sight
look direction vector from the ground point to the spacecraft and the distance
between spacecraft and the ground point. Slew quality can be calculated given
the threshold for the angle between the LOS direction and the
instrument pointing.

.. toctree::
  :maxdepth: 2   
  LA_inst
  LA_sat
  calculate_los_vec
  calculate_reception_angle
  find_instrument_attitude
  find_slew_inst
  find_slew_rri
  rotate_inst
  spacecraft_distance_from_a_point

@author: ceren
"""
import datetime
import numpy as np
import math as m
from geopy.distance import geodesic
import cavsiopy.ephemeris_importer as ei
import cavsiopy.use_rotation_matrices as urm

def LA_sat(plon, plat, palt, slon, slat, salt):
    """
    Calculates the satellite look angle: elevation and azimuth angle that is 
    needed for the satellite to look towards a point in GEO frame.

    Ref: ASD5: Earth Station Look Angles By Prof. Gregory D. Durgin,
                Astrodynamics, Georgia Tech, 2009.  

    Parameters
    ----------
    plon : float
        Longitude of the target on the ground (degrees).
    plat : float
        Latitude of the target on the ground (degrees).
    palt : float
        altitude of the point in GEO frame (km).
    slon : numpy.ndarray[float]
        spacecraft longitude (degrees).
    slat : numpy.ndarray[float]
        spacecraft latitude (degrees).
    salt : numpy.ndarray[float]
        altitude of satellite in GEO frame (km).

    Returns
    -------
    dict : dict
        Look angles of the spacecraft: azimuth and elevation.
        azimuth: numpy.ndarray[float]
            Satellite location in terms of azimuth (degrees).
        elevation: numpy.ndarray[float]
            Elevation angle of the spacecraft (degrees).
            
    Examples
    --------
    sat_look_angles = LA_sat(plon, plat,  palt, slon, slat, salt) 
    """
 
    rp= 6371 + palt
    la_elev=np.empty(len(slat)) # elevation angle
    la_azim=np.empty(len(slon)) # azimuthal angle

    for i in range(0, len(slat)):
        rslat = slat[i]*(m.pi/180) #convert to radians
        rslon = slon[i]*(m.pi/180) #convert to radians
        rplat = plat*(m.pi/180) #convert to radians
        rplon = plon*(m.pi/180) #convert to radians
        rs = salt[i] + 6371

        a = m.sin(rslat) * m.sin(rplat)
        b = m.cos(rplat) * m.cos(rslat) * m.cos (rslon-rplon)

        gamma = m.acos(a+b)

        c = m.sin(gamma)
        d = m.sqrt(1 + ( (rp/rs)**2 ) - 2*( (rp/rs) * m.cos(gamma)))
        la_elev[i] = m.acos(c/d) * (180 / m.pi)


        e = m.sin(abs(rplon-rslon)) * (m.cos(rslat)/c)

        if (rslat<rplat) and (rslon<rplon):
            # when passage is SW of Earth station
            la_azim[i] = (m.asin(e) + m.pi) * (180 / m.pi)
        elif (rslat<rplat) and (rslon>rplon):
            # when passage is SE of Earth station
            la_azim[i] = (m.pi - m.asin(e)) * (180 / m.pi)
        elif (rslat>rplat) and (rslon<rplon):
            # when passage is NW of Earth station
            la_azim[i] = (2*m.pi - m.asin(e)) * (180 / m.pi)
        else:
            la_azim[i] = m.asin(e) * (180 / m.pi)

    return {'azimuth': la_azim, 'elevation': la_elev}

def LA_inst(vec, e, n, lon, lat, pLon):
    """
    Parameters
    ----------
    vec : Array of float
        Boresight pointing vector. 
    e : int
        Column of the east vector in vec. 
    n : int
        Column of the north vector in vec.
    lon : numpy.ndarray[float]
        spacecraft longitude (degrees).
    lat : numpy.ndarray[float]
        spacecraft latitude (degrees).
    pLon : float
        Longitude of the target on the ground (degrees).

           
    Returns
    -------
    dict : dict
        Look angles of the instrument: azimuth and depression. 
        vec_azim: Array of float
            Look angle of the instrument in horizontal plane (degrees)
        vec_dep: Array of float
            Depression angle of the instrument (degrees)
    
    Examples
    --------
    inst_look_angles = LA_rri(RRI_enu, 0, 1, Lon, Lat) 

    """
    
    size = len(lat)
    vec_dep=np.empty(size); vec_dep1=np.empty(size);
    vec_azim=np.empty(size); vec_azim1=np.empty(size);
    
    #  Define the horizontal plane
    if e == 0:
        # vertical is positive upward in ENU, ours looks down
        horizontal_plane = [0, 0, -1]
    elif e == 1: 
        # vertical is positive downward in NED, ours looks down
        horizontal_plane = [0, 0, 1]
    
    # initialize X-Z look vector of the instrument
    XZ = np.empty(size); XZ.fill(np.nan)
    
    for i in range(0,size):
        #  calculate X-Z look vector of the instrument for azimuth calculations
        XZ[i] = np.sqrt(vec[i,e]**2 + vec[i,2]**2)
    
    #  calculate azimuth
    for i in range(0,size):
        if np.isnan(XZ[i]) == False:
            vec_azim1[i] = np.arctan2(vec[i][e], vec[i][n])*180.0/np.pi  
            
            if vec_azim1[i] <= 0:
                vec_azim1[i] = vec_azim1[i] +360
    
            if lon[i] < pLon: # west of ottawa
                if 0 < vec_azim1[i] <= 180: # north of Ottawa
                    vec_azim[i] = vec_azim1[i]+180
                else:
                    vec_azim[i] = vec_azim1[i]+360
            
            else: # east of Ottawa
                if 0 < vec_azim1[i] <= 180:
                    vec_azim[i] = 180 - vec_azim1[i]
                    # vec_azim[i] = vec_azim1[i] # new
                else:
                    vec_azim[i] = 180 + vec_azim1[i]
           
            horizontal_plane_dot_RRI = np.dot(vec[i][:], horizontal_plane)
            vec_dep1 = np.arccos(horizontal_plane_dot_RRI)*180.0/np.pi
            vec_dep[i] = 90-vec_dep1 # instrument depression angle      
        else:
            vec_azim[i] = np.nan 
            vec_dep[i] = np.nan
            
    return {'azimuth': vec_azim, 'depression': vec_dep}

def spacecraft_distance_from_a_point(pLon, pLat, pAlt, slon, slat, salt):
    """
    calculates spacecraft distance from a point using the 
    geopy.dist.geodesic

    Parameters
    ----------
    plon : float
        Longitude of the target on the ground (degrees).
    plat : float
        Latitude of the target on the ground (degrees).
    palt : float
        altitude of the point in GEO frame (km).
    slon : numpy.ndarray[float]
        spacecraft longitude (degrees).
    slat : numpy.ndarray[float]
        spacecraft latitude (degrees).
    salt : numpy.ndarray[float]
        altitude of satellite in GEO frame (km).

    Returns
    -------
    distance : np.ndarray[float]
        slant range between the ground point and spacecraft (km).

    """

    SizeArr = len(salt)

    for i in range(0, SizeArr):
        # Line-of-sight distance between the point and satellite
        distance = np.empty(salt.shape);
        # positive down vector (D in NED)
        dz_down= pAlt-salt[i]

        # find the changes in latitude by keeping longitude constant
        point1 = (pLat, pLon)
        sat1 = (slat[i], pLon)
        # Distance due to changes in latitude (y axis in this case)
        dy = geodesic(point1, sat1).kilometers

        # find the changes in latitude by keeping spacecraft and ground point
        # latitude the same
        point2 = (slat[i], pLon)
        sat2 = (slat[i], slon[i])
        # Distance due to changes in longitude (x axis in this case)
        dx = geodesic(point2, sat2).kilometers
        distance[i] = np.sqrt(dx**2 + dy**2 + dz_down**2)

        return distance
    
def calculate_los_vec(pLon, pLat, pAlt, slon, slat, salt):
    """
    calculates line-of-sight vector from a point to the spacecraft location
    using the geopy.dist.geodesic

    Parameters
    ----------
    plon : float
        Longitude of the target on the ground (degrees).
    plat : float
        Latitude of the target on the ground (degrees).
    palt : float
        altitude of the point in GEO frame (km).
    slon : numpy.ndarray[float]
        spacecraft longitude (degrees).
    slat : numpy.ndarray[float]
        spacecraft latitude (degrees).
    salt : numpy.ndarray[float]
        altitude of satellite in GEO frame (km).

    Returns
    -------
    los_enu_arr : numpy.ndarray[float]
        line-of-sight vector from the point to the spacecraft in ENU system.
    los_ned_arr : numpy.ndarray[float]
        line-of-sight vector from the point to the spacecraft in NED system.

    """       
    SizeArr = len(salt)

    # Line-of-sight distance between the point and satellite
    distance = np.empty(salt.shape);
    # Ray vector
    los_ned = []; los_enu = [];

    for i in range(0, SizeArr):

    # ========================================================================
    # calculate the ray vector
    # ========================================================================
        # positive down vector (D in NED)
        dz_down = pAlt-salt[i]
        # positive up vector (E in ENU)
        dz_up = salt[i]-pAlt

        # find the changes in latitude by keeping longitude constant
        point1 = (pLat, pLon)
        sat1 = (slat[i], pLon)
        # Distance due to changes in latitude (y axis in this case)
        dy = geodesic(point1, sat1).kilometers

        # find the changes in latitude by keeping spacecraft and ground point
        # latitude the same
        point2 = (slat[i], pLon)
        sat2 = (slat[i], slon[i])
        # Distance due to changes in longitude (x axis in this case)
        dx = geodesic(point2, sat2).kilometers
        distance[i] = np.sqrt(dx**2 + dy**2 + dz_down**2)

        # Y is positive northward, X is positive eastward. For other directions 
        # change of sign is needed:

        # if spacecraft is in the south of Ottawa
        if slat[i] < pLat:
            dy = -dy
        # if spacecraft is in the west of Ottawa
        if 360 + slon[i] < 360 + pLon:
            dx = -dx

        # NED: down is positive
        los_ray_ned = [dy/distance[i], dx/distance[i], dz_down/distance[i]]
        # ENU: up is positive
        los_ray_enu = [dx/distance[i], dy/distance[i], dz_up/distance[i]]
        # append the ray vector values
        los_ned.append(los_ray_ned)
        los_enu.append(los_ray_enu)

    los_enu_arr=np.array(los_enu)
    los_ned_arr=np.array(los_ned)

    return los_enu_arr, los_ned_arr   

def rotate_inst(body_vec, Roll, Pitch, Yaw):
    """
    Rotation of body vectors in orbital frame (ORF) of Swarm-E
    using roll, pitch, and yaw angles.
    
    The rotation sequence is X->Y->Z; roll->pitch->yaw according to the RRI
    data description.
    
    Ref: University of Calgary, e-POP Radio Receiver Instrument (RRI)
    User’s Manual, Doc. no. ePOP-5024, Rev. D (2018)
    
    https://epop.phys.ucalgary.ca/data-handbook/coordinate-systems/
    
    Parameters
    ----------
    body_vec : tuple
        x, y, z : initial position vectors for antenna in ORF
    Roll : float
        roll angle in degrees
    Pitch : float
        pitch angle in degrees
    Yaw : float
        yaw angle in degrees
    
    Returns
    -------
    rotated_body : numpy.ndarray
        rotated body vectors of instruments onboard Swarm-E in ORF
    """
    # converting degrees to radians
    rr = np.deg2rad(Roll)
    pr = np.deg2rad(Pitch)
    yr = np.deg2rad(Yaw)
    
    # ORF (Local vertical local horizontal) 3 axis unit vectors
    orf = np.array([body_vec,] *len(Roll)).T
    
    # Initialize body Vectors
    rot_vec = np.empty(orf.shape)
    
    # transform initial definition of LVLH to rotated frame
    for i in range(len(rr)):
    
        # direction cosine matrice
        R = urm.RX_r2i(rr[i]) @ urm.RY_r2i(pr[i]) @ urm.RZ_r2i(yr[i])
    
        # rotation
        rot_vec[:, i] = R @ orf[:, i]
    
    rotated_body = np.array(rot_vec)
    
    return rotated_body
    
def find_instrument_attitude(rotated_body, geiX, geiY, geiZ, 
                             geiVx, geiVy, geiVz, 
                             geoX, geoY, geoZ,
                             time_array, start_date, end_date, 
                             lat, lon, path_to_sofa_files, method1='ephemeris', 
                             frame2 = 'itrf', frame3 = 'nec' ):
    """
    Takes in the rotated body vector in orbital frame and outputs instrument
    look direction in NEC or NED coordinate systems by utilizing the functions
    in use_rotation_matrices module.
    
    Parameters
    ----------
    rotated_body : TYPE
        DESCRIPTION.
    geiX : numpy.ndarray[float]
        X position in GEIJ2K (km).
    geiY : numpy.ndarray[float]
        Y position in GEIJ2K (km).
    geiZ : numpy.ndarray[float]
        Z position in GEIJ2K (km).
    geiVx : numpy.ndarray[float]
        X component of velocity in GEIJ2K (km/s).
    geiVy : numpy.ndarray[float]
        Y component of velocity in GEIJ2K (km/s).
    geiVz : numpy.ndarray[float]
        Z component of velocity in GEIJ2K (km/s).
    geoX : numpy.ndarray[float]
        X position in GEO (km).
    geoY : numpy.ndarray[float]
        Y position in GEO (km).
    geoZ : numpy.ndarray[float]
        Z position in GEO (km).
    time_array : datetime.datetime
        time.
    start_date : datetime.datetime
        start of the experiment.
    lat : numpy.ndarray[float]
        Geodetic latitude in degrees.
    lon : numpy.ndarray[float]
        Geodetic longitude in degrees.
    path_to_sofa_files : str
        path_to_initialization files (IERS and EOP).
    method1 : str, optional
        Transformation method to ICRF/GEI J2K.
        Can be 'ephemeris' or 'orbital_elements'. The default is 'ephemeris'.
    frame2 : str, optional
        Terrestrial frame: 'ecef' or 'itrf'. The default is 'itrf'.
    frame3 : str, optional
        Final coordinate system: 'nec' or 'ned'. The default is 'nec'.

    Returns
    -------
    inst_geo : numpy.ndarray[float]
        instrument look direction in the requested coordinate system:NEC or NED

    """
    
    #  first frame to transform is the icrf
    # if spacecraft ephemeris method is chosen to transform to icrf/gei j2k
    if method1 == 'ephemeris':
        inst_GEI = \
            urm.orf_to_j2k_use_spacecraft_ephemeris(rotated_body, 
                                                    geiX, geiY, geiZ, 
                                                    geiVx, geiVy, geiVz)
        # if itrf is chosen as the second frame
        if frame2 == 'itrf':
            inst_ITRF = \
                urm.icrf2itrf(inst_GEI, path_to_sofa_files, time_array)
            # if nec is chosen as the third frame
            if frame3 == 'nec':
                inst_geo = urm.terrestrial2nec(inst_ITRF, geoX, geoY, geoZ)
                
            # if ned is chosen as the third frame
            else:
                inst_geo = urm.terrestrial2ned(inst_ITRF, lat, lon)
                
        # if ecef is chosen as the second frame
        else:
            inst_ECEF = urm.gei2ecef(inst_GEI, time_array)
            # if nec is chosen as the third frame
            if frame3 == 'nec':
                inst_geo = urm.terrestrial2nec(inst_ECEF, geoX, geoY, geoZ)
            # if ned is chosen as the third frame
            else:
                inst_geo = urm.terrestrial2ned(inst_ECEF, lat, lon)           
    # if spacecraft ephemeris is chosen to transform to icrf/gei j2k    
    elif method1 == 'orbital_elements':
        inst_GEI = \
            urm.orf_to_j2k_use_orbital_elements(rotated_body, 
                                                geiX, geiY, geiZ, 
                                                geiVx, geiVy, geiVz)
        # if itrf is chosen as the second frame
        if frame2 == 'itrf':
            inst_ITRF = \
                urm.icrf2itrf(inst_GEI, path_to_sofa_files, time_array)
            # if nec is chosen as the third frame
            if frame3 == 'nec':
                inst_geo = urm.terrestrial2nec(inst_ITRF, geoX, geoY, geoZ)
                
            # if ned is chosen as the third frame
            else:
                inst_geo = urm.terrestrial2ned(inst_ITRF, lat, lon)
                
        # if ecef is chosen as the second frame
        else:
            inst_ECEF = urm.gei2ecef(inst_GEI, time_array)
            # if nec is chosen as the third frame
            if frame3 == 'nec':
                inst_geo = urm.terrestrial2nec(inst_ECEF, geoX, geoY, geoZ)
            # if ned is chosen as the third frame
            else:
                inst_geo = urm.terrestrial2ned(inst_ECEF, lat, lon)                
    
    return inst_geo
    
def calculate_reception_angle(inst_ned, pLat, pLon, pAlt, lat, lon, alt, 
                              inst = 'boresight'):
    """
    function to calculate the reception angle of an instrument.
    reception angle: angle between the instrument look direction vector and 
    the line-of-sight vector from the target
    
    Parameters
    ----------
    inst_ned : numpy.ndarray[float]
        instrument look direction vector in North-East-Down.
    pLat : float
        geodetic latitude of the target (degrees).
    pLon : float
        geodetic longitude of the target (degrees).
    pAlt : float
        Altitude of the target (km).
    lat : numpy.ndarray[float]
        spacecraft position in geodetic latitude (degrees).
    lon : numpy.ndarray[float]
        spacecraft position in geodetic latitude(degrees).
    alt : numpy.ndarray[float]
        spacecraft altitude (km).
    inst : str, optional
        cra is calculated as (90-ra) for dipoles, (180-ra) for boresight. 
        input can be boresight or dipole. The default is 'boresight'.

    Returns
    -------
    ra_los : numpy.ndarray(float)
        reception angle of the instrument (degrees).
    cra_los : numpy.ndarray(float)
        complementary reception angle of the instrument (degrees).

    """
    
    ra_los=np.empty(alt.shape); ra_los.fill(np.nan)
    cra_los=np.empty(alt.shape); cra_los.fill(np.nan)

    SizeArr = len(alt)

    # Line-of-sight distance between the point and satellite
    distance=np.empty(alt.shape)

    for i in range(0, SizeArr):
        # =====================================================================
        # calculate the line-of-sight ray vector
        # =====================================================================
        # positive down vector (D in NED)
        dz_down= pAlt-alt[i]
    
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
        distance[i] = np.sqrt(dx**2 + dy**2 + dz_down**2)  
        
        # Y is positive northward, X is positive eastward.In other cases change
        # of sign is needed:
        # if spacecraft is in the south of Ottawa
        if lat[i] < pLat:
            dy = -dy
        # if spacecraft is in the west of Ottawa
        if 360+lon[i] < 360+pLon:
            dx = -dx   
       
        # NED: down is positive
        los_ray_ned = [dy/distance[i], dx/distance[i], dz_down/distance[i]]
        
        # Find the angle between the LOS ray vector and RRI pointing direction
        dot_productXR_los = np.dot(inst_ned[i][:], los_ray_ned)
        ra_los[i] = np.arccos(dot_productXR_los)*180.0/np.pi
        
        if inst == 'boresight':
            cra_los[i] = 180 - ra_los[i]
        else:
            cra_los[i] = 90 - ra_los[i]
        
    return ra_los, cra_los


def find_slew_rri(ra_los, ra_D1, ra_D2, criteria):
    """
    function to classify slew according to the given criteria
        2: Front face slew 
        1: 1 dipole back, 1 dipole front, boresight front slew
        0.5 : 1 dipole front slew
        0: no dipoles slewed
        -0.5: 1 dipole back slew
        -1: 1 dipole back, 1 dipole front, boresight back slew
        -2: Back face slew

    Parameters
    ----------
    ra_los : numpy.ndarray[float]
        Boresight reception angle (degrees).
    ra_D1 : numpy.ndarray[float]
        Dipole1 reception angle (degrees).
    ra_D2 : numpy.ndarray[float]
        Dipole2 reception angle (degrees).
    criteria : numpy.ndarray[float]
        criteria/threshold for slew (degrees).

    Returns
    -------
    slew : int
        Classification of slew according to the given criteria

    """ 
    cra_los = 180 - ra_los
    cra_D1 = 90-ra_D1
    cra_D2 = 90-ra_D2

    slew = np.empty(len(ra_los)); slew.fill(np.nan)
    
    # slew conditions
    for i in range(0, len(cra_los)):
        cond2 = 0 <= cra_los[i] <= criteria # front slew
        cond2_2 = \
            abs(cra_D1[i]) < criteria and abs(cra_D2[i]) < criteria

        ## back slew
        # cond2_2 also holds true here
        cond_neg2 = -criteria <= cra_los[i] <= -0.001 

        ## one dipole front slew (slew = 1)===================================
        cond1 = \
            abs(cra_los[i]) > criteria or abs(cra_los[i]) <= criteria
        # dipole 1 front slew
        cond1_1 = \
            abs(cra_D2[i]) > criteria and (0 <= cra_D1[i] < criteria)
        # dipole 2 front slew
        cond1_2 = \
            (0 <= cra_D2[i] < criteria) and abs(cra_D1[i]) > criteria

        ##  one dipole back slew (slew = -1)==================================
        # dipole 1 back slew
        cond1_3 = \
            abs(cra_D2[i]) > criteria and \
                (-criteria < cra_D1[i] < 0)
        # dipole 2 back slew
        cond1_4 = \
            (-criteria < cra_D2[i] < 0) and \
                abs(cra_D1[i]) > criteria

        ## one dipole front slew with data containing nans (slew = 1)=========
        # dipole 1 front slew
        cond1_5 = \
            np.isnan(cra_D2[i])==True and (0 <= cra_D1[i]< criteria)

        # dipole2 2 front slew
        cond1_6 = \
            (0 <= cra_D2[i] < criteria) and np.isnan(cra_D1[i]) ==True

        ## one dipole back slew with data containing nans (slew = -1)---------
        # dipole 1 back slew
        cond1_7 = np.isnan(cra_D2[i])==True and (-criteria < cra_D1[i] < 0)
        # dipole 2 back slew
        cond1_8 = \
            (-criteria < cra_D2[i] < 0) and np.isnan(cra_D1[i]) == True
            
         #  one dipole front, 1 dipole back slew: slew
        cond1_9 = (cond1_7 ==True and cond1_6 == True)
        cond1_10 = (cond1_8 ==True and cond1_5 == True)

        ## no slew
        cond0_1 = cra_los[i] > criteria or cra_los[i] < -criteria
        cond0_2 = \
            abs(cra_D1[i]) > criteria and abs(cra_D2[i]) > criteria
            
        # complete front slew
        if cond2 ==True and cond2_2 ==True:
            slew[i] = 2
        #  one dipole front, 1 dipole back slew, los = front : slew
        elif cond2 == True and (cond1_9==True or cond1_10 == True):
            slew[i] = 1
        # one dipole front slew
        elif cond1 == True and (cond1_1==True or cond1_2 == True):
            slew[i] = 0.5
        # one dipole front slew
        elif cond1 == True and (cond1_5==True or cond1_6 == True):
            slew[i] = 0.5
        # one dipole back slew
        elif cond1 == True and (cond1_7==True or cond1_8 == True):
            slew[i] = -0.5
        # one dipole back slew
        elif cond1 == True and (cond1_3==True or cond1_4 == True):
            slew[i] = -0.5
        #  one dipole front, 1 dipole back slew, los = back : slew
        elif cond_neg2 == True and (cond1_9==True or cond1_10 == True):
            slew[i] = -1
        # complete back slew
        elif cond_neg2 == True and cond2_2==True:
            slew[i] = -2 # back slew
        elif cond0_1==True and cond0_2== True:
            slew[i]= 0 # no slew
        else:
            slew[i]= 0 # no slew

    return slew

def find_slew_inst(ra_los, criteria):
    """
    function to classify slew according to the given criteria
        1 : Front face slew
        0: no dipoles slewed
        -1: Back face slew

    Parameters
    ----------
    ra_los : numpy.ndarray[float]
        Boresight reception angle (degrees).
    criteria : float
        criteria/threshold for slew (degrees).

    Returns
    -------
    slew : int
        Classification of slew according to the given criteria

    """ 
    
    # slew conditions
    if type(ra_los) == float or type(ra_los) == int:
        
        cra_los = 180 - ra_los
        cond_front = 0 <= cra_los <= criteria # front face slew
    
        ## back face slew
        cond_back = -criteria <= cra_los <= -0.0001 
    
        ## no slew
        cond0 = cra_los > criteria or cra_los < -criteria
    
        # complete front slew
        if cond_front == True:
            slew = 1
        # complete back slew
        elif cond_back == True:
            slew = -1 
        # no slew
        elif cond0 == True:
            slew = 0 
        else:
            slew = np.nan # unclassified
    else:
        slew = np.empty(len(ra_los)); slew.fill(np.nan)
        cra_los= 180 - np.array(ra_los)
        
        for i in range(0, len(cra_los)):
            
            cond_front = 0 <= cra_los[i] <= criteria # front face slew
    
            ## back face slew
            cond_back = -criteria <= cra_los[i] <= -0.0001 
    
            ## no slew
            cond0 = cra_los[i] > criteria or cra_los[i] < -criteria
    
            # complete front slew
            if cond_front == True:
                slew[i] = 1
            # complete back slew
            elif cond_back == True:
                slew[i] = -1 # back slew
            elif cond0 == True:
                slew[i]= 0 # no slew
            else:
                slew[i]= np.nan # no slew

    return slew

