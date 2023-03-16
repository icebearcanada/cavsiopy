#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 14:32:10 2021

1. Calculates satellite orbital elements using X, Y, Z coordinates in GEI
2. Imports satellite parameters from Celestrak TLE file
3. Compares calculated orbital elements with the TLE file
4. Handles Euler rotations in LVLH frame.
5. Includes rotation matrices for coordinate transformations.

@author: ceren
"""

import numpy as np

import cavsiopy.ephemeris_importer as ei

import datetime

import math as m

from geopy.distance import geodesic

import cavsiopy.use_rotation_matrices as RM

def inst_direction_in_ned_using_icrf_to_ecef(body_inst, Roll, Pitch, Yaw,
                         GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz,
                         time_array, lat, lon, method_name):
        """
        Created on Thu Jul  8 13:22:03 2021

        This code calculates the direction of pointing vectors in ENU
            when given the body vectors of monopoles or dipoles in LVLH.

        Parameters
        ----------
            body_inst: [x, y, z] location of instrument wrt SC frame
            Roll, Pitch, Yaw: Euler angles
            lvlh_x, lvlh_y, lvlh_z: antenna x, y, z vector components in LVLH
            GEIx, GEIy, GEIz: Position of the satellite in GEI
            GEIVx, GEIVy, GEIVz: Spacecraft Velocity Components in GEI
            time_array: datetime
            lat, lon: spacecraft location in GEO_spherical
            method_name: 
                'oe' for rotations from SC frame to GEI/ECI using orbital
                elements calculated from GEI r and V
                'other' for rotations from SC frame to GEI/ECI 
                using r and V in GEI by means ephemeris to Nadir rotation matrice

        Returns
        -------
            ENUx, ENUy, ENUz: pointing direction of 1 monopole in ENU coordinates

        @author: ceren
        """
        # values from the database
        GEI = [GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz]
        # 1- Rotation of SC vectors in LVLH frame
        rotated_inst = RM.rotate_antenna(Roll, Pitch, Yaw, body_inst)

        if method_name == 'oe':
            # 2- rotation of SC frame to GEI/ECI frame using orbital elements in GEI
            inst_GEI= RM.ORF2J2K_use_orbital_elements(rotated_inst, GEI)
            # 3- rotation of GEI/ECI frame to Earth Centered Earth Fixed frame (ECEF):
            inst_ECEF = RM.GEI2ECEF(time_array, lon, inst_GEI)
            NED_vec = RM.ECEF2NED(lat, lon, inst_ECEF)
        else:
            # 2- rotation of SC frame to GEI/ECI frame using r and V in GEI
            inst_GEI = RM.ORF2J2K_use_spacecraft_ephemeris(rotated_inst, \
                                        GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz)
            # 3- rotation of GEI/ECI frame to Earth Centered Earth Fixed frame (ECEF):
            inst_ECEF = RM.GEI2ECEF(time_array, lon, inst_GEI)
            NED_vec = RM.ECEF2NED(lat, lon, inst_ECEF)

        # 4- rotation of ECEF frame to NED:
        inst_ned = NED_vec.T
        ENU = RM.NED2ENU(inst_ned)
        inst_enu = ENU.T

        return inst_ned, inst_enu
# =============================================================================
# function for finding the instrument look direction using ICRF to ITRF
# =============================================================================
def inst_direction_in_nec_using_icrf_to_itrf(path_to_files, \
                                   filedate, time_array, start_date, end_date,\
                                GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz, \
                                MET, GEO, roll, pitch, yaw):

    # specify the location of the instrument wrt [x, y, z].
    a=1
    body_RRI = [a, 0, 0]

    # rotate the antenna
    rri_lvlh = RM.rotate_antenna(roll, pitch, yaw, body_RRI)
    # find the instrument look direction in GEI
    icrf_RRI = RM.ORF2J2K_use_spacecraft_ephemeris(rri_lvlh, GEIx, GEIy, GEIz,\
                                                    GEIVx, GEIVy, GEIVz)

    # transform from ICRF to ITRF
    RRI_ITRF, RRI_nec = RM.icrf2itrf(path_to_files, filedate, MET, GEO, \
                                     time_array, start_date, end_date, icrf_RRI)
    # transform nec to enc
    RRI_enc = RM.NED2ENU(RRI_nec).T

    return RRI_ITRF, RRI_enc, RRI_nec

def inst_direction_in_nec_using_oe_icrf_to_itrf(path_to_files,  \
                                  file_sat, filedate, time_array,\
                            start_date, end_date, Roll, Pitch, Yaw, body_inst):
    """ This function uses orbital elements to rotate the instrument direction
    in LVLH frame
    Parameters
    ----------
    path_to_files: str
    file_sat: file that contains the ephemeris information
    filedate: datetime
    timearray: datetime
    start_date: datetime
    Roll, Pitch, Yaw: arrays of float
    body_inst: instrument body vector

    Returns
    -------
    inst_NEC, inst_ENC: float
        instrument body vector in NEC and ENC

    """
    # 1- Rotation of SC vectors in LVLH frame
    rotated_inst = RM.rotate_antenna(Roll, Pitch, Yaw, body_inst)

    inst_GEI= RM.ORF2J2K_use_orbital_elements(rotated_inst, file_sat)
    # handle the transformation
    RRI_ITRF, RRI_nec = RM.icrf2itrf('uc', path_to_files, filedate, \
                             time_array, start_date, end_date, inst_GEI)
    # transpose matrices such that they are ready for future calculations
    RRI_enc = RM.NED2ENU(RRI_nec).T
    return RRI_nec, RRI_enc

def inst_direction_in_nec_using_icrf_to_ecef(file_SP3, Roll, Pitch, Yaw,
                     body_inst, GEIx, GEIy, GEIz, Vx, Vy, Vz,
                     time_array, start_date, end_date, lat, lon):
    """
    Parameters
    ----------
    Roll, Pitch, Yaw: Euler angles
    lvlh_x, lvlh_y, lvlh_z: antenna x, y, z vector components in LVLH
    GEIx, GEIy, GEIz: Position of the satellite in GEI
    Vx, Vy, Vz: Speed of satellite in GEI
    time_array: datetime
    lat, lon: spacecraft location in GEO

    Returns
    -------
    ENUx, ENUy, ENUz: pointing direction of 1 monopole in ENU coordinates

    @author: ceren
    """
    # 1- Rotation of SC vectors in LVLH frame
    rotated_inst = RM.rotate_antenna(Roll, Pitch, Yaw, body_inst)
    # %% 2- rotation of SC frame to GEI/ECI frame
    inst_GEI = RM.ORF2J2K_use_spacecraft_ephemeris(rotated_inst, GEIx, GEIy, GEIz, Vx, Vy, Vz)
    # %% 3- rotation of GEI/ECI frame to Earth Centered Earth Fixed frame (ECEF):
    inst_ECEF = RM.GEI2ECEF(time_array, lon, inst_GEI)

    import use_quaternions as uq
    time_sp3, GEO_X, GEO_Y, GEO_Z = ei.sp3_ephemeris_short(file_SP3, start_date, end_date)

    inst_nec=np.empty((len(GEO_X),3));

    for i in range(0, len(GEO_X)):

        rNEC_from_ITRF = RM.Build_NEC_from_ITRF_Rotation_Matrix(GEO_X[i], \
                                                                GEO_Y[i], \
                                                                    GEO_Z[i])

        inst_nec[i, :]  = np.dot( rNEC_from_ITRF,  inst_ECEF[: , i])

    # transform nec to enc
    inst_enc= RM.NED2ENU(inst_nec).T

    return inst_nec, inst_enc

def inst_direction_in_ned_from_itrf(Lat, Lon, inst_ITRF):

    itrf2ned = RM.ECEF2NED(Lat, Lon, inst_ITRF)
    j2k_itrf_ned = itrf2ned.T

    itrf2enu = RM.NED2ENU(j2k_itrf_ned)
    j2k_itrf_enu = itrf2enu.T

    return j2k_itrf_ned, j2k_itrf_enu

def LA_sat(plat, plon, slat, slon, palt, rsat):

    """Calculates the satellite look angle: elevation and azimuth angle that is needed
    for the satellite to look towards a point in GEO frame.

    Parameters
    ----------
    plat, plon: latitude and longitude of the desired pointing location on Earth
                in GEO frame: degrees, float
    slat, slon: latitude and longitude of the subsatellite point
                in GEO frame: degrees, float
    rp, rsat: altitude of the point, altitude of satellite
            in GEO frame: km, float

    Returns
    -------
    la_azim: Azimuth angle for a satellite to look at the given point on Earth
    la_elev: Elevation angle for a satellite to look at the given point on Earth


    References: ASD5: Earth Station Look Angles By Prof. Gregory D. Durgin,
            Astrodynamics, Georgia Tech, 2009.
            """

    rp= 6371 + palt
    la_elev=np.empty(len(slat)) # elevation angle
    la_azim=np.empty(len(slon)) # azimuthal angle

    for i in range(0, len(slat)):
        rslat = slat[i]*(m.pi/180) #convert to radians
        rslon = slon[i]*(m.pi/180) #convert to radians
        rplat = plat*(m.pi/180) #convert to radians
        rplon = plon*(m.pi/180) #convert to radians
        rs = rsat[i] + 6371

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

    return la_elev, la_azim


def LA_inst(vec, n, m, lon, lat, OLon, OLat):
    size = len(lon)
    vec_elev=np.empty(size); vec_elev1=np.empty(size)
    vec_azim=np.empty(size); vec_azim1=np.empty(size)
    horizontal_plane =[0, 0, -1]
    for i in range(0,size):

        vec_azim1[i] = np.arctan2(vec[i][n],vec[i][m])*180.0/np.pi  # azimuth angle
        if vec_azim1[i]<=0:
            vec_azim1[i] = vec_azim1[i] +360

        if lon[i]<OLon: # west of ottawa
            if 0<vec_azim1[i]<=180: # north of Ottawa
                vec_azim[i]=vec_azim1[i]+180
            else:
                vec_azim[i]=vec_azim1[i]+360
        # West of Ottawa
        else: # east of Ottawa
            if 0<vec_azim1[i]<=180:
                vec_azim[i] = 180 - vec_azim1[i]
                # vec_azim[i] = vec_azim1[i] # new
            else:
                vec_azim[i] = 180 + vec_azim1[i]

        horizontal_plane_dot_RRI=np.dot(vec[i][:], horizontal_plane)
        vec_elev1 = np.arccos(horizontal_plane_dot_RRI)*180.0/np.pi
        vec_elev[i]= 90-vec_elev1 # RRI dip angle

    return vec_azim, vec_elev

def calculate_los(pLat, pLon, pAlt, alt, lat, lon):
    SizeArr = len(alt)

    # Line-of-sight distance between the point and satellite
    distance=np.empty(alt.shape);
    # Ray vector
    los_ned=[]; los_enu=[];

    for j in range(1, SizeArr):
        i= j-1

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

    los_enu_arr=np.array(los_enu)
    los_ned_arr=np.array(los_ned)

    return los_enu_arr, los_ned_arr

def spacecraft_distance_from_apoint(time_array, start_time, \
                               pLat, pLon, pAlt, alt, lat, lon):

    SizeArr = len(alt)

    for i in range(0, SizeArr):
    # Line-of-sight distance between the point and satellite
        distance=np.empty(alt.shape);
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

        return distance

def calculate_angles(time_array, file_ray, start_time, \
                                pLat, pLon, pAlt, alt, lat, lon, inst_NED,\
                                    Bottom_ned, ray_model_input):
    SizeArr = len(inst_NED)

    DOY=start_time.timetuple().tm_yday

    # Line-of-sight distance between the point and satellite
    distance=np.empty(alt.shape); distance.fill(np.nan)
    # Ray vector
    los_ned=[]; los_enu=[];

    # angle between RRI and LOS ray path
    obs_los=np.empty(alt.shape); obs_los.fill(np.nan)

    # angle between MF and spacecraft
    offnadir_los=np.empty(alt.shape)

    for i in range(0, SizeArr):
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
        distance[i] = m.sqrt(dx**2 + dy**2 + dz_down**2)

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
    # Find the angle between the LOS ray vector and RRI pointing direction
        dot_productXR_los = np.dot(inst_NED[i][:], los_ray_ned)
        obs_los[i] = m.acos(dot_productXR_los)*180.0/np.pi
    # =============================================================================
    # Find the angle between the nadir vector and RRI los ray vector
        dot_productnadir_los = np.dot(Bottom_ned[i][:], los_ray_ned)
        offnadir_los[i] = m.acos(dot_productnadir_los)*180.0/np.pi

    return offnadir_los, obs_los

def find_slew(offset_los, offset_D1, offset_D2, slew_angle):
    slew = np.empty(len(offset_los))
    slew.fill(np.nan)
    # slew conditions
    for i in range(0, len(offset_los)):
        cond2 = 0 <= offset_los[i] <= slew_angle # slew
        cond2_2 = abs(offset_D1[i])< slew_angle or abs(offset_D2[i])< slew_angle

        ## back slew
        cond_neg2 = -slew_angle <= offset_los[i] <= -0.001 # cond2_2 also holds true here


        ## one dipole front slew (slew = 1)-----------------------------------
        cond1 = -slew_angle <= offset_los[i] <= slew_angle
        cond1_1 = abs(offset_D2[i])> slew_angle and (0 <= offset_D1[i] < slew_angle)
        # dipole 2 front slew
        cond1_2 = (0 <= offset_D2[i] < slew_angle) and abs(offset_D1[i])> slew_angle

        ##  one dipole back slew (slew = -1)----------------------------------
        # dipole 1 back slew
        cond1_3 = abs(offset_D2[i])> slew_angle and (-slew_angle < offset_D1[i] < -0.001)
        # dipole 2 back slew
        cond1_4 = (-slew_angle < offset_D2[i] < -0.001) and abs(offset_D1[i])> slew_angle

        ## one dipole front slew with data containing nans (slew = 1)---------
        # dipole 1 front slew
        cond1_5 = np.isnan(offset_D2[i])==True and (0 < offset_D1[i]< slew_angle)
        # dipole2 2 front slew
        cond1_6 = (0 < offset_D2[i]< slew_angle) and np.isnan(offset_D1[i])==True

        ## one dipole back slew with data containing nans (slew = -1)---------
        # dipole 1 back slew
        cond1_7 = np.isnan(offset_D2[i])==True and (offset_D1[i] < -slew_angle)
        # dipole 2 back slew
        cond1_8 = (offset_D2[i] < -slew_angle) and np.isnan(offset_D1[i])==True

        ## no slew
        cond0_1 = offset_los[i] > slew_angle or offset_los[i] < -slew_angle
        cond0_2 = abs(offset_D1[i])> slew_angle and abs(offset_D2[i])> slew_angle

        # complete front slew
        if cond2 ==True and cond2_2 ==True:
            slew[i] = 2
        # one dipole fromt slew
        elif cond1 == True and (cond1_1==True or cond1_2 == True):
            slew[i] = 1
        # one dipole front slew
        elif cond1 == True and (cond1_5==True or cond1_6 == True):
            slew[i] = 1
        # one dipole back slew
        elif cond1 == True and (cond1_7==True or cond1_8 == True):
            slew[i] = -1
        # one dipole back slew
        elif cond1 == True and (cond1_3==True or cond1_4 == True):
            slew[i] = -1
        # complete back slew
        elif cond_neg2 == True and cond2_2==True:
            slew[i] = -2 # back slew
        elif cond0_1==True and cond0_2== True:
            slew[i]= 0 # no slew
        else:
            slew[i]= 0 # no slew

    return slew
