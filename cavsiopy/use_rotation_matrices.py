#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 14:32:10 2021

1. Calculates satellite orbital elements using X, Y, Z coordinates and \
    Vx, Vy, Vz in GEI
2. Calculates GMST time from noon or midnight
3. Rotates the instrument in orbital frame using roll, pitch, and yaw angles:
    rotate_antenna
4. Includes modules for transformations from
    Spacecraft orbital frame to GEI J2000:
        ORF2J2K_use_orbital_elements, ORF2J2K_use_spacecraft_ephemeris
    J2K to ECEF:
        GEI2ECEF
    J2K to ITRF:
        icrf2itrf
    ITRF to NED:
        ECEF2NED
    NED to ENU or NEC to ENC:
        NED2ENU
5. Computes the rotation matrice from ITRF to NEC:
        Build_NEC_from_ITRF_Rotation_Matrix

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
"""
import numpy as np
from astropy.time import Time
import cavsiopy.ephemeris_importer as ei
import pysofa
import cavsiopy.misc as misc

def orbital_elements(X, Y, Z, Vx, Vy, Vz):
    """
    Parameters
    ----------
        X, Y, Z: float
            X, Y, Z in GEI
        Vx, Vy, Vz: float
            Velocity components in GEI

    Returns
    -------
        raan: float
            Right ascension of ascending node
        inc: float
            Spacecraft inclination
        ap: float
            Argument of periapsis
        TA: float
            True anomaly

    Examples
    --------
    e, TA, raan, inc, ap = orbital_elements(file_RRI)

    Reference: Curtis, H. D. (2014). Orbits in Three Dimensions in
    Orbital mechanics for engineering students. Butterworth-Heinemann."""

    inc=np.empty(X.shape); inc.fill(np.nan)
    raan=np.empty(X.shape); raan.fill(np.nan)
    ap=np.empty(X.shape); ap.fill(np.nan)
    TA=np.empty(X.shape); TA.fill(np.nan)
    E=np.empty(X.shape); E.fill(np.nan)
    Praan=np.empty(X.shape); Praan.fill(np.nan)
    for i in range (0,len(X)):
        # distance
        d = np.sqrt(X[i]**2 + Y[i]**2 + Z[i]**2)
        d_vec = [X[i], Y[i], Z[i]]
        # speed
        V = np.sqrt(Vx[i]**2 + Vy[i]**2 + Vz[i]**2)
        V_vec = [Vx[i], Vy[i], Vz[i]]

        # Vrad < 0: satellite is flying towards perigee
        # Vrad > 0: satellite is flying away from perigee
        Vrad = (X[i]*Vx[i] + Y[i]* Vy[i] + Z[i]* Vz[i]) / d

        # components of angular momentum vector
        hx = Y[i] * Vz[i] - Vy[i] * Z[i]
        hy = Vx[i] * Z[i] - X[i] * Vz[i]
        hz = X[i] * Vy[i] - Vx[i] * Y[i]
        h_vec=[hx, hy, hz]

        h= np.linalg.norm(h_vec)

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

        e_vec_x = e_d * X[i] + e_v * Vx[i]
        e_vec_y = e_d * Y[i] + e_v * Vy[i]
        e_vec_z = e_d * Z[i] + e_v * Vz[i]
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

    return e, TA, raan, inc, ap


def calculate_GMST_midnight(utc_dt_array):
  """
   Calculates GMST in radians.

   Includes day fraction corrected with respect to midnight.

   Follows: Curtis, H. D. (2013). Orbital mechanics for engineering students.
    Butterworth-Heinemann. This formula gives GMST in degrees.

   Parameters
   ----------
   utc_dt_array: datetime
       time array in UTC

   Returns
   -------
   GMST: float
       GMST in radians

   Examples
   --------
   GMST = calculate_GMST_midnight(utc_dt_array)
   """

  # Astropy allows JD/MJD calculations to be done using UTC, oddly.
  #   Time().utc.jd != Time().tt.jd
  GMST = np.zeros(len(utc_dt_array), dtype=np.float)

  for i in range(len( utc_dt_array )):
    JD = Time( utc_dt_array[i] ).ut1.jd
    T  = (JD - 2451545)/36525 # Julian Centuries since J2K epoch in TT.
    UT = JD % 1 - 0.5   # Day fraction, correcting to midnight rather than noon.

    gmst = 100.460618375 + 36000.770053608336*T + 0.0003879333*T**2 - 2.583e-8*T**3

    E_deg = 360.98564724 * UT # /24 built-in, UT is just the fractional day.

    Theta_deg = gmst + E_deg
    GMST[i] = ( Theta_deg * np.pi/180 ) % (2*np.pi)

  return GMST

def calculate_GMST_noon(time_array):
    """ Follows: Curtis, H. D. (2013). Orbital mechanics for engineering students.
     Butterworth-Heinemann. This formula gives GMST in degrees.

     Includes day fraction corrected with respect to noon.

     Parameters
     ----------
     utc_dt_array: datetime
         time array in UTC

     Returns
     -------
     GMST: float
         GMST in radians

     Examples
     --------
     GMST = calculate_GMST_noon(time_array)
     """
    SizeArr=len(time_array)
    Theta_rad=np.empty(SizeArr)

    for i in range(0, SizeArr):
        # time of epoch
        JD=Time(time_array[i]).ut1.jd  # jd: Julian date
        T = (JD - 2451545)/36525;  # Time elapsed since JD 2451545 (Jan 1, 2000, 12 UT)

        # Reference: Curtis, H. D. (2013). Orbital mechanics for engineering students.
        # Butterworth-Heinemann. This formula gives GMST in degrees.

        GMST_deg = 100.460618375 + 36000.770053608336 * T + 0.0003879333 * (T**2) - \
            2.583 * (10**(-8)) * (T**3)

        # Formula above can yield values greater than 360. If that's the case:
        if GMST_deg < 0:
            GMST = GMST_deg + np.fix(GMST_deg/360)*360
        elif GMST_deg >= 360:
            GMST = GMST_deg - (np.fix(GMST_deg/360) - 1)*360
        else:
            GMST = GMST_deg

        # from Earth's rotation
        UT_H = time_array[i].timetuple().tm_hour +\
            time_array[i].timetuple().tm_min/60 +\
                ( time_array[i].timetuple().tm_sec / 3600 )

        E_deg = 360.98564724 *  (UT_H/24)

        Theta_deg = GMST + E_deg

        Theta_rad[i] = Theta_deg * np.pi/180

    return Theta_rad
# =============================================================================
# SET UP ROTATION MATRICES
# =============================================================================
# SET-1: According to Coordinate Transformations via Euler Angles, Riggs, 2019. Rev. E.
# expresses the inertial frame vectors in terms of rotated frame vectors
# used for transformation from a rotated frame to inertial frame
# gives i, j, k in terms of i', j', k'; abbreviated as r2i in the code
# definition of the reference: x (out), y (right), and z (up)
# ROLL: positive is from +Y towards +Z
# PITCH: positive is from +Z towards +X
# YAW: positive is from +X towards +Y
# =============================================================================
def RX_r2i(x):
    """frame rotation about x axis-roll"""
    return np.array([[1,      0,       0],
                      [0, np.cos(x), -np.sin(x)],
                      [0, np.sin(x),  np.cos(x)]])

def RY_r2i(y):
    """
    Frame rotation about y axis-pitch

    Parameters
    ----------
    y: float
       pitch angle in radians

    Returns
    -------
    Ry: float
        rotation matrice along the Y axis
    Examples
    -------
    Ry = RY_r2i(pitch_angle)
    """
    return np.array([[np.cos(y), 0, np.sin(y)],
                      [0,      1,        0],
                      [-np.sin(y), 0, np.cos(y)]])

def RZ_r2i(z):
    """
    Frame rotation about z axis-Yaw

    Parameters
    ----------
    z: float
       yaw angle in radians

    Returns
    -------
    Rz: float
        rotation matrice along the Z axis

    Examples
    --------
    Rz = RZ_r2i(yaw_angle)
    """
    return np.array([[np.cos(z), -np.sin(z), 0],
                      [np.sin(z), np.cos(z), 0],
                      [0,      0,       1]])

# =============================================================================
# SET-2: According to Curtis, H. D. (2013). Orbital mechanics for engineering students.
# expresses the rotated frame vectors in terms of inertial frame vectors
# for transformation from an inertial frame to rotated frame
# gives i', j', k' in terms of i, j, k; abbreviated as i2r in the code
# RRI description: x (out), y (left), and z (down)
# ROLL: positive is from +Y towards +Z
# PITCH: positive is from +Z towards +X
# YAW: positive is from +X towards +Y
# =============================================================================
def RX_i2r(x):
    """frame rotation about x axis-roll"""
    return np.array([[1,      0,       0],
                      [0, np.cos(x), np.sin(x)],
                      [0, -np.sin(x),  np.cos(x)]])

def RY_i2r(y):
    """frame rotation about y axis-pitch"""
    return np.array([[np.cos(y), 0, -np.sin(y)],
                      [0,      1,        0],
                      [np.sin(y), 0, np.cos(y)]])

def RZ_i2r(z):
    """frame rotation about z axis-yaw"""
    return np.array([[np.cos(z), np.sin(z), 0],
                      [-np.sin(z), np.cos(z), 0],
                      [0,      0,       1]])

# =============================================================================
# functions for antenna body rotations
# =============================================================================
def rotate_antenna(Roll, Pitch, Yaw, x):
   """
   rotate_monopole(Roll,  Pitch, Yaw, x, y, z)

   Parameters
   ----------
   Roll, Yaw, Pitch: in degrees, float
   x, y, z : initial position vectors for antenna in LVLH, tuple [xx, xy, xz]

   Returns
   -------
   rotated_body: rotated body vectors wrt LVLH,
                           numpy array of (3 x body vector length)

   Rotation of body vectors in LVLH(Local vertical Local horizontal) frame
   using roll, yaw, pitch angles provided by RRI data.

   The rotation sequence is X->Y->Z; roll->pitch->yaw according to the RRI
   data description.

   Reference: University of Calgary, e-POP Radio Receiver Instrument (RRI)
   User’s Manual, Doc. no. ePOP-5024, Rev. D (2018)

   """

   # converting degrees to radians
   rr = Roll*(np.pi/180)
   pr = Pitch*(np.pi/180)
   yr = Yaw*(np.pi/180)

   # LVLH (Local vertical local horizontal) 3 axis unit vectors
   lvlh= np.array([x,] *len(Roll)).T

   # Initialize body Vectors
   rot_vec = np.empty(lvlh.shape)

   # transform initial definition of LVLH to rotated frame
   """ body_vectors = (RzRyRx)(lvlh_vectors)- roll->pitch->yaw sequence"""
   for i in range(Roll.shape[0]):

       # direction cosine matrice
       R = RX_r2i(rr[i]) @ RY_r2i(pr[i]) @ RZ_r2i(yr[i])

       # rotation
       rot_vec[:, i] = R @ lvlh[:, i]

   rotated_body=np.array(rot_vec)


   return rotated_body

# =============================================================================
# # %% function for PQW to GEI/ECI transformations
# # GEI: GEocentric Equatorial Inertial, ECI: Earth Centered Inertial
# =============================================================================
def ORF2J2K_use_orbital_elements(body, GEI):
    """ Uses the inverse of the GEI2LVLH matrice, however,
    rotated to inertial frame matrices can also be used.

    Ref: Frame rotations and quaternions, Pedro A. Capó-Lugo, Peter M. Bainum,
    in Orbital Mechanics and Formation Flying, 2011

    Parameters
    ----------
    body: rotated body frame vectors
    file_RRI: name of the RRI file

    Returns
    -------
    SC_GEIx, SC_GEIy, SC_GEIz: body vectors in GEI coordinates,
                           numpy array of (3 x body vector length)
    raan: Right ascension of ascending node, radians, float
    inc: Satellite inclination, radians, float
    ap: Argument of Perigee, radians, float

    Examples
    --------
    inst_GEI = ORF2J2K_use_orbital_elements(inst_body_vector, GEI)

    """

    # Initialize GEI/ECI Vectors
    SC_OE2GEI = np.empty(body.shape)
    SC_OE2GEI2 = np.empty(body.shape)

    # P= np.array([[0, 1, 0],
    #               [-1, 0, 0],
    #               [0, 0, 1]])
    #  obtain orbital elements
    e, TA, raan, inc, ap= orbital_elements(GEI)

    for i in range(0, body.shape[1]):

        # Canuto2018, Chapter 3
        P= np.array([[0, 0, -1],
                      [1, 0, 0],
                      [0,-1, 0]])

        LVLH2OE = RZ_r2i(TA[i]) @ P @  body[:, i]

        OE2GEI = (RZ_r2i(raan[i]) @ RX_r2i(inc[i]) @ RZ_r2i(ap[i]))

        SC_OE2GEI[:, i] = OE2GEI @ (LVLH2OE)
        #  normalize
        denominator = np.sqrt(SC_OE2GEI[0, i]**2 + SC_OE2GEI[1, i]**2 + \
                              SC_OE2GEI[2, i]**2)
        SC_OE2GEI2[:, i]= SC_OE2GEI[:, i] / denominator

    inst_GEI=np.array(SC_OE2GEI2)

    return inst_GEI

# =============================================================================
# # %% function for GEI to ECEF transformations
# =============================================================================
def GEI2ECEF(time_array, lon, SC_GEI):
    """
    GEI2ECEF(time_array, lon, SC_GEIx, SC_GEIy, SC_GEIz)

    Parameters
    ----------
    time_array: time array, datetime
    SC_GEIx, SC_GEIy, SC_GEIz: X, Y, Z unit vectors in GEI

    Returns
    -------
    SC_ECEFx, SC_ECEFy, SC_ECEFz: X, Y, Z unit vectors in ECEF

    For this, one rotation is needed in xy frame (yaw rotation)
    yaw angle is the azimuth angle between ECI and ECEF in this case.
    Azimuth= Greenwich sidereal time + Earth's rotation speed * UTC.
    References:
    1. J. Riggs, Coordinate transformations via Euler Angle Rotations,
    Rev. E - 10 April 2019
    2. A. Sofyali, Orbital Mechanics notes, 7 February 2019. """

    SizeArr = np.size(time_array)

     # Initialize GEI/ECI Vectors
    SC_GEI2ECEF = np.empty(SC_GEI.shape)
    SC_GEI2ECEF2 = np.empty(SC_GEI.shape)

    Theta_rad = calculate_GMST_midnight(time_array)


    for i in range(0, SizeArr):

        SC_GEI2ECEF[:, i] = RZ_i2r(Theta_rad[i]) @ SC_GEI[:, i]

        below= np.sqrt(SC_GEI2ECEF[0, i]**2+SC_GEI2ECEF[1, i]**2+SC_GEI2ECEF[2, i]**2)

        SC_GEI2ECEF2[:, i]= SC_GEI2ECEF[:, i]/below


    SC_ECEF=np.array(SC_GEI2ECEF2)


    return SC_ECEF

# =============================================================================
# # %% function for ECEF to NED: (North, East, Down) transformations
# # z is positive downwards
# =============================================================================
def ECEF2NED(lat, lon, SC_ECEF):
    """
    ECEF2NED(lat, lon, SC_ECEF)

    Parameters
    ----------
    lat: latitude, in degrees, float
    lon: longitude, in degrees, float
    SC_ECEF: [SC_ECEFx, SC_ECEFy, SC_ECEFz] X, Y, Z unit vectors in ECEF, float

    Returns
    -------
    SC_NEDx, SC_NEDy, SC_NEDz: X, Y, Z unit vectors in NED

    Y-Z (pitch-yaw) sequence with longitude and latitude are needed from
    transformation from ECEF to NED
    1. Aligning the x-y plane of ECEF and NED using longitude: yaw(longitude)
    2. Aligning the x-z plane using latitude: pitch(-(pi/2+latitude))
    Reference: J. Riggs, Coordinate transformations via Euler Angle Rotations,
    Rev. E - 10 April 2019

    """
    # Initialize NED Vectors
    SC_ECEF2NED = np.empty(SC_ECEF.shape)
    SC_ECEF2NED2=np.empty(SC_ECEF.shape)

    for i in range(0, SC_ECEF.shape[1]):
        # Open form ----------------------------------------------------------
        lamda = lon[i]*(np.pi/180)
        psi = lat[i]*(np.pi/180)
        a = -np.cos(lamda)*np.sin(psi)
        b = -np.sin(lamda)
        c = -np.cos(lamda)*np.cos(psi)
        d = -np.sin(lamda)*np.sin(psi)
        e = np.cos(lamda)
        f = -np.sin(lamda)*np.cos(psi)
        g = np.cos(psi)
        h = 0
        m = -np.sin(psi)

        RE2N = np.array([[a, d, g],
                      [b, e, h],
                      [c, f,  m]])
        # -------------------------------------------------------------------

        SC_ECEF2NED[:, i] = RE2N @ SC_ECEF[:, i]

        denominator = np.sqrt(SC_ECEF2NED[0, i]**2 + SC_ECEF2NED[1, i]**2 + \
                              SC_ECEF2NED[2, i]**2)

        SC_ECEF2NED2[:, i]= SC_ECEF2NED[:, i] / denominator

    return SC_ECEF2NED2

# =============================================================================
# # %% function for NED to ENU transformations
# =============================================================================
def NED2ENU(inst_NED):
    """
    NED2ENU(inst_NED)

    Parameters
    ----------
    inst_NED: body vector in NED

    Returns
    -------
    inst_ENU: body vector in ENU

    The transformation matrice from NED to ENU and ENU to NED are the same.
    [0 1 0 ; 1 0 0 , 0 0 -1]
    in ENU, z increases in upward direction.

    Reference:  Grewal, M. S., Weill, L. R., & Andrews, A. P. (2008).
    Appendix C: Coordinate Transformations. In Global Positioning Systems,
    Inertial Navigation, and Integration (pp. 456–501). John Wiley & Sons, Inc.
    https://doi.org/10.1002/9780470099728.app3
"""
    # # Initialize ECEF Vectors
    SC_NED=np.array(inst_NED).T

    # Initialize ENU Vectors
    SC_NED2ENU = np.empty(SC_NED.shape)

    # rotation matrice from NED to ENU
    RN2E= np.array([[0, 1, 0],
                    [1, 0, 0],
                    [0, 0, -1]])

    SC_NED2ENU = RN2E @ SC_NED

    return SC_NED2ENU

def J2K_Ephemeris_to_J2K_to_Nadir_Rotation_Matrix(pX,pY,pZ,vX,vY,vZ):
  '''
  Given 6 floats of the J2K/ICRF position and velocity data,
  Returns a 3x3 rotation matrix describing the rotation of J2K->S/C (or Nadir).
  '''
  # Rebuild the Position and Velocity vectors in ICRF/J2K
  Position = np.array([pX,pY,pZ])
  Velocity = np.array([vX,vY,vZ])

  zVector = -Position  # Invert vector, as Z is satellite-down rather than center-of-earth-up as given in the position.
  yVector = np.cross(zVector, Velocity)
  xVector = np.cross(yVector, zVector)

  # Normalize the vectors:
  xVector = xVector / np.linalg.norm(xVector)
  yVector = yVector / np.linalg.norm(yVector)
  zVector = zVector / np.linalg.norm(zVector)

  # numpy.dstack stacks side-by-side, instead of numpy.array's stack, annoyingly.
  # But it's still convenient for rotation-matricies.
  # Build the rotation matrix, which due to the side-by-side construction, is Nadir->J2K
  rNadir_to_J2K = np.dstack([xVector, yVector, zVector])[0]
  # Transpose to get the J2K->Nadir rotation matrix.
  rJ2K_to_Nadir = rNadir_to_J2K.T

  return rJ2K_to_Nadir, rNadir_to_J2K

def Build_NEC_from_ITRF_Rotation_Matrix( pX, pY, pZ ):
  """ Given the Geographic Cartesian Position ephemeris data from get_geographic_cartesian_coordinates() below,
   Return a 3d notation matrix describing the rotation from NEC<-ITRF
   This is given in the SWARM documentation:
   https://earth.esa.int/documents/10174/1514862/Swarm_Level-1b_Processor_Algorithms
   Page 137, "C.3 North East Center (NEC)"""

  Position_Vector = np.array([ pX,pY,pZ ])

  eCenter = -Position_Vector
  eEast   = np.cross( eCenter, [0,0,1] ) # The cross product of geographic north-pole vector and the center vector.
  eNorth  = np.cross( eEast, eCenter )

  # And then normalize the vectors:
  eCenter /= np.linalg.norm(eCenter)
  eEast   /= np.linalg.norm(eEast  )
  eNorth  /= np.linalg.norm(eNorth )

  return np.stack(( eNorth, eEast, eCenter)) # Vertically stack them, and return.

def ORF2J2K_use_spacecraft_ephemeris(lvlh, GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz):
    """
    ORF2J2K_use_spacecraft_ephemeris(lvlh, GEIx, GEIy, GEIz, GEIVx, GEIVy, GEIVz)
    Converts from the orbital frame to GEI J2000.

    Parameters
    ----------
    GEIx, GEIy, GEIz: x, y, z positions in GEI coordinates (row x) numpy array
    Vx, Vy, Vz: X, Y, Z components of velocity in GEI coordinates (row x) numpy array

    Returns
    -------
    SC_GEIx, SC_GEIy, SC_GEIz: body vectors in GEI coordinates,
                           numpy array of (3 x body vector length)

    """

    # Initialize GEI/ECI Vectors
    SC_LVLH_J2K = np.empty(lvlh.shape)
    SC_LVLH_J2K2 = np.empty(lvlh.shape)

    # Z-X-Z sequence

    for i in range(0, lvlh.shape[1]):
        rJ2K_to_Nadir, rNadir_to_J2K= \
            J2K_Ephemeris_to_J2K_to_Nadir_Rotation_Matrix(GEIx[i], GEIy[i],
                                                          GEIz[i], GEIVx[i],
                                                          GEIVy[i], GEIVz[i])
        R = rNadir_to_J2K

        SC_LVLH_J2K[:, i] = R @ (lvlh[:, i])
        #  normalize the matrix
        denominator= np.sqrt(SC_LVLH_J2K[0, i]**2 + SC_LVLH_J2K[1, i]**2 + \
                             SC_LVLH_J2K[2, i]**2)

        SC_LVLH_J2K2[:, i]= SC_LVLH_J2K[:, i]/denominator


    SC_J2K=np.array(SC_LVLH_J2K2)

    return SC_J2K

def icrf2itrf(path_to_files, filedate, MET, GEO, \
              time_array, start_date, end_date, inst_ICRF):
    """
    icrf2itrf(path_to_files, filedate, MET, GEO, \
                  time_array, start_date, end_date, inst_ICRF)
    Parameters
    ----------
    path_to_files: str
    Vx, Vy, Vz: X, Y, Z components of velocity in GEI coordinates (row x) numpy array

    Returns
    -------
    SC_GEIx, SC_GEIy, SC_GEIz: body vectors in GEI coordinates,
                           numpy array of (3 x body vector length)
    Notes
    -----
    using Standards of Fundamental Astronomy library
    https://www.iausofa.org/

    """

    GEOx, GEOy, GEOz = GEO

    row, col = np.shape(inst_ICRF)
    inst_ITRF_Vector=np.empty((col,3)); inst_ITRF_Vector.fill(np.nan)

    year = time_array[0].year

    # calculate day-of-year = DOY
    DOY=start_date.timetuple().tm_yday
    date_year=start_date.year
    if (date_year%4 == 0):
        doy = DOY / 366
        # calculate year + doy to match with IERS data
        ydoy = year + doy
    # year in IERS tabulations for TT-UT1 changes with 0.05
        inc = 0.05 * 366
    else:
        doy = DOY / 365.25
    # calculate year + doy to match with IERS data
        ydoy = year + doy
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

    for i in range(0, col):

        JD=Time(time_array[i]).jd  # jd: Julian date
        dT_JD = (JD - 2451545)/36525;  # Time elapsed since JD 2451545 (Jan 1, 2000, 12 UT)

        utc1 = np.floor(JD)
        utc2= JD-utc1

        MJD = int(np.floor(Time(time_array[i]).mjd))

        #  IERS tabulations
        #  delta UT (UT1-UTC) can be obtained from the link below
        # link = "https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt"
        filename_IERS = 'IERS_daily_deltaT.txt'
        file_IERS = path_to_files + filename_IERS
        IER=np.loadtxt(fname=file_IERS, skiprows=14, comments='#')

        # find where MJD equals the MJD in the text file
        index_day = np.where(IER== MJD)[0]
        # delta UT from IERS file
        dut = IER[index_day, 7]

        # find the terrestrial time from coordinated universal time--------
        #  dut = ut1-utc
        uta, utb = pysofa.utcut1(utc1, utc2, dut)

    #----------------------------------------------------------------------
        #  dt = ut-tt
        tta, ttb = pysofa.ut1tt(uta, utb, dt)

        #  find the x_p, y_p, s: celestial pole coordinates----------------
        x_p, y_p, s = pysofa.xys00b(tta, ttb)

        #  find the celestial to terrestrial matrix------------------------
        c2t = pysofa.c2t00b(tta, ttb, uta, utb, x_p, y_p)

        inst_ITRF_Vector[i, :] = c2t @ inst_ICRF[:,i]

    inst_ITRF = inst_ITRF_Vector.T

    # transform inst_ITRF to inst_NEC using quaternions
    inst_NEC_Vector=np.empty((col,3)); inst_NEC_Vector.fill(np.nan)

    for i in range(0, col):

        rNEC_from_ITRF = Build_NEC_from_ITRF_Rotation_Matrix(GEOx[i], GEOy[i], GEOz[i])

        inst_NEC_Vector[i, :]  = np.dot( rNEC_from_ITRF,  inst_ITRF[: , i])

    return inst_ITRF, inst_NEC_Vector
