# cavsiopy: Calculation and Visualization of Spacecraft Instrument Orientation

## Description:
cavsiopy imports the ephemerides and attitude information (roll, pitch, yaw) of the spacecraft and determines the pointing direction of an instrument onboard.

## Table of Contents:
This package contains routines for
1. Finding the orientation of a spacecraft (spacecraft_attitude_xyz)
2. Rotating the look direction of an instrument on-board the spacecraftin orbital frame (rotate_instrument)
3. Calculation of the look angles of the spacecraft (elevation and azimuth) (LA_sat)
4. Calculation of the look angles of the instrument (elevation and azimuth) (LA_inst)
5. Calculation of the distance between the spacecraft and a designated point on the ground (spacecraft_distance_from_a_point)
6. Calculation of the line-of-sight direction vector from the target to the spacecraft location (calculate_los_vec)
7. Transformations between GEI J2K, ECEF, NED, NEC, ICRF, ITRF reference frames (find_instrument_attitude)
8. Visualization of spacecraft and instrument direction in 2D and 3D (simple or overlaid on geographical regions of the Earth below the satellite trajectory) (attitude_plotter module)
9. Rotation matrices for rotations by x, y, z axes

Instrument pointing direction vectors can be obtained in GEI J2K, ECEF, NED, NEC, ICRF, ITRF.

## Requirements:
numpy, matplotlib, astropy, cartopy, geopy, h5py, pysofa

## Installation:
before installing cavsiopy some of the dependencies need special installation instructions:

### cartopy:
* If you are installing cartopy via pip (on Ubuntu/Debian) first libgeos-dev has to be on your system:

  sudo apt -y install libgeos-dev

  pip3 install cartopy 

* Follow the installation guide for other systems: https://scitools.org.uk/cartopy/docs/latest/installing.html

* for installation of cartopy with conda:

  conda install -c conda-forge cartopy

### pysofa:
pysofa is the Python wrapper for the Standards of Fundamental Astronomy (SOFA)  written by Brian Kloppenborg.
cavsiopy depends on pysofa for determining the ICRF to ITRF rotation matrix.

Here we follow (https://kloppenborg.net/blog/building-sofa-for-pysofa/) for pysofa installation:

sudo apt-get install python python-all-dev build-essential cmake

Download and extract SOFA from the official download page (https://www.iausofa.org/2010_1201_C/sofa_c-20101201.tar.gz)

After extracting SOFA, cd into the main directory and create a CMakeLists.txt file with the following content:

cmake_minimum_required(VERSION 2.6)

project(sofa_c C)
  
####  Set a few variables:
set(LIBS ${LIBS} m)

#### Extract all of the source files:
file(GLOB_RECURSE C_SOURCE . src/*.c)

#### Build a shared library:
add_library(sofa_c SHARED ${C_SOURCE})
  
#### Now define the installation options:
install(TARGETS sofa_c LIBRARY DESTINATION lib)

cmake .

make

make install

sudo ldconfig

#### you may need some more tweaking for pysofa

1. ModuleNotFoundError: No module named 'pysofa_ctypes'

  copy the contents of pysofa_ctypes to __init__.py

2. if __sofa_version < (2010, 12, 01):
                                   ^
  SyntaxError: leading zeros in decimal integer literals are not permitted; use an 0o prefix for octal integers

  find and replace all 2010, 12, 01 with 2010, 12, 1 in __init__.py

### other packages and requirements:
pip install --no-deps astropy

pip install geopy

pip install -r requirements.txt

### now, we are good to go!

pip install cavsiopy


Credits: C. Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley Martin

Contributing: Glenn C. Hussey, Robert Gillies, Andrew W. Yau

License: GNU LESSER GENERAL PUBLIC LICENSE Version 3, 29 June 2007
