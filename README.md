# cavsiopy: Calculation and Visualization of Spacecraft Instrument Orientation

## Description:
cavsiopy imports the ephemerides and attitude information of the spacecraft and calculates the pointing direction of an instrument onboard.

## Table of Contents:
This package contains routines for
1. Finding the orientation of a spacecraft
2. Finding the look direction of an instrument on-board the spacecraft
3. Calculation of the look angles of the spacecraft (elevation and azimuth)
4. Calculation of the look angles of the instrument (elevation and azimuth)
5. Calculation of the distance between the spacecraft and a designated point on the ground
6. Calculation of the line-of-sight direction vector from the spacecraft to the ground point
7. Transformation routines for the transformations between GEI J2K, ECEF, NED, NEC, ICRF, ITRF reference frames.
8. Visualization of spacecraft and instrument direction in 2D and 3D (simple or overlaid on geographical regions of the Earth below the satellite trajectory).
9. Rotation matrices for rotations by x, y, z axes

The pointing direction vectors can be obtained in GEI J2K, ECEF, NED, NEC, ICRF, ITRF.

## Requirements:
numpy, matplotlib, astropy, cartopy, geopy, h5py

## Installation:
before installing requirements: cartopy requires the below for pip
sudo apt -y install libgeos-dev
conda install -c conda-forge cartopy
======================================================================================
before installing pysofa (https://kloppenborg.net/blog/building-sofa-for-pysofa/):
compile shared c library

sudo apt-get install python python-all-dev build-essential cmake

Download and extract SOFA from the official download page

After extracting SOFA, cd into the main directory and create a CMakeLists.txt file with the following content:

----------------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6)

project(sofa_c C)
  
# Set a few variables:
set(LIBS ${LIBS} m)

# Extract all of the source files.
file(GLOB_RECURSE C_SOURCE . src/*.c)

# Build a shared library
add_library(sofa_c SHARED ${C_SOURCE})
  
# Now define the installation options
install(TARGETS sofa_c LIBRARY DESTINATION lib)
----------------------------------------------------------------------------------------

cmake .
make
make install

sudo ldconfig
========================================================================================
pysofa related problems:
========================================================================================

ModuleNotFoundError: No module named 'pysofa_ctypes'
copy the contents of pysofa_ctypes to __init__.py

-----------------------------------------------------------------------------------------

if __sofa_version < (2010, 12, 01):
                                   ^
SyntaxError: leading zeros in decimal integer literals are not permitted; use an 0o prefix for octal integers

find and replace all 2010, 12, 01 with 2010, 12, 1 in __init__.py
=========================================================================================

pip install -r requirements.txt

finally:
pip install cavsiopy

## Usage:
Given the information of instrument body vector with respect to spacecraft in spacecraft body frame X, Y and Z axes, users can obtain the pointing direction of an instrument on-board spacecraft with this python package. For antenna on spacecraft, this information is useful for finding the deviation of the antenna boresight from the line-of-sight signal transmitted from a ground transmitter. For imagers, the user can easily find the ground coverage of the imager if the field of view of the imager is known.

Credits: C. Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley Martin
Contributing: Glenn C. Hussey, Robert Gillies, Andrew W. Yau

License: Apache License V2.0
