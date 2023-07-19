# cavsiopy: Calculation and Visualization of Spacecraft Instrument Orientation

Welcome to cavsiopy!
====================================
'cavsiopy' is an open source Python package, which was specifically developed to determine the look direction of the Radio Receiver Instrument (RRI) on the e-POP/ CASSIOPE/ Swarm-E for HF radio wave propagation studies, but can be applied to any satellite mission requiring accurate pointing information. Given the spacecraft position, the roll, yaw, and pitch angles, the body vector of the instrument, and a target location, the pointing direction and orientation of the instrument onboard the spacecraft with respect to the target can be determined in various reference frames. 'cavsiopy' supports the reference frames: Geocentric Equatorial Inertial J2000 (GEI J2K), Earth Centred Earth Fixed (ECEF), International Terrestrial Reference Frame (ITRF), geodetic North-East-Down, and geocentric North-East-Center. 

cavsiopy contains routines for 

1. Rotating the look direction of an instrument on-board the spacecraft in orbital frame 

2. Calculation of the look angles of the spacecraft (elevation and azimuth) 

3. Calculation of the look angles of the instrument (elevation and azimuth) 

4. Calculation of the distance between the spacecraft and a designated point on the ground 

5. Calculation of the line-of-sight direction vector from the target to the spacecraft location

6. Transformations between GEI J2K, ECEF, NED, NEC, ICRF, ITRF reference frames 

7. Visualization of spacecraft and instrument direction in 2D and 3D (simple or overlaid on geographical regions of the Earth below the satellite trajectory)

8. Rotation matrices for rotations by x, y, z axes

Instrument pointing direction vectors can be obtained in GEI J2K, ECEF, NED, NEC, ICRF, ITRF.

Dependencies
------------
numpy, matplotlib, astropy, cartopy, geopy, h5py, pysofa

Installation
------------
Before installing cavsiopy some of the dependencies need special installation instructions:

# cartopy

- If you are installing cartopy via pip (on Ubuntu/Debian) first libgeos-dev has to be on your system:

```
$ sudo apt -y install libgeos-dev
$ pip install -r requirements.txt
```

- Follow the installation guide for other systems: https://scitools.org.uk/cartopy/docs/latest/installing.html

- for installation of cartopy with conda:

```
$ conda install -c conda-forge cartopy
```
# pysofa

pysofa is the Python wrapper for the Standards of Fundamental Astronomy (SOFA) written by Brian Kloppenborg.
cavsiopy depends on pysofa for determining the ICRF to ITRF rotation matrix.

Here we follow <https://kloppenborg.net/blog/building-sofa-for-pysofa/> for pysofa installation.

```
$ sudo apt-get install python python-all-dev build-essential cmake
```
* Download and extract SOFA from the official download page (https://www.iausofa.org/2010_1201_C/sofa_c-20101201.tar.gz)

* After extracting SOFA, cd into the main directory and create a CMakeLists.txt file with the following content

```
cmake_minimum_required(VERSION 2.6)
project(sofa_c C)
	
# Set a few variables
set(LIBS ${LIBS} m)
	
# Extract all of the source files
file(GLOB_RECURSE C_SOURCE . src/"\\*".c)

# Build a shared library
add_library(sofa_c SHARED ${C_SOURCE})
	
# Now define the installation options
install(TARGETS sofa_c LIBRARY DESTINATION lib)
```
   
Save and close the file.

Execute the following comments:

```
$ cmake .

$ make

$ make install

$ sudo ldconfig
```

* you may need some more tweaking for the installation of pysofa. Possible solutions are listed below:

1. ModuleNotFoundError: No module named 'pysofa_ctypes'

  copy the contents of pysofa_ctypes to __init__.py

2. if __sofa_version < (2010, 12, 01):
                                   ^
  SyntaxError: leading zeros in decimal integer literals are not permitted; use an 0o prefix for octal integers

  find and replace all 2010, 12, 01 with 2010, 12, 1 in __init__.py

Other packages and requirements
-------------------------------
```
pip install --no-deps astropy

pip install geopy

pip install -r requirements.txt
```
** Now, we are good to go! **

```
pip install cavsiopy
```

Credits: C. Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley Martin

Contributing: Glenn C. Hussey, Robert Gillies, Andrew W. Yau

License: GNU LESSER GENERAL PUBLIC LICENSE Version 3, 29 June 2007
