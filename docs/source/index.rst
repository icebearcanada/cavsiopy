.. cavsiopy documentation master file, created by
   sphinx-quickstart on Wed Mar  8 10:32:43 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to cavsiopy's documentation!
====================================
'cavsiopy' is an open source Python package, which was specifically developed to determine the look direction of the The Radio Receiver Instrument (RRI) on the e-POP/ CASSIOPE/ Swarm-E for HF radio wave propagation studies, but can be applied to any satellite mission requiring accurate pointing information. Given the spacecraft position, the roll, yaw, and pitch angles, the body vector of the instrument, and a target location, the pointing direction and orientation of the instrument onboard the spacecraft with respect to the target can be determined in various reference frames. 'cavsiopy' supports the reference frames: Geocentric Equatorial Inertial J2000 (GEI J2K), Earth Centred Earth Fixed (ECEF), International Terrestrial Reference Frame (ITRF), geodetic North-East-Down, and geocentric North-East-Center. 

cavsiopy contains routines for 
1. Rotating the look direction of an instrument on-board the spacecraft in orbital frame  
2. Calculation of the look angles of the spacecraft (elevation and azimuth) 
3. Calculation of the look angles of the instrument (elevation and azimuth)  
4. Calculation of the distance between the spacecraft and a designated point on the ground  
5. Calculation of the line-of-sight direction vector from the target to the spacecraft location
6. Transformations between GEI J2K, ECEF, NED, NEC, ICRF, ITRF reference frames 
7. Visualization of spacecraft and instrument direction in 2D and 3D (simple or overlaid on geographical regions of the Earth below the satellite trajectory)
8. Rotation matrices for rotations by x, y, z axes

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   installation
   ephemeris_importer
   use_rotation_matrices
   attitude_analysis
   attitude_plotter
   miscellaneous
   
ephemeris_importer: module to import RRI ephemeris and Celestrak TLE
====================================================================
.. automodule:: ephemeris_importer
   :members:   
   
use_rotation_matrices: Module to rotate instrument and transform between reference frames
=========================================================================================
.. automodule:: use_rotation_matrices
   :members:   

attitude_analysis: Module to obtain instrument pointing direction in various reference frames
=============================================================================================
.. automodule:: attitude_analysis
   :members:

attitude_plotter: Module to visualize satellite instrument pointing direction
=============================================================================
.. automodule:: attitude_plotter
   :members:

miscellaneous: Utilities
========================
.. automodule:: miscellaneous
   :members:




