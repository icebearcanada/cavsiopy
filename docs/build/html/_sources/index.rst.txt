.. cavsiopy documentation master file, created by
   sphinx-quickstart on Wed Mar  8 10:32:43 2023.

Welcome to cavsiopy's documentation!
====================================

'cavsiopy' is an open source Python package, which was specifically developed to determine the look direction of the The Radio Receiver Instrument (RRI) on e-POP/ CASSIOPE/ Swarm-E for HF radio wave propagation studies, but can be applied to any satellite mission requiring accurate pointing information. Given the spacecraft position, the roll, yaw, and pitch angles, the body vector of the instrument, and a target location, the pointing direction and orientation of the instrument onboard the spacecraft with respect to the target can be determined in various reference frames. 'cavsiopy' supports the reference frames: Geocentric Equatorial Inertial J2000 (GEI J2K), Earth Centred Earth Fixed (ECEF), International Terrestrial Reference Frame (ITRF), geodetic North-East-Down, and geocentric North-East-Center. 

cavsiopy contains routines for 

1. Rotating the look direction of an instrument on-board the spacecraft in orbital frame  

2. Calculation of the look angles of the spacecraft (elevation and azimuth) 

3. Calculation of the look angles of the instrument (elevation and azimuth)  

4. Calculation of the distance between the spacecraft and a designated point on the ground  

5. Calculation of the line-of-sight direction vector from the target to the spacecraft location

6. Transformations between GEI J2K, ECEF, NED, NEC, ICRF, ITRF reference frames 

7. Visualization of spacecraft and instrument direction in 2D and 3D (simple or overlaid on geographical regions of the Earth below the satellite trajectory)

8. Rotation matrices for rotations by x, y, z axes

For citation information please look at : 'Citation <https://cavsiopy.readthedocs.io/en/latest/citation.html>'_

Detailed information on the current and previous releases: 'Releases <https://cavsiopy.readthedocs.io/en/latest/releases.html>'_

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   getting-started
   citation
   releases
   installation
   ephemeris_importer
   use_rotation_matrices
   attitude_analysis
   attitude_plotter
   miscellaneous
