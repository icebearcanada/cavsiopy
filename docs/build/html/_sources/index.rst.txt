.. cavsiopy documentation master file, created by
   sphinx-quickstart on Wed Mar  8 10:32:43 2023.

Welcome to cavsiopy's documentation!
====================================
.. image:: https://github.com/icebearcanada/cavsiopy/tree/master/logos/cavsiopy.png

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

Version 1.1.0 Release Highlights: 'pysofa2' Integration and Updates
-------------------------------------------------------------------

In previous versions, cavsiopy used the 'pysofa' package developed by Frederic Grollier in 2010.

Starting from version 1.1.0, we have transitioned to using 'pysofa2.'

To address missing functions in 'pysofa2,' we introduced 'complement_missing_sofa.py,' which utilizes the SOFA C Library compiled by 'pysofa2.'

In addition, the following were updated:

1. In the 'rri\_example,' the 'rotate\_rri' function, which rotated the instrument's body vector, has been replaced with 'rotate\_inst.'
2. 'attitude\_3d\_ground\_quiver' has been enhanced to display a line connecting the subsatellite point with the ground target on the ground map.
3. name changes for several functions in auxiliary\preliminary\_data\_analysis.py module : plot\_data\_validity, import\_quaternions.

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   getting-started
   installation
   ephemeris_importer
   use_rotation_matrices
   attitude_analysis
   attitude_plotter
   miscellaneous
