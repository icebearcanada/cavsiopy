# cavsiopy: Calculation and Visualization of Spacecraft Instrument Orientation

![cavsiopy](logos/cavsiopy.png)

Welcome to cavsiopy!
====================================
'cavsiopy' is an open source Python package, which was specifically developed to determine the look direction of the Radio Receiver Instrument (RRI) on the e-POP/ CASSIOPE/ Swarm-E for HF radio wave propagation studies, but can be applied to any satellite mission requiring accurate pointing information. Given the spacecraft position, the roll, yaw, and pitch angles, the body vector of the instrument, and a target location, the pointing direction and orientation of the instrument onboard the spacecraft with respect to the target can be determined in various reference frames. 'cavsiopy' supports the reference frames: Geocentric Equatorial Inertial J2000 (GEI J2K), Earth Centered Earth Fixed (ECEF), International Terrestrial Reference Frame (ITRF), geodetic North-East-Down, and geocentric North-East-Center. 

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

# Citation
----------

If you use any of the cavsiopy functions/modules/plots please cite the software using DOI: https://doi.org/10.5281/zenodo.8361256 as in:

Kalafatoglu Eyiguler, E. C., Holley, W., Howarth, A. D., Danskin, D. W., Pandey, K., Martin, C. J., et al. (2023b). icebearcanada/cavsiopy: v1.1.1. doi:10.5281/zenodo.8361256 

cavsiopy in-action:

under review: 

E. Ceren Kalafatoglu Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley J. Martin, Robert G. Gillies, Andrew W. Yau, Glenn C. Hussey (2023). cavsiopy: A Python package to calculate and visualize spacecraft instrument orientation. Manuscript submitted for publication.

published:
1. Kalafatoglu Eyiguler, E. C., Danskin, D., Hussey, G., Pandey, K., Gillies, R., and Yau, A. (2022). Satellite attitude effects on the reception of transionospheric hf signals: Examples from the radio receiver instrument onboard e-pop/swarm-e. In 2022 3rd URSI Atlantic and Asia Pacific Radio Science Meeting(AT-AP-RASC) (IEEE), 1–4

2. Kalafatoglu Eyiguler, E. C., Danskin, D. W., Howarth, A. D., Holley, W., Pandey, K., Gillies, R. G., et al. (2023). Attitude effects on the observed orientation angle of hf waves from the radio receiver instrument on e-pop/swarm-e. In Ionospheric Effects Symposium 2023 (IES2023), ed. K. Groves (Alexandria, VA, USA: IES 2023), 1–6

3. Pandey, K., Eyiguler, E. K., Gillies, R. G., Hussey, G. C., Danskin, D. W., and Yau, A. W. (2022). Polarization characteristics of a single mode radio wave traversing through the ionosphere: A unique observation from the rri on epop/swarm-e. Journal of Geophysical Research: Space Physics 127. doi:10.1029/2022JA030684

Documentation
-------------
Detailed documentation on [readthedocs.org](http://cavsiopy.readthedocs.io/)

Dependencies
------------
numpy, matplotlib, astropy, cartopy, geopy, h5py, pysofa2

Installation using pip:
-----------------------
Before installing cavsiopy some of the dependencies may need special installation instructions:

cartopy
-------

- If you are installing cartopy via pip (on Ubuntu/Debian) first libgeos-dev has to be on your system:

```
$ sudo apt -y install libgeos-dev
```

- Please follow the installation guide of cartopy for other systems from here: https://scitools.org.uk/cartopy/docs/latest/installing.html

- installation of cartopy using conda:

```
$ conda install -c conda-forge cartopy
```

Other packages and requirements
-------------------------------
Download the requirements.txt from our [GitHub repository] (https://github.com/icebearcanada/cavsiopy)

```
$ pip install -r requirements.txt
```
** Now, we are good to go! **

```
$ pip install cavsiopy
```

Version 1.2.2
-------------

Patch release:

1. slew_example.py: added to examples

2. in attitude\_analysis.py: functions in utils.py and complement\_missing\_sofa.py are now embedded in attitude\_analysis.py

This structural change does not affect any end-user experience / call to functions etc.

3. attitude\_analysis.spacecraft\_distance\_from\_a\_point: fixed a minor bug, which caused the distance array to return empty.

4. use\_rotation\_matrices.py: utils.py and missing\_complement\_sofa.py imports removed

Version 1.2.1
-------------
Minor release: change in the return parameters of find\_attitude and cas\_ephemeris. 

1. attitude\_analysis.find\_attitude: East-North-Up: enu and Easth-North-Center(up): enc_u have been added to the returns.

2. utils.py: added to find the location of coefficient files for pysofa, and removed the declaration of coefficient files for pysofa for the find\_attitude function.

3. ephemeris\_importer.cas\_ephemeris: Eclipse parameter can now be imported using ephemeris\_importer.cas\_ephemeris function from the CAS_ephemeris data files.

4. \_\_init\_\_.py: Citation information has been added.

5. requirements.txt: file has been updated.

Version 1.1.1
-------------
patch release: 

1. 'attitude\_plotter.attitude\_3d\_ground\_quiver': has been enhanced to display a line connecting the subsatellite point with the ground target on the ground map.

2. rri\_example.py: rotate\_rri is renamed as 'rotate\_inst.'

3. auxiliary\preliminary\_data\_analysis.py: name changes for plot\_data\_validity, import\_quaternions.

Version 1.1.0
-------------
'pysofa2' integration.

In previous versions, cavsiopy used the 'pysofa' package developed by Frederic Grollier from 2010.

Starting from version 1.1.0, we have transitioned to using 'pysofa2.'

To address missing functions in 'pysofa2,' we introduced 'complement\_missing\_sofa.py,' which utilizes the SOFA C Library compiled by 'pysofa2.'

Version 1.0.0
-------------
Initial major release: Uses pysofa (https://pypi.org/project/pysofa/) to compute the ICRF to ITRF rotation matrix

Before v1.0.0
--------------
Test releases

Credits: E. Ceren Kalafatoglu Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley J. Martin

Contributing: Glenn C. Hussey, Robert G. Gillies, Andrew W. Yau
