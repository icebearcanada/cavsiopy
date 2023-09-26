Releases
========

Current release: v1.2.2

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
