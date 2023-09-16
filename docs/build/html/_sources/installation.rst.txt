Installation
============

Dependencies
------------
numpy, matplotlib, astropy, cartopy, geopy, h5py, pysofa2

Installation
------------
Before installing cavsiopy some of the dependencies need special installation instructions:

cartopy
^^^^^^^
- If you are installing cartopy via pip (on Ubuntu/Debian) first libgeos-dev has to be on your system:

.. code-block::

       $ sudo apt -y install libgeos-dev
       $ pip install -r requirements.txt

- Follow the installation guide for other systems: https://scitools.org.uk/cartopy/docs/latest/installing.html

- for installation of cartopy with conda:

.. code-block::

       $ conda install -c conda-forge cartopy

Other packages and requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block::

	pip install --no-deps astropy
	
	pip install pysofa2

	pip install geopy

	pip install -r requirements.txt

** Now, we are good to go! **

.. code-block::

	pip install cavsiopy


Credits: C. Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley Martin

Contributing: Glenn C. Hussey, Robert Gillies, Andrew W. Yau

License: GNU LESSER GENERAL PUBLIC LICENSE Version 3, 29 June 2007
