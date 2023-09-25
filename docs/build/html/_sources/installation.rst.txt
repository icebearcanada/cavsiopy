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

- Please follow the installation guide of cartopy for other systems from here: https://scitools.org.uk/cartopy/docs/latest/installing.html

- installation of cartopy using conda:

.. code-block::

       $ conda install -c conda-forge cartopy

Proceed to installation
^^^^^^^^^^^^^^^^^^^^^^^
Download the requirements.txt from https://github.com/icebearcanada/cavsiopy

Using pip:

.. code-block::

       $ pip install -r requirements.txt

** Now, we are good to go! **

.. code-block::

       $ pip install cavsiopy


Credits: E. Ceren Kalafatoglu Eyiguler, Warren Holley, Andrew D. Howarth, Donald W. Danskin, Kuldeep Pandey, Carley J. Martin

Contributing: Glenn C. Hussey, Robert G. Gillies, Andrew W. Yau
