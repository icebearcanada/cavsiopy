#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 15:20:54 2023

@author: ceren
"""

import os
import pkg_resources

def determine_installation_directory(package_name):
    try:
        distribution = pkg_resources.get_distribution(package_name)
        installation_dir = distribution.location
        return installation_dir
    except pkg_resources.DistributionNotFound:
        # Handle the case where the package is not found
        print(package_name, 'not found!')
        return None

def get_txt_file_path(filename):
    installation_dir = determine_installation_directory('cavsiopy')  # get the installation directory
    if installation_dir:
        txt_file_path = os.path.join(installation_dir, 
                                     "cavsiopy/coefficient_files/",
                                     filename)  # Adjust the path as needed
        return txt_file_path
    else:
        return None
    