#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script to rotate instrument pointing direction in orbital frame.

This specific example uses RRI data and rotates the instrument using RPY.

Users can change the ordering to their liking.
"""
import numpy as np
import cavsiopy.use_rotation_matrices as rm
import cavsiopy.ephemeris_importer as ei

file_RRI = "RRI_20160418_222759_223156_lv1_12.0.0.h5" 

body_vec = [1, 0, 0]

dict_rri  = ei.rri_ephemeris(file_RRI)

# converting degrees to radians
rr = np.deg2rad(dict_rri['roll'])
pr = np.deg2rad(dict_rri['pitch'])
yr = np.deg2rad(dict_rri['yaw'])

# ORF (Local vertical local horizontal) 3 axis unit vectors
orf = np.array([body_vec,] *len(rr)).T

# Initialize body Vectors
rot_vec = np.empty(orf.shape)

# transform initial definition of ORF to rotated frame
for i in range(rr.shape[0]):

    # direction cosine matrice
    R = rm.RX_r2i(rr[i]) @ rm.RY_r2i(pr[i]) @ rm.RZ_r2i(yr[i])

    # rotation
    rot_vec[:, i] = R @ orf[:, i]

#  rotated instrument look direction in ORF
rotated_body=np.array(rot_vec)
