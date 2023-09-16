# -*- coding: utf-8 -*-

#
# Copyright 2010 Frédéric Grollier
#
# Distributed under the terms of the MIT license
#

"""
This function complements the undeclared c2t00b in the 'pysofa2' module.
It locates the SOFA C library used by 'pysofa2' and uses it to construct
a celestial-to-terrestrial rotation matrix using SOFA routines.

The 'c2t00b' function is adapted from 'pysofa' by Frédéric Grollier (2010).
"""

from ctypes import  c_double
from numpy.ctypeslib import ndpointer
from numpy import zeros, asmatrix
from pysofa2 import sofa_wrapper


_sofa = sofa_wrapper._sofa

# iauC2t00b
_sofa.iauC2t00b.argtypes = [c_double, #tta
                            c_double, #ttb
                            c_double, #uta
                            c_double, #utb
                            c_double, #xp
                            c_double, #yp
                            ndpointer(shape=(3,3), dtype=float)]

def c2t00b(tta, ttb, uta, utb, xp, yp):
    """ Form the celestial-to-terrestrial matrix given the date, the UT1 and
    the polar motion, using IAU 2000B nutation model.

    :param tta, ttb: TT as a two-part Julian date.
    :type tta, ttb: float

    :param uta, utb: UT1 as a two-part Julian date.
    :type uta, utb: float

    :param xp, yp: coordinates of the pole in radians.
    :type xp, yp: float

    :returns: the celestial-to-terrestrial matrix, as a numpy.matrix of shape
        3x3.

    .. seealso:: |MANUAL| page 42
    """
    rc2t = asmatrix(zeros(shape=(3,3), dtype=float))
    _sofa.iauC2t00b(tta, ttb, uta, utb, float(xp), float(yp), rc2t)
    return rc2t

