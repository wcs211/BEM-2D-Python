# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:59:14 2015

@author: biofluids_1
"""

import numpy as np

class SwimmerParameters(object):
    def __init__(self, CE, SW_GEOMETRY, SW_KUTTA):
        self.CE = CE
        self.SW_GEOMETRY = SW_GEOMETRY
        self.SW_KUTTA = SW_KUTTA

class GeoVDVParameters(object):
    def __init__(self, N, S, C, K, EPSILON):
        self.N = N
        self.S = S
        self.C = C
        self.K = K
        self.EPSILON = EPSILON
    
class MotionParameters(object):
    def __init__(self, V0, THETA_MAX, H_C, F, PHI):
        self.V0 = V0
        self.THETA_MAX = THETA_MAX
        self.H_C = H_C
        self.F = F
        self.PHI = PHI
    
class BodyBFC(object):
    def __init__(self, x, z, x_col, z_col):
        self.x = x
        self.z = z
        self.x_col = x_col
        self.z_col = z_col
    
class BodyAFC(object):
    def __init__(self, N):
        self.x = np.empty(N+1)
        self.z = np.empty(N+1)
        self.x_col = np.empty(N)
        self.z_col = np.empty(N)
        self.x_mid = np.zeros((3,N))
        self.z_mid = np.zeros((3,N))
        self.x_neut = np.empty(N/2)
        self.z_neut = np.empty(N/2)
        self.x_le = 0.
        self.z_le = 0.