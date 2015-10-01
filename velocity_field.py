#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
from functions_general import panel_vectors, transformation

def induced_velocity(Swimmers, i):
    NT = i # Number of targets (wake panel points that are rolling up)
    for SwimT in Swimmers:
        SwimT.Wake.vx = np.zeros(NT)
        SwimT.Wake.vz = np.zeros(NT)
        DELTA_CORE = SwimT.DELTA_CORE
        for SwimI in Swimmers:
            # Coordinate transformation for body panels influencing wake
            (xp1, xp2, zp) = transformation(SwimT.Wake.x[1:i+1], SwimT.Wake.z[1:i+1], SwimI.Body.AF.x, SwimI.Body.AF.z)
    
            # Angle of normal vector with respect to global z-axis
            (nx, nz) = panel_vectors(SwimI.Body.AF.x, SwimI.Body.AF.z)[2:4]
            beta = np.arctan2(-nx, nz)
    
            # Katz-Plotkin eqns 10.20 and 10.21 for body source influence
            dummy1 = np.log((xp1**2+zp**2)/(xp2**2+zp**2))/(4*np.pi)
            dummy2 = (np.arctan2(zp,xp2)-np.arctan2(zp,xp1))/(2*np.pi)
    
            # Rotate back to global coordinates
            dummy3 = dummy1*np.cos(beta) - dummy2*np.sin(beta)
            dummy4 = dummy1*np.sin(beta) + dummy2*np.cos(beta)
    
            # Finish eqns 10.20 and 10.21 for induced velocity by multiplying with sigma
            SwimT.Wake.vx += np.dot(dummy3, SwimI.Body.sigma)
            SwimT.Wake.vz += np.dot(dummy4, SwimI.Body.sigma)
    
            # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
            NI = SwimI.Body.N+1
            xp = np.repeat(SwimT.Wake.x[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Body.AF.x[:,np.newaxis].T, NT, 0)
            zp = np.repeat(SwimT.Wake.z[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Body.AF.z[:,np.newaxis].T, NT, 0)
    
            # Find distance r_b between each influence/target
            r_b = np.sqrt(xp**2+zp**2)
    
            # Katz-Plotkin eqns 10.9 and 10.10 for body doublet (represented as point vortices) influence
            dummy1 = zp/(2*np.pi*(r_b**2+DELTA_CORE**2))
            dummy2 = -xp/(2*np.pi*(r_b**2+DELTA_CORE**2))
    
            # Finish eqns 10.9 and 10.10 by multiplying with Body.gamma, add to induced velocity
            SwimT.Wake.vx += np.dot(dummy1, SwimI.Body.gamma)
            SwimT.Wake.vz += np.dot(dummy2, SwimI.Body.gamma)
    
            # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
            NI = SwimI.Edge.N+1
            xp = np.repeat(SwimT.Wake.x[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Edge.x[:,np.newaxis].T, NT, 0)
            zp = np.repeat(SwimT.Wake.z[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Edge.z[:,np.newaxis].T, NT, 0)
    
            # Find distance r_e between each influence/target
            r_e = np.sqrt(xp**2+zp**2)
    
            # Katz-Plotkin eqns 10.9 and 10.10 for edge (as point vortices) influence
            dummy1 = zp/(2*np.pi*(r_e**2+DELTA_CORE**2))
            dummy2 = -xp/(2*np.pi*(r_e**2+DELTA_CORE**2))
    
            # Finish eqns 10.9 and 10.10 by multiplying with Edge.gamma, add to induced velocity
            SwimT.Wake.vx += np.dot(dummy1, SwimI.Edge.gamma)
            SwimT.Wake.vz += np.dot(dummy2, SwimI.Edge.gamma)
    
            # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
            NI = i+1
            xp = np.repeat(SwimT.Wake.x[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Wake.x[:i+1,np.newaxis].T, NT, 0)
            zp = np.repeat(SwimT.Wake.z[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Wake.z[:i+1,np.newaxis].T, NT, 0)
    
            # Find distance r_w between each influence/target
            r_w = np.sqrt(xp**2+zp**2)
    
            # Katz-Plotkin eqns 10.9 and 10.10 for wake (as point vortices) influence
            dummy1 = zp/(2*np.pi*(r_w**2+DELTA_CORE**2))
            dummy2 = -xp/(2*np.pi*(r_w**2+DELTA_CORE**2))
    
            # Finish eqns 10.9 and 10.10 by multiplying with Wake.gamma, add to induced velocity
            SwimT.Wake.vx += np.dot(dummy1, SwimI.Wake.gamma[:i+1])
            SwimT.Wake.vz += np.dot(dummy2, SwimI.Wake.gamma[:i+1])