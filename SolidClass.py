#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np

class solid(object):
    'Toolkit for Finite Element structural analysis'
    def __init__(self, Nnodes,xp_0,zp_0,tmax):
        """Iniitalizes object related variables needed for other class methods."""
        self.Nnodes = Nnodes
        self.Nelements = Nnodes - 1
        self.xp_0 = xp_0
        self.zp_0 = zp_0
        self.pivotPoint = (0.5*tmax) / (max(xp_0) - min(xp_0))
        self.nodes = np.zeros((Nnodes,3))
        self.nodes_0 = self.nodes
        self.tBeam = np.zeros((self.Nelements,1))
        self.ttemp = self.tBeam
        self.tBeamStruct = self.tBeam
        
    def initThinPlate(self,tmax,c,constThickBeam,tConst,flexionRatio):
        """
        This function initializes the element nodal positions.
        
        Keyword arguments:
        tmax -- array of solid thicknesses for each element
        c -- undeformed/initial chord length
        constThickBeam -- flag argument for constant thickness properties 
        tConst -- constant thickness position
        flexionRatio -- fraction of rigid body
        """
        beamCounter = 0
        fixedCounter = 0
        for i in xrange(self.Nelements):
            if self.nodes[i,0] <= 0.5*tmax:
                self.tBeam[i,0] = tmax
                self.tBeamStruct[i,0] = self.tBeam[i,0]
                self.ttemp[i,0] = self.tBeam[i,0]
                beamCounter += 1
                fixedCounter += 1
            elif self.nodes[i,0] >= c-0.5*tmax:
                self.tBeam[i,0] = 2*np.sqrt((0.5*tmax)**2 -(self.nodes[i,0]-(c-0.5*tmax))**2 )
                self.tBeamStruct[i,0] = self.tBeam[i,0]
                self.ttemp[i,0] = self.tBeam[i,0]
            else:
                self.ttemp[i,0] = tmax
                self.tBeam[i,0] = tmax
                self.tBeamStruct[i,0] = self.tBeam[i,0]
                if ( constThickBeam == 1 and  self.nodes[i,2] >= tConst):
                    self.tBeamStruct[i,0] = self.tBeamStruct[i-1,0]
                else:
                    self.tBeamStruct[i,0] = self.tBeam[i,0]
                if (self.nodes[i,2] <= flexionRatio):
                    fixedCounter += 1

    def meanline(self, xc_0, xp_0):
        """
        Function to calculate the meanline fraction position along the body.
        0 corresponds to the leading edge 
        1 corresponds to the trailing edge
        
        Keyword arguments:
        xc_0 -- Initial colocation point (x-component)
        xp_0 -- Iniital element endpoint (x-component)
        """
        
        meanline_p0 = xp_0 / (np.max(xp_0) - np.min(xp_0))
        meanline_c0 = xc_0 / (np.max(xp_0) - np.min(xp_0))
        
        return meanline_p0, meanline_c0

    def initMesh(self):
        """
        Initializes the finite element mesh based on the object's __init__
        values. This is only valid for the undeformed structure at time t = 0.
        """
        self.nodes[:,0] = np.arange(min(self.xp_0),max(self.xp_0),\
                                (max(self.xp_0)-min(self.xp_0))/self.Nelements)
        self.nodes[:,1] = np.zeros((self.Nnodes,1))
        self.nodes[:,2] = self.nodes[:,0] / (max(self.nodes[:,0])-min(self.nodes[:,1]))
        self.nodes_0 = self.nodes
        
    def rotatePts(x0, y0, theta):
        x = x0 * np.cos(theta) - y0 * np.sin(theta)
        y = x0 * np.sin(theta) + y0 * np.cos(theta)
        return x, y
        