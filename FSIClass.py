#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
import scipy as Sci


class FSI(object):
    'Toolkit for Boundary Elelment Method Fluid Structure Interaction'
    def __init__(self, Npanels, Nelements):
        """
        Iniitalizes object related variables needed for other class methods.
        
        Keyword arguments:
        Npanels -- Number of fluid body panels
        Nelements -- Number of solid body elements
        """
        
        self.Npanels = Npanels
        self.Nelements = Nelements
        self.fluidNodeDispl = np.zeros((Npanels+1, 2))
        self.fluidNodeDisplOld = np.zeros((Npanels+1, 2))
        self.solidNodeDispl = np.zeros((Npanels+1, 2))
        self.nodeDispl = np.zeros((Nelements+1, 2))
        self.nodeDisplOld = np.zeros((Nelements+1, 2))
        self.fsiResidual = np.zeros((Npanels+1, 2))
        self.fsiResidualOld = np.zeros((Npanels+1, 2))
        self.nodeResidual = np.zeros((Nelements+1, 2))
        self.nodeResidualOld = np.zeros((Nelements+1, 2))
        self.initialFsiResidualNorm = 0
        self.maxInitialFsiResidualNorm = 0
        self.fsiResidualNorm = 0
        self.maxFsiResidualNorm = 0
        
    def rotatePts(x0, y0, theta):
        """Rotates a pair of points a specified angle"""
        x = x0 * np.cos(theta) - y0 * np.sin(theta)
        y = x0 * np.sin(theta) + y0 * np.cos(theta)
        return x, y
        
    def setInterfaceDisplacemet(self, displ, relaxationFactor, 
                                residual, outerCorr, couplingScheme):
        """
        Determines the relaxed solid/fluid body position based on a relaxation 
        factor and residual between solid and fluid domains.
        
        Keyword arguments:
        displ --
        relaxationFactor --
        residual --
        outerCorr --
        souplingScheme --
        """
        
        if (outerCorr < 3 or couplingScheme == 'FixedRelaxation'):
            # Use fixed-point relaxation
            self.fluidNodeDisplOld = self.fluidNodeDispl
            self.nodeDisplOld = self.nodeDispl
            self.fluidNodeDispl = self.fluidNodeDispl + self.fsiRelaxationFactor * self.fsiResidual
            self.nodeDispl = self.nodeDispl + self.fsiRelaxationFactor * self.nodeResidual
        elif (outerCorr >= 3 and couplingScheme == 'Aitken'):
            # Determine the new relaxation factor using the Aitken Method
            self.fsiRelaxationFactor = self.fsiRelaxationFactor * (np.transpose(self.fsiResidualOld) * (self.fsiResidualOld - self.fsiResidual)) / (Sci.linalg.norm(self.fsiResidualOld - self.fsiResidual))**2
            self.fsiRelaxationFactor = Sci.linalg.norm(self.fsiRelaxationFactor)
            if (self.fsiRelaxationFactor > 1.):
                self.fsiRelaxationFactor = 1
            
            self.fluidNodeDisplOld = self.fluidNodeDispl
            self.nodeDisplOld = self.nodeDispl
            self.fluidNodeDispl = self.fluidNodeDispl + self.fsiRelaxationFactor * self.fsiResidual
            self.nodeDispl = self.nodeDispl + self.fsiRelaxationFactor * self.nodeResidual
        else:
            # Not sure how to make this throw an exception and hault exec.
            print 'ERROR! Invalid coupling scheme "%s"' % couplingScheme
            print 'Valid coupling schemes are:'
            print '    "Fixed Relaxation"'
            print '    "Aitken"'
                   
    def readFsiControls(self, fixedPtRelax, nOuterCorrMax):
        """
        Initializes FSI relaxation coupling variables.
        
        Keyword arguments:
        fixedPtRelax -- Fixed-point relaxation value for Newton iteration
        nOutCorrMax -- Maximum allowed FSI coupling loops
        """
        self.fsiRelaxationFactorMin = fixedPtRelax
        self.fsiRelaxationFactor = self.fsiRelaxationFactorMin
        self.nOuterCorr = nOuterCorrMax
        
    def setInterFaceForce(self, outerCorr, nodes, nodesNew, theta, heave, 
                          x_b, z_b, xp, zp, xc, zc, P_b, ViscDrag, vn, delFs, 
                          interpMtd, meanline_c0, tBeamStruct, fixedCounter, c,
                          U_nPlus, Udot_nPlus, i_t):
        """
        Updates the structural mesh position, calculates the traction forces on
        the free nodes, and determines the initial condisitons for the timestep.
        Returns the node positions, traction forces, and initial conditions for
        the structural solver.        
        
        Keyword arguments:
        """
        # Superposing the structural displacements
        if (outerCorr > 1):
            nodes[:,0] = nodes[:,0] + self.nodeDispl[:,0] - self.nodeDisplOld[:,0]
            nodes[:,1] = nodes[:,1] + self.nodeDispl[:,1] - self.nodeDisplOld[:,1]
            
        if (outerCorr <= 1):
            # Updating the new kinematics
            nodes[:,0] = (nodesNew[:,0] - nodesNew[0,0])*np.cos(theta)
            nodes[:,1] = heave + (nodesNew[:,0] - nodesNew[0,0])*np.sin(theta)
            
            # Calculating the shift in node positions with the swimming velocity
            nodeDelxp = x_b[1]*np.ones((self.Nelements+1,1))
            nodeDelzp = z_b[1]*np.ones((self.Nelements+1,1))
            
            #Superposiing the kinematics and swimming translations
            nodes[:,0] = nodes[:,0] + nodeDelxp
            nodes[:,1] = nodes[:,1] + nodeDelzp
            
        # Determine the load conditons from the fluid solver
        # Calculate the panel lengths
        lp = np.zeros((self.Npanels,1))
        for i in xrange(self.Npanels):
            lp[i,0] = np.sqrt((xp[i+1] - xp[i])**2 + (zp[i+1] - zp[i])**2)
        
        # Calculate the force magnitude acting on the panel due to pressure,
        # then calculate the x-z components of this force
        magPF = P_b * lp * 1.
        pF = np.zeros((self.Npanels,2))
        if (ViscDrag == 1):
            pF[:,0] = (magPF * vn[:,0] * -1.) + delFs[:,0]
            pF[:,1] = (magPF * vn[:,1] * -1.) + delFs[:,1]
        else:
            pF[:,0] = magPF * vn[:,0] * -1.
            pF[:,1] = magPF * vn[:,1] * -1.
            
        # Determine the moment arm between top and bottom panel points, and
        # collapse force and moments to the camber line
        colM = np.zeros((0.5*self.Npanels,1))
        colPF = np.zeros((0.5*self.Npanels,2))
        meanPt = np.zeros((0.5*self.Npanels,2))
        for i in xrange(0.5*self.Npanels):
            meanPt[i,0] = 0.5*(xc[i] - xc[-i])
            meanPt[i,1] = 0.5*(zc[i] - zc[-i])
            colPF = pF[i,:] + pF[-i,:]
            colM = -1. * pF[i,0] * (zc[i,0] - meanPt[i,1]) + \
                   pF[i,1] * (xc[i,0] - meanPt[i,0]) + \
                   -1. * pF[-i,0] * (zc[-i] - meanPt[i,1]) + \
                   pF[-i,1] * (xc[-i] - meanPt[i,0])
                   
        colPF = np.flipud(colPF)
        colM = np.flipud(colM)
        relXp = xp - np.min(xp)
        relXc = xc - np.min(xp)
        
        # Interpolate the collapsed forces and moments onto the structural mesh
        nodalInput = np.zeros((self.Nelements,5))
        if (interpMtd == 1):
            nodalInput[:,0] = np.interp(nodes[:,2], meanline_c0[0.5*self.Npanels+1:self.Npanels,0], colPF[:,0], left=True, right=True)
            nodalInput[:,1] = np.interp(nodes[:,2], meanline_c0[0.5*self.Npanels+1:self.Npanels,0], colPF[:,1], left=True, right=True)
            nodalInput[:,5] = np.interp(nodes[:,2], meanline_c0[0.5*self.Npanels+1:self.Npanels,0], colM[:,0], left=True, right=True)
        else:
            nodalInput[:,0] = Sci.interpolate.spline(meanline_c0[0.5*self.Npanels+1:self.Npanels,0], colPF[:,0], nodes[:,2])
            nodalInput[:,1] = Sci.interpolate.spline(meanline_c0[0.5*self.Npanels+1:self.Npanels,0], colPF[:,1], nodes[:,2])
            nodalInput[:,5] = Sci.interpolate.spline(meanline_c0[0.5*self.Npanels+1:self.Npanels,0], colM[:,0], nodes[:,2])

        # Rotate force components into the relative cooridnate system
        nodalInput[:,0], nodalInput[:,1] = self.rotatePts(nodalInput[:,0], nodalInput[:,1], theta)
        
        # Create the load matrix
        Fload = np.zeros((3*(self.Nelements+1),1))
        Fload[0:-3:3,0] = nodalInput[:,0]
        Fload[1:-2:3,0] = nodalInput[:,1]
        Fload[2:-1:3,0] = nodalInput[:,5]
        
        # Create element area matrix
        A = tBeamStruct[:,0]
        
        # Create area moment of inertia matrix
        I = 1. * tBeamStruct[:,0]**3 / 12
        
        # Set fixed nodes
        fixedNodes = fixedCounter
        
        # Initial element length
        l_0 = c / self.Nelements
        
        # Initial displacements and velocities
        if (i_t <= 2 and outerCorr <= 1):
            U_n = np.zeros((3*(self.Nelements+1),1))
            Udot_n = np.zeros((3*(self.Nelements+1),1))
            #UdotDot_n = np.zeros((3*(self.Nelements+1),1))
            #U_nPlus = np.zeros((3*(self.Nelements+1),1))
            #Udot_nPlus = np.zeros((3*(self.Nelements+1),1))
            #UdotDot_nPlus = np.zeros((3*(self.Nelements+1),1))
        if (i_t > 1 and outerCorr <= 1):
            U_n = np.zeros((3*(self.Nelements+1),1))
            Udot_n = np.zeros((3*(self.Nelements+1),1))
            U_n[3*fixedNodes+1:-1,0] = U_nPlus
            Udot_n[3*fixedNodes+1:-1,0] = Udot_nPlus
            
        return relXp, relXc, Fload, A, I, l_0, U_n, Udot_n
            
    def getDisplacements(self, theta, heavePos, xp, zp, nodalDelxp, nodalDelzp, tBeam, nodes_0, U_nPlus, interpMtd, meanline_p0, fixedNodes, flexionRatio):
        """
        Returns the calculated displacements of the fluid body based on the 
        displacements calculated by the structural body. This is used to 
        calculate the FSI coupling residual used in an interative strong 
        coupling between the fluid and solid solutions.
        
        Keyword arguments:
        theta -- 
        heavePos -- 
        xp -- 
        zp -- 
        nodalDelxp -- 
        nodalDelzp -- 
        tBeam -- 
        nodes_0 -- 
        U_nPlus --
        interpMtd -- 
        meanline_p0 --
        fixedNodes -- 
        flexionRatio --
        """        
        # Get the absolute x and z displacements
        nodeDisplacements = np.zeros((self.Nelements+1,2))
        nodeDisplacements[:,1] =  U_nPlus[2:-2:3,0]
        nodeDisplacements[:,0] =  (nodes_0[fixedNodes+1:-1,1] + \
                                  nodeDisplacements[:,1]) * np.sin(-U_nPlus[3:-1:3,0] + \
                                  U_nPlus[1:-3:3,0])
        
        # Calculate the new structural locations
        tempNodes = nodes_0;
        tempNodes[fixedNodes+1:-1,0] = nodes_0[fixedNodes+1:-1,0] + nodeDisplacements[:,0] # New x-position
        tempNodes[fixedNodes+1:-1,1] = nodes_0[fixedNodes+1:-1,1] + nodeDisplacements[:,1] # New z-position
        
        #frameNodes = tempNodes
        
        tempNodes[:,0], tempNodes[:,1] = self.rotatePts(tempNodes[:,0], tempNodes[:,1], theta)
        
        tempNodes[:,1] = tempNodes[:,1] + heavePos
        
        tempNodes[:,0] = tempNodes[:,0] + nodalDelxp
        tempNodes[:,1] = tempNodes[:,1] + nodalDelzp;
              
        normal = np.zeros((self.Nelements,2))
        for i in xrange(self.Nelemenets):
            normal[i,0] = -1 * (tempNodes[i+1,1] - tempNodes[i,1])
            normal[i,1] =  1 * (tempNodes[i+1,0] - tempNodes[i,0])
            normal[i,:] = normal[i,:] / np.sqrt(normal[i,0]**2 + normal[i,1]**2)
        
        topNodes = np.zeros((self.Nelements+1, 3))
        bottomNodes = np.zeros((self.Nelements+1, 3))        
        topNodes[:,0] = tempNodes[0:-2,0] + 0.5 * tBeam[:,0] * normal[:,0]
        topNodes[:,1] = tempNodes[0:-2,1] + 0.5 * tBeam[:,1] * normal[:,1]
        topNodes[-1,0:1] = tempNodes[-1,0:1]
        topNodes[:,2] = tempNodes[:,2]
        bottomNodes[:,1] = tempNodes[0:-2,0] - 0.5 * tBeam[:,0] * normal[:,0]
        bottomNodes[:,2] = tempNodes[0:-2,1] - 0.5 * tBeam[:,0] * normal[:,1]
        bottomNodes[-1,0:1] = tempNodes[-1,0:1]
        bottomNodes[:,2] = tempNodes[:,2]
        
        # Interpolate the structual results back to the fluid domain
        i = np.rint(0.5 * np.shape(meanline_p0)[0])
        bottomXp = np.zeros_like(meanline_p0)
        bottomZp = np.zeros_like(meanline_p0)
        topXp = np.zeros_like(meanline_p0)
        topZp = np.zeros_like(meanline_p0)
        if (interpMtd == 1):
            bottomXp[:,0] = np.interp(meanline_p0[0:i-1,0], bottomNodes[:,2], bottomNodes[:,0], left=True, right=True)
            bottomZp[:,0] = np.interp(meanline_p0[0:i-1,0], bottomNodes[:,2], bottomNodes[:,1], left=True, right=True)
            topXp[:,0] = np.interp(meanline_p0[i:-1,0], topNodes[:,2], topNodes[:,0], left=True, right=True)
            topZp[:,0] = np.interp(meanline_p0[i:-1,0], topNodes[:,2], topNodes[:,1], left=True, right=True)
        else:
            bottomXp[:,0] = Sci.interpolate.spline(bottomNodes[:,2], bottomNodes[:,0], meanline_p0[0:i-1,0])      
            bottomZp[:,0] = Sci.interpolate.spline(bottomNodes[:,2], bottomNodes[:,1], meanline_p0[0:i-1,0])
            topXp[:,0] = Sci.interpolate.spline(topNodes[:,2], topNodes[:,0], meanline_p0[i:-1,0])
            topZp[:,0] = Sci.interpolate.spline(topNodes[:,2], topNodes[:,1], meanline_p0[i:-1,0])
        
        newxp = np.zeros_like(meanline_p0)
        newzp = np.zeros_like(meanline_p0)
        newxp[0:i-1,1] = bottomXp
        newxp[i:-1,1] = topXp
        newzp[0:i-1,1] = bottomZp
        newzp[i:-1,1] = topZp
        
        for i in xrange(self.Npanels+1):
            if (meanline_p0[i,0] <= flexionRatio):
                newxp[i,0] = xp[i,0]
                newzp[i,0] = zp[i,0]

        return newxp - xp, newzp - zp, tempNodes
        
    def calcFSIResidual(self, DU, nodes, tempNodes, outerCorr):
        self.solidNodeDispl = DU
        self.fsiResidualOld = self.fsiResidual
        self.nodeResidualOld = self.nodeResidual
        self.fsiResidual = self.solidNodeDispl - self.fluidNodeDispl
        self.nodeResidual = (tempNodes[:,0:1]-nodes[:,0:1]) - self.nodeDispl
        magFsiResidual = np.sqrt(self.fsiResidual[:,0]**2 + self.fsiResidual[:,1]**2)
        np.linalg.norm()
        self.fsiResidualNorm = np.linalg.norm(self.fsiResidual, ord=2)
        self.maxFsiResidualNorm = np.linalg.norm(self.fsiResidual, ord='inf')
        
        if (outerCorr == 1 ):
            self.initialFsiResidualNorm = self.fsiResidualNorm
            self.maxInitialFsiResidualNorm = self.maxFsiResidualNorm

        self.fsiResidualNorm = self.fsiResidualNorm / self.initialFsiResidualNorm
        self.maxFsiResidualNorm = self.maxFsiResidualNorm / self.maxInitialFsiResidualNorm
        
        return magFsiResidual, self.fsiResidualNorm, self.maxFsiResidualNorm

            