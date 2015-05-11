#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
import scipy as Sci
from functions_general import panel_vectors


class FSI(object):
    'Toolkit for Boundary Elelment Method Fluid Structure Interaction'
    def __init__(self, Body, Solid):
        """
        Iniitalizes object related variables needed for other class methods.
        
        Keyword arguments:
        Npanels -- Number of fluid body panels
        Nelements -- Number of solid body elements
        """
        self.fluidNodeDispl = np.zeros((Body.N+1, 2))
        self.fluidNodeDisplOld = np.zeros((Body.N+1, 2))
        self.solidNodeDispl = np.zeros((Body.N+1, 2))
        self.nodeDispl = np.zeros((Solid.Nelements+1, 2))
        self.nodeDisplOld = np.zeros((Solid.Nelements+1, 2))
        self.fsiResidual = np.zeros((Body.N+1, 2))
        self.fsiResidualOld = np.zeros((Body.N+1, 2))
        self.nodeResidual = np.zeros((Solid.Nelements+1, 2))
        self.nodeResidualOld = np.zeros((Solid.Nelements+1, 2))
        self.initialFsiResidualNorm = 0
        self.maxInitialFsiResidualNorm = 0
        self.fsiResidualNorm = 0
        self.maxFsiResidualNorm = 0
        self.maxMagFsiResidual = 0
        self.DU = np.zeros((Body.N+1,2))
        self.maxDU = 0
     
    def rotatePts(self, x0, y0, theta):
        """Rotates a pair of points a specified angle"""
        x = x0 * np.cos(theta) - y0 * np.sin(theta)
        y = x0 * np.sin(theta) + y0 * np.cos(theta)
        return (x, y)
        
    def setInterfaceDisplacemet(self, outerCorr, couplingScheme):
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
#            print self.nodeResidual[-5:,:]
            self.nodeDispl = self.nodeDispl + self.fsiRelaxationFactor * self.nodeResidual
        elif (outerCorr >= 3 and couplingScheme == 'Aitken'):
            # Determine the new relaxation factor using the Aitken Method
            self.fsiRelaxationFactor = self.fsiRelaxationFactor * (self.fsiResidualOld * (self.fsiResidualOld - self.fsiResidual)) / (np.linalg.norm(self.fsiResidualOld - self.fsiResidual))**2
            self.fsiRelaxationFactor = np.linalg.norm(self.fsiRelaxationFactor)
            if (self.fsiRelaxationFactor > 1.):
                self.fsiRelaxationFactor = 1
#            print self.nodeResidual[-5:,:]
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
        
    def setInterfaceForce(self, Solid, Body, PyFEA, t, TSTEP, outerCorr, 
                          SWITCH_VISC_DRAG, delFs, SWITCH_INTERP_MTD, C, i_t):
        """
        Updates the structural mesh position, calculates the traction forces on
        the free nodes, and determines the initial condisitons for the timestep.
        Returns the node positions, traction forces, and initial conditions for
        the structural solver.        
        
        Keyword arguments:
        """
        # Determine current pitching angle and heave position
        # TODO: Add heaving functionality to kinematics
        theta = Body.THETA_MAX * np.sin(2 * np.pi * Body.F * (t + TSTEP) + Body.PHI)  
        heave = 0
        
        # Superposing the structural displacements
        if (outerCorr > 1):
            Solid.nodes[:,0] = Solid.nodes[:,0] + self.nodeDispl[:,0] - self.nodeDisplOld[:,0]
            Solid.nodes[:,1] = Solid.nodes[:,1] + self.nodeDispl[:,1] - self.nodeDisplOld[:,1]          
            
        if (outerCorr <= 1):
            # Updating the new kinematics
            Solid.nodes[:,0] = (Solid.nodesNew[:,0] - Solid.nodesNew[0,0])*np.cos(theta)
            Solid.nodes[:,1] = heave + (Solid.nodesNew[:,0] - Solid.nodesNew[0,0])*np.sin(theta)
            
            # Calculating the shift in node positions with the swimming velocity
            nodeDelxp = Body.xLE[1] * np.ones((Solid.Nelements + 1,1))
            nodeDelzp = Body.zLE[1] * np.ones((Solid.Nelements + 1,1))
            
            #Superposiing the kinematics and swimming translations
            Solid.nodes[:,0] = Solid.nodes[:,0] + nodeDelxp.T
            Solid.nodes[:,1] = Solid.nodes[:,1] + nodeDelzp.T
            
        # Determine the load conditons from the fluid solver
        # Calculate the panel lengths and normal vectors
        (nx,nz,lp) = panel_vectors(Body.x,Body.z)[2:5]
        
        # Calculate the force magnitude acting on the panel due to pressure,
        # then calculate the x-z components of this force
        magPF = Body.p * lp * 1.
        pF = np.zeros((Body.N,2))
        if (SWITCH_VISC_DRAG == 1):
            pF[:,0] = (magPF * nx * -1.) + delFs[:,0]
            pF[:,1] = (magPF * nz * -1.) + delFs[:,1]
        else:
            pF[:,0] = magPF * nx * -1.
            pF[:,1] = magPF * nz * -1.
            
        # Determine the moment arm between top and bottom panel points, and
        # collapse force and moments to the camber line
        colM = np.zeros((0.5*Body.N,1))
        colPF = np.zeros((0.5*Body.N,2))
        meanPt = np.zeros((0.5*Body.N,2))
        for i in xrange(int(0.5*Body.N)):
            meanPt[i,0] = 0.5*(Body.x_mid[0,i] - Body.x_mid[0,-i])
            meanPt[i,1] = 0.5*(Body.z_mid[0,i] - Body.z_mid[0,-i])
            colPF[i,:] = pF[i,:] + pF[-i,:]
            colM[i,:] = -1. * pF[i,0] * (Body.z_mid[0,i] - meanPt[i,1]) + \
                   pF[i,1] * (Body.x_mid[0,i] - meanPt[i,0]) + \
                   -1. * pF[-i,0] * (Body.z_mid[0,-i] - meanPt[i,1]) + \
                   pF[-i,1] * (Body.x_mid[0,-i] - meanPt[i,0])
                   
        colPF = np.flipud(colPF)
        colM = np.flipud(colM)
        relXp = Body.x - np.min(Body.x)
        relXc = Body.x_mid - np.min(Body.x)
        
        # Interpolate the collapsed forces and moments onto the structural mesh
        nodalInput = np.zeros((Solid.Nnodes,6))
        if (SWITCH_INTERP_MTD == 1):
            nodalInput[:,0] = np.interp(Solid.nodes[:,2], Solid.meanline_c0[0.5*Body.N:Body.N+1], colPF[:,0], left=True, right=True)
            nodalInput[:,1] = np.interp(Solid.nodes[:,2], Solid.meanline_c0[0.5*Body.N:Body.N+1], colPF[:,1], left=True, right=True)
            nodalInput[:,5] = np.interp(Solid.nodes[:,2], Solid.meanline_c0[0.5*Body.N:Body.N+1], colM[:,0], left=True, right=True)
        else:
            nodalInput[:,0] = Sci.interpolate.spline(Solid.meanline_c0[0.5*Body.N:Body.N+1], colPF[:,0], Solid.nodes[:,2])
            nodalInput[:,1] = Sci.interpolate.spline(Solid.meanline_c0[0.5*Body.N:Body.N+1], colPF[:,1], Solid.nodes[:,2])
            nodalInput[:,5] = Sci.interpolate.spline(Solid.meanline_c0[0.5*Body.N:Body.N+1], colM[:,0], Solid.nodes[:,2])

        # Rotate force components into the relative cooridnate system
        (nodalInput[:,0], nodalInput[:,1]) = self.rotatePts(nodalInput[:,0], nodalInput[:,1], theta)

        
        # Create the load matrix
        Fload = np.zeros((3*(Solid.Nnodes),1))
        Fload[0::3,0] = nodalInput[:,0]
        Fload[1::3,0] = nodalInput[:,1]
        Fload[2::3,0] = nodalInput[:,5]
        
        # Create element area matrix
        A = Solid.tBeamStruct[:,0]
        
        # Create area moment of inertia matrix
        I = 1. * Solid.tBeamStruct[:,0]**3 / 12
        
        # Set fixed nodes
        fixedNodes = Solid.fixedCounter
        
        # Initial element length
        l_0 = C / Solid.Nelements
        
        # Initial displacements and velocities
        if (i_t <= 2 and outerCorr <= 1):
            PyFEA.U_n = np.zeros((3*(Solid.Nnodes),1))
            PyFEA.Udot_n = np.zeros((3*(Solid.Nnodes),1))
            #UdotDot_n = np.zeros((3*(Solid.Nelements+1),1))
            #U_nPlus = np.zeros((3*(Solid.Nelements+1),1))
            #Udot_nPlus = np.zeros((3*(Solid.Nelements+1),1))
            #UdotDot_nPlus = np.zeros((3*(Solid.Nelements+1),1))
        if (i_t > 2 and outerCorr <= 1):
            PyFEA.U_n.resize((3*(Solid.Nnodes),1))
            PyFEA.Udot_n.resize((3*(Solid.Nnodes),1))
#            PyFEA.U_n = np.zeros()
#            PyFEA.Udot_n = np.zeros((3*(Solid.Nnodes),1))
            PyFEA.U_n[3*fixedNodes:,0] = PyFEA.U_nPlus.T
            PyFEA.Udot_n[3*fixedNodes:,0] = PyFEA.Udot_nPlus.T
        
        # Resize matricies to acount for all nodes after first subiteration
        PyFEA.Fload.resize((3*Solid.Nnodes,1))
        
        PyFEA.Fload = Fload
        PyFEA.A = A
        PyFEA.I = I
        PyFEA.l = l_0 * np.ones(Solid.Nelements)
            
        return relXp, relXc
            
    def getDisplacements(self, Solid, Body, PyFEA, t, TSTEP, SWITCH_INTERP_MTD, FLEX_RATIO):
        """
        Returns the calculated displacements of the fluid body based on the 
        displacements calculated by the structural body. This is used to 
        calculate the FSI coupling residual used in an interative strong 
        coupling between the fluid and solid solutions.
        
        Keyword arguments:
        Solid --
        Body --
        PyFEA --
        t --
        TSTEP --
        SWITCH_INTERP_MTD -- 
        FLEX_RATIO --
        """
        # Determine current pitching angle and heave position
        # TODO: Add heaving functionality to kinematics
        theta = Body.THETA_MAX * np.sin(2 * np.pi * Body.F * (t + TSTEP) + Body.PHI)  
        heave = 0      
        
        # Get the absolute x and z displacements
        nodeDisplacements = np.zeros((Solid.Nnodes-Solid.fixedCounter,2))
        nodeDisplacements[:,1] =  PyFEA.U_nPlus[1::3]
        nodeDisplacements[:,0] =  (Solid.nodes_0[Solid.fixedCounter:,1] + \
                                  nodeDisplacements[:,1]) * np.sin(-PyFEA.U_nPlus[2::3] + \
                                  PyFEA.U_nPlus[::3])
        
        # Calculate the new structural locations
        tempNodes = Solid.nodes_0;
        tempNodes[Solid.fixedCounter:,0] = Solid.nodes_0[Solid.fixedCounter:,0] + nodeDisplacements[:,0] # New x-position
        tempNodes[Solid.fixedCounter:,1] = Solid.nodes_0[Solid.fixedCounter:,1] + nodeDisplacements[:,1] # New z-position
        
        #frameNodes = tempNodes
        
        tempNodes[:,0], tempNodes[:,1] = self.rotatePts(tempNodes[:,0], tempNodes[:,1], theta)
        tempNodes[:,1] = tempNodes[:,1] + heave
        
        # Calculating the shift in node positions with the swimming velocity
        nodeDelxp = Body.xLE[1] * np.ones((Solid.Nnodes,1))
        nodeDelzp = Body.zLE[1] * np.ones((Solid.Nnodes,1))

        tempNodes[:,0] = tempNodes[:,0] + nodeDelxp.T
        tempNodes[:,1] = tempNodes[:,1] + nodeDelzp.T
              
        normal = np.zeros((Solid.Nelements,2))
        for i in xrange(Solid.Nelements):
            normal[i,0] = -1 * (tempNodes[i+1,1] - tempNodes[i,1])
            normal[i,1] =  1 * (tempNodes[i+1,0] - tempNodes[i,0])
            normal[i,:] = normal[i,:] / np.sqrt(normal[i,0]**2 + normal[i,1]**2)
        
        topNodes = np.zeros((Solid.Nnodes, 3))
        bottomNodes = np.zeros((Solid.Nnodes, 3)) 
        
        topNodes[:-1,0] = tempNodes[:-1,0] + 0.5 * Solid.tBeam[:,0] * normal[:,0]
        topNodes[:-1,1] = tempNodes[:-1,1] + 0.5 * Solid.tBeam[:,0] * normal[:,1]
        topNodes[-1,0:2] = tempNodes[-1,0:2]
        topNodes[:,2] = tempNodes[:,2]
        bottomNodes[:-1,0] = tempNodes[:-1,0] - 0.5 * Solid.tBeam[:,0] * normal[:,0]
        bottomNodes[:-1,1] = tempNodes[:-1,1] - 0.5 * Solid.tBeam[:,0] * normal[:,1]
        bottomNodes[-1,0:2] = tempNodes[-1,0:2]
        bottomNodes[:,2] = tempNodes[:,2]
        
        # Interpolate the structual results back to the fluid domain
        i = np.rint(0.5 * np.shape(Solid.meanline_p0)[0])
        if (SWITCH_INTERP_MTD == 1):
            bottomXp = np.interp(Solid.meanline_p0[0:i], bottomNodes[:,2], bottomNodes[:,0], left=True, right=True)
            bottomZp = np.interp(Solid.meanline_p0[0:i], bottomNodes[:,2], bottomNodes[:,1], left=True, right=True)
            topXp = np.interp(Solid.meanline_p0[i:], topNodes[:,2], topNodes[:,0], left=True, right=True)
            topZp = np.interp(Solid.meanline_p0[i:], topNodes[:,2], topNodes[:,1], left=True, right=True)
        else:
            bottomXp = Sci.interpolate.spline(bottomNodes[:,2], bottomNodes[:,0], Solid.meanline_p0[0:i])      
            bottomZp = Sci.interpolate.spline(bottomNodes[:,2], bottomNodes[:,1], Solid.meanline_p0[0:i])
            topXp = Sci.interpolate.spline(topNodes[:,2], topNodes[:,0], Solid.meanline_p0[i:])
            topZp = Sci.interpolate.spline(topNodes[:,2], topNodes[:,1], Solid.meanline_p0[i:])
        
        newxp = np.zeros_like(Solid.meanline_p0)
        newzp = np.zeros_like(Solid.meanline_p0)
        newxp[0:i] = bottomXp
        newxp[i:] = topXp
        newzp[0:i] = bottomZp
        newzp[i:] = topZp       
        
        for i in xrange(Body.N+1):
            if (Solid.meanline_p0[i] <= FLEX_RATIO):
                newxp[i] = Body.x[i]
                newzp[i] = Body.z[i]
        
#       Store the absolute displacements and temporary nodes.        
        self.DU[:,0] = newxp - Body.x
        self.DU[:,1] = newzp = Body.z
        self.maxDU = np.max(np.sqrt(self.DU[:,0]**2 + self.DU[:,1]**2))
        Solid.tempNodes = tempNodes
        
    def calcFSIResidual(self, Solid, outerCorr):
        self.solidNodeDispl = self.DU
        self.fsiResidualOld = self.fsiResidual
        self.nodeResidualOld = self.nodeResidual
        self.fsiResidual = self.solidNodeDispl - self.fluidNodeDispl       
        self.nodeResidual = (Solid.tempNodes[:,0:2]-Solid.nodes[:,0:2]) - self.nodeDispl        
        magFsiResidual = np.sqrt(self.fsiResidual[:,0]**2 + self.fsiResidual[:,1]**2)
        self.fsiResidualNorm = np.linalg.norm(self.fsiResidual, ord=2)
        self.maxFsiResidualNorm = np.linalg.norm(self.fsiResidual, ord=np.inf)
        if (outerCorr == 1 ):
            self.initialFsiResidualNorm = self.fsiResidualNorm
            self.maxInitialFsiResidualNorm = self.maxFsiResidualNorm

        self.fsiResidualNorm = self.fsiResidualNorm / self.initialFsiResidualNorm
        self.maxFsiResidualNorm = self.maxFsiResidualNorm / self.maxInitialFsiResidualNorm
        self.maxMagFsiResidual = np.max(magFsiResidual)