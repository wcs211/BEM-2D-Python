# -*- coding: utf-8 -*-
"""Module for the Swimmer class and its methods."""

import numpy as np
from general_functions import panel_vectors, transformation
from swimmer_subclasses import Body, Edge, Wake

class Swimmer(object):
    """A single swimmer, consisting of Body, Edge, and Wake objects.
    
    Currently the swimmer is made of a single Body, Edge, and Wake. Separation
    occurs at the trailing edge only, hence the single Edge and single Wake
    behind it. In the future this class will likely need to possess multiple
    Edges where separation occurs, so the methods would have to be
    changed accordingly.
    
    Attributes:
        CE: Constant that determines the length of Edge panels.
        DELTA_CORE: Constant used in wake_rollup to avoid singularities near
            wake panels.
        SW_GEOMETRY: Switch for Body geometry type (currently VDV only).
        SW_KUTTA: Switch for Kutta condition (explicit or unsteady).
        V0: The freestream velocity.
    """
    def __init__(self, SwimmerParameters, GeoParameters, MotionParameters, N_WAKE):
        """Inits Swimmer with all necessary parameters.
        
        This is also where Body, Edge, and Wake are created.
        GeoParameters and MotionParameters get passed along into Body creation.
        """
        self.V0 = MotionParameters.V0
        self.CE = SwimmerParameters.CE
        self.DELTA_CORE = SwimmerParameters.DELTA_CORE
        self.SW_GEOMETRY = SwimmerParameters.SW_GEOMETRY
        self.SW_KUTTA = SwimmerParameters.SW_KUTTA
        
        if self.SW_GEOMETRY == 'VDV':
            self.Body = Body.from_van_de_vooren(GeoParameters, MotionParameters)
            
        self.Edge = Edge(self.CE)
        self.Wake = Wake(N_WAKE)
    
    def edge_shed(self, DEL_T, i):
        """Updates the position of the Edge panel.
        
        The edge panel's length is based on Edge.CE.
        
        Args:
            DEL_T: Time step size.
            i: Time step number.
            Body: Body of the swimmer.
            Edge: Edge panel of separation.
        """
        Body = self.Body
        Edge = self.Edge
        
        if i == 0:
            pass
        
        else:
            Edge.x[0] = Body.AF.x[0]
            Edge.z[0] = Body.AF.z[0]
            Edge.x[1] = Body.AF.x[0] + Edge.CE*panel_vectors(Body.AF.x_neut,Body.AF.z_neut)[0][0]*Body.V0*DEL_T
            Edge.z[1] = Body.AF.z[0] + Edge.CE*panel_vectors(Body.AF.x_neut,Body.AF.z_neut)[1][0]*Body.V0*DEL_T
        
    def wake_shed(self, DEL_T, i):
        """Updates the positions of the Wake panels.
        
        Args:
            DEL_T: Time step size.
            i: Time step number.
            Edge: Edge panel of separation.
            Wake: Wake panels.
            V0: Freestream velocity.
        """
        Edge = self.Edge
        Wake = self.Wake
        V0 = self.V0
        
        # First timestep is before t==0 (geometry initialization timestep)
        if i == 0:
            pass
        
        # Initialize wake coordinates when i==1
        elif i == 1:
            
            Wake.x[0] = Edge.x[-1]
            Wake.z[0] = Edge.z[-1]
            
            Wake.x[1:] = Wake.x[0] + np.arange(1,np.size(Wake.x))*(-V0)*DEL_T
            Wake.z[1:] = Wake.z[0]
        
        else:
            Wake.x[1:] = Wake.x[:-1]
            Wake.z[1:] = Wake.z[:-1]
            Wake.mu[1:] = Wake.mu[:-1]
            
            Wake.x[0] = Edge.x[-1]
            Wake.z[0] = Edge.z[-1]
            Wake.mu[0] = Edge.mu
            
            Wake.gamma[0] = -Wake.mu[0]
            Wake.gamma[1:-1] = Wake.mu[:-1]-Wake.mu[1:]
            Wake.gamma[-1] = Wake.mu[-1]

    def wake_rollup(self, DEL_T, i):
        """Performs wake rollup on a swimmer's wake panels.
        
        Args:
            DEL_T: Time step size.
            i: Time step number.
            Body: Body of the swimmer.
            Edge: Edge panel of separation.
            Wake: Wake panels.
            DELTA_CORE: Constant used to avoid singularities when getting too close
                to other panels.
        """
        Body = self.Body
        Edge = self.Edge
        Wake = self.Wake
        DELTA_CORE = self.DELTA_CORE
        
        # Wake panels initialize when i==2
        if i < 2:
            pass
        
        else:
            
            NT = i-1 # Number of targets (wake panel points that are rolling up)
            
            vx = np.zeros(NT)
            vz = np.zeros(NT)
            
            # Coordinate transformation for body panels influencing wake
            (xp1,xp2,zp) = transformation(Wake.x[1:i],Wake.z[1:i],Body.AF.x,Body.AF.z)
            
            # Angle of normal vector with respect to global z-axis
            (nx,nz) = panel_vectors(Body.AF.x,Body.AF.z)[2:4]
            beta = np.arctan2(-nx,nz)
            
            # Katz-Plotkin eqns 10.20 and 10.21 for body source influence
            dummy1 = np.transpose(np.log((xp1**2+zp**2)/(xp2**2+zp**2))/(4*np.pi))
            dummy2 = np.transpose((np.arctan2(zp,xp2)-np.arctan2(zp,xp1))/(2*np.pi))
            
            # Rotate back to global coordinates
            dummy3 = dummy1*np.cos(beta) - dummy2*np.sin(beta)
            dummy4 = dummy1*np.sin(beta) + dummy2*np.cos(beta)
            
            # Finish eqns 10.20 and 10.21 for induced velocity by multiplying with Body.sigma
            vx = np.dot(dummy3,Body.sigma)
            vz = np.dot(dummy4,Body.sigma)
            
            # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
            NI = Body.N+1
            xp = np.repeat(Wake.x[1:i,np.newaxis].T,NI,0) - np.repeat(Body.AF.x[:,np.newaxis],NT,1)
            zp = np.repeat(Wake.z[1:i,np.newaxis].T,NI,0) - np.repeat(Body.AF.z[:,np.newaxis],NT,1)
            
            # Find distance r_b between each influence/target
            r_b = np.sqrt(xp**2+zp**2)
            
            # Katz-Plotkin eqns 10.9 and 10.10 for body doublet (represented as point vortices) influence
            dummy1 = np.transpose(zp/(2*np.pi*(r_b**2+DELTA_CORE**2)))
            dummy2 = np.transpose(-xp/(2*np.pi*(r_b**2+DELTA_CORE**2)))
            
            # Finish eqns 10.9 and 10.10 by multiplying with Body.gamma, add to induced velocity
            vx += np.dot(dummy1,Body.gamma)
            vz += np.dot(dummy2,Body.gamma)
            
            # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
            NI = Edge.N+1
            xp = np.repeat(Wake.x[1:i,np.newaxis].T,NI,0) - np.repeat(Edge.x[:,np.newaxis],NT,1)
            zp = np.repeat(Wake.z[1:i,np.newaxis].T,NI,0) - np.repeat(Edge.z[:,np.newaxis],NT,1)
            
            # Find distance r_e between each influence/target
            r_e = np.sqrt(xp**2+zp**2)
            
            # Katz-Plotkin eqns 10.9 and 10.10 for edge (as point vortices) influence
            dummy1 = np.transpose(zp/(2*np.pi*(r_e**2+DELTA_CORE**2)))
            dummy2 = np.transpose(-xp/(2*np.pi*(r_e**2+DELTA_CORE**2)))
            
            # Finish eqns 10.9 and 10.10 by multiplying with Edge.gamma, add to induced velocity
            vx += np.dot(dummy1,Edge.gamma)
            vz += np.dot(dummy2,Edge.gamma)
            
            # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
            NI = i
            xp = np.repeat(Wake.x[1:i,np.newaxis].T,NI,0) - np.repeat(Wake.x[:i,np.newaxis],NT,1)
            zp = np.repeat(Wake.z[1:i,np.newaxis].T,NI,0) - np.repeat(Wake.z[:i,np.newaxis],NT,1)
            
            # Find distance r_w between each influence/target
            r_w = np.sqrt(xp**2+zp**2)
            
            # Katz-Plotkin eqns 10.9 and 10.10 for wake (as point vortices) influence
            dummy1 = np.transpose(zp/(2*np.pi*(r_w**2+DELTA_CORE**2)))
            dummy2 = np.transpose(-xp/(2*np.pi*(r_w**2+DELTA_CORE**2)))
            
            # Finish eqns 10.9 and 10.10 by multiplying with Wake.gamma, add to induced velocity
            vx += np.dot(dummy1,Wake.gamma[:i])
            vz += np.dot(dummy2,Wake.gamma[:i])
            
            # Modify wake with the total induced velocity
            Wake.x[1:i] += vx*DEL_T
            Wake.z[1:i] += vz*DEL_T
            
    def influence_matrices(self, i):
        """Constructs the influence coefficient matrices.
        
        Args:
            i: Time step number.
            Body: Body of the swimmer.
            Edge: Edge panel of separation.
            Wake: Wake panels.
        """
        Body = self.Body
        Edge = self.Edge
        Wake = self.Wake
        
        if i > 0:
            # Tangential and normal body panels and panel length calculations
            (xp1,xp2,zp) = transformation(Body.AF.x_col,Body.AF.z_col,Body.AF.x,Body.AF.z)
        
            # Body source singularities influencing the body
            # Transpose so that row elements represent the effect on the (row number)th panel
            Body.phi_s = np.transpose((xp1 * np.log(xp1**2 + zp**2) - xp2 * np.log(xp2**2 + zp**2) \
                                      + 2*zp*(np.arctan2(zp,xp2) - np.arctan2(zp,xp1)))/(4*np.pi))
            
            # Body source strength calculations
            (nx,nz) = panel_vectors(Body.AF.x,Body.AF.z)[2:4]
            Body.sigma = nx*(Body.V0 + Body.vx) + nz*Body.vz   # normal vector pointing outward (overall sigma pointing outward)
            
            # Body doublet singularities influencing body itself
            # Transpose similar to phi_s
            Body.phi_db = np.transpose(-(np.arctan2(zp,xp2)\
                                      - np.arctan2(zp,xp1))/(2*np.pi))
            
            # Edge doublet influencing the body
            (xp1,xp2,zp) = transformation(Body.AF.x_col,Body.AF.z_col,Edge.x,Edge.z)
            
            if i > 1: # No wake panels until i==2
                # Wake doublets influencing the body
                (xp1_w,xp2_w,zp_w) = transformation(Body.AF.x_col,Body.AF.z_col,Wake.x[:i],Wake.z[:i])
                # Join edge and wake doublet influences into a single matrix
                xp1 = np.insert(xp1, Edge.N, xp1_w, axis=0)
                xp2 = np.insert(xp2, Edge.N, xp2_w, axis=0)
                zp = np.insert(zp, Edge.N, zp_w, axis=0)
        
            Body.phi_dw = np.transpose(-(np.arctan2(zp,xp2) \
                                       -np.arctan2(zp,xp1))/(2*np.pi))
            
    def kutta(self, RHO, SWITCH_KUTTA, DEL_T, i):
        """Applies Kutta condition to the swimmer's trailing edge.
        
        Args:
            RHO: Fluid density.
            SWITCH_KUTTA: Switch for either explicit (0) or unsteady (1) Kutta
                condition.
            DEL_T: Time step size.
            i: Time step number.
            Body: Body of the swimmer.
            Edge: Edge panel of separation.
            Wake: Wake panels.
        """
        Body = self.Body
        Edge = self.Edge
        Wake = self.Wake
        
        if i == 0:
            pass
        else:
            n_iter = 0
            mu_guess = np.empty(2) # [0] is current guess, [1] is previous
            delta_cp = np.empty(2) # [0] is current delta_cp, [1] is previous
            
            while True:
                n_iter += 1
                
                if n_iter == 1:
                    # Begin with explicit Kutta condition as first guess
                    # Construct the augmented body matrix by combining body and trailing edge panel arrays
                    c = np.zeros((Body.N,Body.N))
                    c[:,0] = -Body.phi_dw[:,0] # c contains the explicit Kutta condition (first and last columns)
                    c[:,-1] = Body.phi_dw[:,0]
                    Body.phi_dinv = np.linalg.inv(Body.phi_db + c)
                    # Get rhs
                    if i == 1:
                        rhs = -np.dot(Body.phi_s,Body.sigma)
                    else:
                        rhs = -np.dot(Body.phi_s,Body.sigma) - np.dot(Body.phi_dw[:,1:],Wake.mu[:i-1])
                    # Solve for body doublet strengths using explicit Kutta
                    Body.mu = np.dot(Body.phi_dinv,rhs)
                    # First mu_guess (from explicit Kutta)
                    mu_guess[0] = Body.mu[-1]-Body.mu[0]
                    
                else:
                    if n_iter == 2: # Make a second initial guess
                        # Update phi_dinv so it no longer includes explicit Kutta condition
                        Body.phi_dinv = np.linalg.inv(Body.phi_db)
                        
                        mu_guess[1] = mu_guess[0]
                        delta_cp[1] = delta_cp[0]
                        mu_guess[0] = 0.8*mu_guess[1] # Multiply first (explicit) guess by arbitrary constant to get second guess
                        
                    else: # Newton method to get delta_cp == 0
                        # Get slope, which is change in delta_cp divided by change in mu_guess
                        slope = (delta_cp[0]-delta_cp[1])/(mu_guess[0]-mu_guess[1])
                        mu_guess[1] = mu_guess[0]
                        delta_cp[1] = delta_cp[0]
                        mu_guess[0] = mu_guess[1] - delta_cp[0]/slope
                    
                    # Form right-hand side including mu_guess as an influence
                    if i == 1:
                        rhs = -np.dot(Body.phi_s,Body.sigma) - np.squeeze(np.dot(Body.phi_dw,mu_guess[0]))
                    else:              
                        rhs = -np.dot(Body.phi_s,Body.sigma) - np.dot(Body.phi_dw,np.insert(Wake.mu[:i-1],0,mu_guess[0]))
                    
                    Body.mu = np.dot(Body.phi_dinv,rhs)
                    
                
                Body.pressure(RHO, DEL_T, i)
                if SWITCH_KUTTA == 0:
                    break
                delta_cp[0] = np.absolute(Body.cp[-1]-Body.cp[0])
                if delta_cp[0] < 0.0001:
                    break
                        
            # mu_past used in differencing for pressure
            Body.mu_past[1,:] = Body.mu_past[0,:]
            Body.mu_past[0,:] = Body.mu
            
            Edge.mu = mu_guess[0]
            Edge.gamma[0] = -Edge.mu
            Edge.gamma[1] = Edge.mu
            
            # Get gamma of body panels for use in wake rollup
            Body.gamma[0] = -Body.mu[0]
            Body.gamma[1:-1] = Body.mu[:-1]-Body.mu[1:]
            Body.gamma[-1] = Body.mu[-1]