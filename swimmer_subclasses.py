#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code. Module for the Body, Edge, and Wake classes.

"""
import numpy as np
from functions_general import panel_vectors
import parameter_classes as PC
    
def finite_diff(mu, dL, stencil):
    """
    Method for generating arbitrary finite difference approximation based on the stencil
    Solving for the coefficients of the finite difference scheme (three-point stencil example):
        1.) dphi/dx = a*phi_m1 + b*phi_0 + c*phi_p1
        2.) Write out a Taylor series expansion of phi_m1 and phi_p1
        3.) a*phi_m1 + b*phi_0 + c*phi_p1 = (a + b + c)*phi_0 + (-a*s_m1 + c*s_p1)*phi_x + ...
                  1/2*(-a*s_m1^2 + c*s_p1^2)*phi_xx + 1/6*(-a*s_m1^3 + c*s_p1^3)*phi_xxx + ...
        4.) Set coefficient of phi and derivative of phi to 1 and 0.  For instance:               
                           (a + b + c) = 0,          (-a*s_m1 + c*s_p1) = 1, 
            1/2*(-a*s_m1^2 + c*s_p1^2) = 0,  1/6*(-a*s_m1^3 + c*s_p1^3) = 0
        5.) Write out equations for the coefficients (a,b,c) in matrix form A*coeffs = b and solve.
        6.) Construct finite difference scheme for the gradient
    
    Calculating the arc length vector from the zero stencil point to the
    other stencil points.  Breaking up the arc length vector into the minus
    and plus arc length vectors relative to the zero stencil point
    """
    # Calculating arc lengths between stencil points
    s = 0.5 * dL[:-1] + 0.5 * dL[1:]
    
    # Determine the number of minus and plus points
    n_m = np.absolute(np.min(stencil))
    n_p = np.absolute(np.max(stencil))
    
    s_m = np.zeros(n_m)
    for i in xrange(n_m):
        s_m[i] = np.sum(s[i:n_m])
        
    s_p = np.zeros(n_p)
    for i in xrange(n_p):
        s_p[i] = np.sum(s[n_m: n_m + i + 1])
        
    # Building submatrices of A for up to a 5-point stencil
    numel = stencil.shape[0]
    A = np.zeros((numel, n_m + n_p + 1))
    if (n_m > 0):
        A_m = np.zeros((numel, n_m))
        A_m[0,:] = np.ones((1,n_m))
        A_m[1,:] = -s_m
        A_m[2,:] = 0.5*s_m**2
        A_m[3,:] = -1./6.*s_m**3
        A_m[4,:] = 1./24.*s_m**4
        A[:numel,:n_m] = A_m[:numel, :]

    A[:numel,n_m ] = np.array([1., 0., 0., 0., 0.])

    if (n_p > 0):
        A_p = np.zeros((numel, n_p))
        A_p[0,:] = np.ones((1,n_p))
        A_p[1,:] = s_p
        A_p[2,:] = 0.5*s_p**2
        A_p[3,:] = 1./6.*s_p**3
        A_p[4,:] = 1./24.*s_p**4
        A[:numel, n_m + 1: n_p + n_m + 1] = A_p[:numel, :]
        
    # Determining the full vector b_full and vector b for specific stencil number of points
    b_full = np.array([[0.], [1.], [0.], [0.], [0.]])
    b = b_full[:numel]
    
    # The A matrix and b vector need to be scaled so that the matrix isn't
    # ill-conditioned.  The minimum arc length will be determined and each row 
    # of A and b will be divided by the power of s_min that will scale that
    # equation (row) to order 1.  This is get around the ill-conditioned matrix
    # warning that we are getting.
    if (s_m.shape[0] > 0):
        s1 = np.min(s_m)
    else:
        s1 = np.nan
    
    if (s_p.shape[0] > 0):
        s2 = np.min(s_p)
    else:
        s2 = np.nan
    
    s_min = np.nanmin(np.array([s1, s2]))
    vec_scaling_full = np.array([[1.], [s_min], [s_min**2], [s_min**3], [s_min**4]])
    b_vec_scaling = vec_scaling_full[:numel]
    A_mat_scaling = np.repeat(b_vec_scaling, 1 + n_m + n_p, axis=1)
    A = A / A_mat_scaling
    b = b / b_vec_scaling
    
    # Solving for finite difference approximation coefficients
    coeffs = np.linalg.solve(A, b)
    
    # Finite difference approximation of the gradient
    dmu_ds = np.dot(coeffs.T, mu)
    
    # Determining the pertubation velocity (Qp)
    return(-dmu_ds)

class Edge(object):
    """An edge doublet panel located where separation occurs on a body.

    Attributes:
        N: Number of edge panels (just one).
        CE: Constant that determines the length of the edge panel.
        x, z: X- and Z-coordinates of the edge panel endpoints.
        mu: Doublet strength of the edge panel.
        gamma: Circulation at the edge panel endpoints.
    """
    def __init__(self, CE):
        """Inits Edge with all necessary parameters."""
        self.N = 1
        self.CE = CE
        self.x = np.zeros(self.N+1)
        self.z = np.zeros(self.N+1)
        self.mu = np.zeros(self.N)
        self.gamma = np.zeros(self.N+1)

class Wake(object):
    """A chain of wake doublet panels.

    Attributes:
        N: Number of wake panels.
        x, z: X- and Z-coordinates of the wake panel endpoints.
        mu: Doublet strengths of the wake panels.
        gamma: Circulations at the wake panel endpoints.
    """
    def __init__(self, N):
        """Inits Wake with all necessary parameters."""
        self.N = N
        self.x = np.zeros(N+1)
        self.z = np.zeros(N+1)
        self.mu = np.zeros(N)
        self.gamma = np.zeros(N+1)

class Body(object):
    """An arrangement of source/doublet panels in the shape of a swimming body.

    Attributes:
        N: Number of body panels.
        S: Parameter for shifting collocation points into the body.
        BF: A collection of various body-frame coordinates.
        AF: A collection of various absolute-frame coordinates.
        MP: A collection of parameters describing the motion of the swimmer.
        V0: The free-stream velocity (included in MP as well).
        vx, vz: X- and Z-components of body-frame surface velocities.
        sigma: Source strengths of the body panels.
        phi_s: Matrix of body source panel influences on the body.
        phi_db: Matrix of body doublet panel influences on the body.
        phi_dw: Matrix of edge and wake panel influences on the body.
        mu: Doublet strengths of the body panels.
        gamma: Circulations at the body panel endpoints.
        p: Surface pressures of the body panels.
        cp: Surface pressure coefficients of the body panels.
        mu_past: mu arrays from previous time steps for backwards differencing.
    """
    def __init__(self, N, S, BodyFrameCoordinates, MotionParameters):
        """Inits Body with all necessary parameters."""
        self.N = N
        self.S = S

        # Body-frame panel coordinates:
        # x, z, x_col, z_col
        self.BF = BodyFrameCoordinates
        # Initialize absolute-frame panel coordinates:
        # x, z, x_col, z_col, x_mid, z_mid, x_neut, z_neut
        self.AF = PC.BodyAFC(N)
        # Prescribed motion
        self.MP = MotionParameters
        self.V0 = MotionParameters.V0
        self.V  = np.copy(self.V0)

        self.vx = np.zeros(N)
        self.vz = np.zeros(N)

        self.sigma = np.zeros(N)
        self.phi_s = np.zeros((N,N))
        self.phi_db = np.zeros((N,N))
        self.phi_dw = np.zeros(N)
        self.mu = np.zeros(N)
        self.gamma = np.zeros(N+1)
        self.zcval = np.copy(self.BF.z_col)

        self.p = np.zeros(N)
        self.cp = np.zeros(N)
        self.mu_past = np.zeros((4,N))
        
        self.Cf = 0.
        self.Cl = 0.
        self.Ct = 0.
        self.Ct_net = 0.
        self.Cd_visc = 0.
        self.D_visc = 0.
        self.Cpow = 0.
        self.forceData = np.zeros((0,7))

    @classmethod
    def from_van_de_vooren(cls, GeoVDVParameters, MotionParameters):
        """Creates a Body object based on a Van de Vooren airfoil geometry.

        MotionParameters are unused here, just getting passed through for the
        creation of the Body.

        Args:
            GeoVDVParameters: A collection of parameters for constructing a
                Van de Vooren geometry. (N, S, C, K, EPSILON)
            MotionParameters: Motion parameters of the swimmer.

        Returns:
            A Body object with the Van de Vooren airfoil geometry.
        """
        N = GeoVDVParameters.N
        S = GeoVDVParameters.S
        C = GeoVDVParameters.C
        K = GeoVDVParameters.K
        EPSILON = GeoVDVParameters.EPSILON

        A = C*((1+EPSILON)**(K-1))*(2**(-K))

        THETA = np.linspace(0,np.pi,N/2+1)

        R1 = np.sqrt((A*np.cos(THETA)-A)**2+(A**2)*np.sin(THETA)**2)
        R2 = np.sqrt((A*np.cos(THETA)-EPSILON*A)**2+(A**2)*np.sin(THETA)**2)

        THETA1 = np.arctan2((A*np.sin(THETA)) , (A*np.cos(THETA)-A))
        THETA2 = np.arctan2(A*np.sin(THETA) ,(A*np.cos(THETA)-EPSILON*A))

        x = ((R1**K)/(R2**(K-1)))*(np.cos(K*THETA1)*np.cos((K-1)*THETA2) + np.sin(K*THETA1)*np.sin((K-1)*THETA2))
        z_top = ((R1**K)/(R2**(K-1)))*(np.sin(K*THETA1)*np.cos((K-1)*THETA2) - np.cos(K*THETA1)*np.sin((K-1)*THETA2))
        z_bot = -((R1**K)/(R2**(K-1)))*(np.sin(K*THETA1)*np.cos((K-1)*THETA2) - np.cos(K*THETA1)*np.sin((K-1)*THETA2))

        x = x-x[-1] # Carrying the leading edge to the origin
        x[0] = C

        z_top[0] = 0
        z_bot[0] = 0
        z_bot[-1] = 0

        # Merge top and bottom surfaces together
        x = np.hstack((x , x[-2::-1]))
        z = np.hstack((z_bot , z_top[-2::-1]))

        x_col = ((x[1:] + x[:-1])/2)
        z_col = ((z[1:] + z[:-1])/2)

        BodyFrameCoordinates = PC.BodyBFC(x, z, x_col, z_col)

        return Body(N, S, BodyFrameCoordinates, MotionParameters)

    @classmethod
    #Flat plate geometry
    def flat_plate(cls, GeoFPParameters, MotionParameters):
        N = GeoFPParameters.N
        S = GeoFPParameters.S
        C = GeoFPParameters.C
        D = GeoFPParameters.D

        # Stepping through each spanwise position to calculate the positions of the
        # fin neutral plane at the given time step.
        start = 0
        stop  = np.pi
        step  = (N+2)/2
        theta = np.linspace(start,stop,step)
        xb = (C*np.cos(theta).T + C)/(2.)

        start = np.pi
        stop  = 0
        step  = (N+2)/2
        theta = np.linspace(start,stop,step)
        xt = (C*np.cos(theta).T + C)/(2.)
        zb = -0.5*D*np.ones((N+2)/2)
        zt =  0.5*D*np.ones((N+2)/2)

        # Circle origins ( z is assumed to be zero)
        oF = 0.5 * D
        oB = C - 0.5 * D

        # Calculate new theta positions for points on the rounded ends
        count1 = np.shape(xb[xb>=oB])[0]
        count2 = np.shape(xb[xb<=oF])[0]
        count3 = np.shape(xt[xt<=oF])[0]
        count4 = np.shape(xt[xt>=oB])[0]

        thetab = np.linspace(0,np.pi,count1+count2).T
        thetat = np.linspace(np.pi,2*np.pi,count3+count4).T

        # Calculate transform leading and trailing edge points
        x1 = oB + 0.5 * D * np.cos(thetab[0:count1])
        z1 = 0  - 0.5 * D * np.sin(thetab[0:count1])
        x2 = oF + 0.5 * D * np.cos(thetab[-1:-1-count1+1:-1])
        z2 = 0  - 0.5 * D * np.sin(thetab[-1:-1-count1+1:-1])
        x3 = oF + 0.5 * D * np.cos(thetat[0:count3])
        z3 = 0  - 0.5 * D * np.sin(thetat[0:count3])
        x4 = oB + 0.5 * D * np.cos(thetat[-1:-1-count3+1:-1])
        z4 = 0  - 0.5 * D * np.sin(thetat[-1:-1-count3+1:-1])

        # Replace x and z transformed points
        xb[:count1] = x1
        xb[-1:-1-count2+1:-1] = x2
        xt[:count3] = x3
        xt[-1:-1-count4+1:-1] = x4

        zb[:count1] = z1
        zb[-1:-1-count2+1:-1] = z2
        zt[:count3] = z3
        zt[-1:-1-count4+1:-1] = z4

        zb[0] = 0
        zt[0] = 0
        zb[-1] = 0

        # Merge top and bottom surfaces together
        xb = np.hstack((xb , xb[-2::-1]))
        zb = np.hstack((zb , zt[-2::-1]))

        xb_col = ((xb[1:] + xb[:-1])/2)
        zb_col = ((zb[1:] + zb[:-1])/2)

        BodyFrameCoordinates = PC.BodyBFC(xb, zb, xb_col, zb_col)

        return Body(N, S, BodyFrameCoordinates, MotionParameters)

    @classmethod
    #Tear-drop geometry
    def tear_drop(cls, GeoTDParameters, MotionParameters):
        N = GeoTDParameters.N
        S = GeoTDParameters.S
        C = GeoTDParameters.C
        D = GeoTDParameters.D
        
        # Stepping through each spanwise position to calculate the positions of
        # the fin neutral plane at the given time step.
        xb = np.linspace(np.pi, 0., (N+2.)/2)
        xt = np.linspace(0., np.pi, (N+2.)/2)

        # Slopes and intersects for the line segments
        m = -D/2/(C - D/2)
        b = D/2 + D**2/4/(C - D/2)

        # Tear drop shape equation.
        x_c = 0.5 * (1 - np.cos(xb))
        xb = x_c * C
        xb1 = xb[xb <= D/2]
        xb2 = xb[xb > D/2]

        zb2 = -m * xb2 - b
        zb1 = -np.sqrt((D/2)**2 - (xb1 - D/2)**2)
        zb = np.hstack((zb2, zb1))

        # Tear drop shape equation.
        x_c = 0.5 * (1 - np.cos(xt))
        xt = x_c * C
        xt1 = xt[xt <= D/2]
        xt2 = xt[xt > D/2]

        zt1 = np.sqrt((D/2)**2 - (xt1 - D/2)**2)
        zt2 = m * xt2 + b
        zt = np.hstack((zt1, zt2))

        zb[0] = 0
        zt[0] = 0
        zb[-1] = 0

        # Merge top and bottom surfaces together
        x = np.hstack((xb , xt[1:]))
        z = np.hstack((zb , zt[1:]))

        x_col = ((x[1:] + x[:-1])/2)
        z_col = ((z[1:] + z[:-1])/2)

        BodyFrameCoordinates = PC.BodyBFC(x, z, x_col, z_col)

        return Body(N, S, BodyFrameCoordinates, MotionParameters)

    def neutral_axis(self, x, T, THETA, HEAVE, DSTEP=0):
        """Finds a body's neutral axis for a given time.

        The neutral axis is the axis which coincides with the chord line and
        divides the symmetric airfoil into two.

        The axis that it finds is in an absolute frame of reference.

        Args:
            x: An array of body-frame x-coordinates to use as reference points.
            DSTEP: Small incremental distance offset (intended for differencing).
            T: Time of the current step.
            THETA: Current pitching angle.

        Returns:
            x_neut and z_neut: X- and Z-coordinates of the neutral axis points.
        """
        X0 = self.MP.X0
        Z0 = self.MP.Z0
        V0 = self.MP.V0

        x_neut = X0 + (x+DSTEP)*np.cos(THETA) + V0*T
        z_neut = Z0 + (x+DSTEP)*np.sin(THETA) + HEAVE

        return(x_neut, z_neut)

    def panel_positions(self, P, i):
        """Updates all the absolute-frame coordinates of the body.

        Args:
            DSTEP: Small incremental distance to pass into neutral_axis().
            T: Time of current step.
            THETA: Current pitching angle.

        """
        T = P['T'][i]
        THETA = P['THETA'][i]
        HEAVE = P['HEAVE'][i]
        bfx = self.BF.x
        bfz = self.BF.z
        
#        afx = bfx * np.cos(THETA) - bfz * np.sin(THETA) + self.V0*T
        afx = bfx * np.cos(THETA) - bfz * np.sin(THETA) + self.AF.x_le
        
        afz = bfx * np.sin(THETA) + bfz * np.cos(THETA) + HEAVE

        (x_neut, z_neut) = self.neutral_axis(bfx, T, THETA, HEAVE)

        # Absolute-frame panel midpoint positions
        x_mid = (afx[:-1]+afx[1:])/2
        z_mid = (afz[:-1]+afz[1:])/2

        # Collocation points are the points where impermeable boundary condition is forced
        # They should be shifted inside or outside of the boundary depending on the dirichlet or neumann condition
        # Shifting surface collocation points some percent of the height from the neutral axis
        # Normal vectors point outward but positive S is inward, so the shift must be subtracted from the panel midpoints
        afx_col = x_mid - self.S*panel_vectors(afx, afz)[2]*np.absolute(self.zcval) #bfz_col
        afz_col = z_mid - self.S*panel_vectors(afx, afz)[3]*np.absolute(self.zcval)

        self.AF.x = afx
        self.AF.z = afz
        self.AF.x_col = afx_col
        self.AF.z_col = afz_col
        self.AF.x_mid[0,:] = x_mid
        self.AF.z_mid[0,:] = z_mid
        self.AF.x_neut = x_neut
        self.AF.z_neut = z_neut
        
    def fsi_panel_positions(self, FSI, P, i):
        T     = P['T'][i]
        THETA = P['THETA'][i]
        HEAVE = P['HEAVE'][i]
        
        self.AF.x = self.AF.x + (FSI.fluidNodeDispl[:,0] - FSI.fluidNodeDisplOld[:,0])
        self.AF.z = self.AF.z + (FSI.fluidNodeDispl[:,1] - FSI.fluidNodeDisplOld[:,1])                 

        self.AF.x_mid[0,:] = (self.AF.x[:-1] + self.AF.x[1:])/2
        self.AF.z_mid[0,:] = (self.AF.z[:-1] + self.AF.z[1:])/2

        self.BF.x = (self.AF.x - self.AF.x_le) * np.cos(-1*THETA) - (self.AF.z - self.AF.z_le) * np.sin(-1*THETA)
        self.BF.z = (self.AF.z - self.AF.z_le) * np.cos(-1*THETA) + (self.AF.x - self.AF.x_le) * np.sin(-1*THETA)
        self.BF.x_col = ((self.BF.x[1:] + self.BF.x[:-1])/2)
        self.BF.z_col = ((self.BF.z[1:] + self.BF.z[:-1])/2)

        (self.AF.x_neut, self.AF.z_neut) = self.neutral_axis(self.BF.x, T, THETA, HEAVE)

        self.AF.x_col = self.AF.x_mid[0,:] - self.S*panel_vectors(self.AF.x, self.AF.z)[2]*np.absolute(self.zcval) #self.BF.z_col
        self.AF.z_col = self.AF.z_mid[0,:] - self.S*panel_vectors(self.AF.x, self.AF.z)[3]*np.absolute(self.zcval)

    def surface_kinematics(self, P, i):
        """Calculates the body-frame surface velocities of body panels.

        Also finds the body panel source strengths based on these surface
        velocities.

        Args:
            DSTEP, TSTEP: Incremental distance/time passed into neutral_axis().
            DEL_T: Time step length.
            T: Time of current step.
            i: Time step number.
            THETA_MINUS: Pitching angle minus a small time difference (TSTEP)
            THETA_PLUS: Pitching angle plus a small time difference (TSTEP)
        """
        TSTEP       = P['TSTEP']
        THETA_MINUS = P['THETA_MINUS'][i]
        THETA_PLUS  = P['THETA_PLUS'][i]
        HEAVE_MINUS = P['HEAVE_MINUS'][i]
        HEAVE_PLUS  = P['HEAVE_PLUS'][i]
        DEL_T       = P['DEL_T']
        T           = P['T'][i]     
        
        if i == 0:
            # Use the prescribed kinematics to do a central difference over a 
            # small period of time
            x_col = self.BF.x_col
            z_col = self.BF.z_col
            
            xp = x_col * np.cos(THETA_PLUS) - z_col * np.sin(THETA_PLUS) + self.V*(T+TSTEP)
            zp = x_col * np.sin(THETA_PLUS) + z_col * np.cos(THETA_PLUS) + HEAVE_PLUS
            
            xm = x_col * np.cos(THETA_MINUS) - z_col * np.sin(THETA_MINUS) + self.V*(T-TSTEP)
            zm = x_col * np.sin(THETA_MINUS) + z_col * np.cos(THETA_MINUS) + HEAVE_MINUS
            
            self.vx = (xp - xm) / (2. * TSTEP) - self.V
            self.vz = (zp - zm) / (2. * TSTEP)

        elif i == 1:
            # First-order backwards differencing of body collocation point positions
            self.vx = (self.AF.x_mid[0,:]-self.AF.x_mid[1,:])/DEL_T - self.V
            self.vz = (self.AF.z_mid[0,:]-self.AF.z_mid[1,:])/DEL_T
            
        elif i == 2 or i == 3:
            # Second-order backwards differencing of body collocation point positions
            self.vx = (3*self.AF.x_mid[0,:]-4*self.AF.x_mid[1,:]+self.AF.x_mid[2,:])/(2*DEL_T) - self.V
            self.vz = (3*self.AF.z_mid[0,:]-4*self.AF.z_mid[1,:]+self.AF.z_mid[2,:])/(2*DEL_T)

        else:
            # Fourth-order backwards differencing of body collocation point positions
            self.vx = (25/12*self.AF.x_mid[0,:] - 4*self.AF.x_mid[1,:] + 3*self.AF.x_mid[2,:] - 4/3*self.AF.x_mid[3,:] + 1/4*self.AF.x_mid[4,:]) / DEL_T - self.V
            self.vz = (25/12*self.AF.z_mid[0,:] - 4*self.AF.z_mid[1,:] + 3*self.AF.z_mid[2,:] - 4/3*self.AF.z_mid[3,:] + 1/4*self.AF.z_mid[4,:]) / DEL_T

        # Body source strengths with normal vector pointing outward (overall sigma pointing outward)
        (nx,nz) = panel_vectors(self.AF.x,self.AF.z)[2:4]
        self.sigma = nx*(self.V + self.vx) + nz*self.vz    

    def pressure(self, RHO, DEL_T, i, stencil_npts=5):
        """Calculates the pressure distribution along the body's surface.

        Args:
            RHO: Fluid density.
            DEL_T: Time step length.
            i: Time step number.
        """

        (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x,self.AF.z)
        
        # Defining stencil
        if (stencil_npts == 5):
            stencil_chordwise = np.zeros((self.N, 5),dtype=int)
            stencil_chordwise[0:2,:] = np.array([[0, 1, 2, 3, 4], [-1, 0, 1, 2, 3]])
            stencil_chordwise[2:-2,:] = np.tile(np.array([-2, -1, 0, 1, 2]), (self.N-4, 1))
            stencil_chordwise[-2:,:] = np.array([[-3, -2, -1, 0, 1], [-4, -3, -2, -1, 0]])
        else:
            if (stencil_npts != 3):
                print '+-----------------------------------------------------------------------------+'
                print '| WARNING! There are only three- and five-point stencils available.           |'
                print '|          Defaulting to the three-point stencil.                             |'
                print '+-----------------------------------------------------------------------------+'
            stencil_chordwise = np.zeros((self.N, 3))
            stencil_chordwise[0,:] = np.array([0, 1, 2])
            stencil_chordwise[1:-1,:] = np.repeat(np.array([-1, 0, 1]), self.N-2, axis=0)
            stencil_chordwise[-1,:] = np.array([-2, -1, 0])
        
        # Tangential panel velocity dmu/dl
        dmu_dl = np.empty(self.N)
        for i in xrange(self.N):
            # Defining stencil depending on chordwise element
            stencil = np.copy(stencil_chordwise[i,:])
            pan_elem = i + stencil
            
            # Calling finite difference approximations based on stencil (3-point and 5-point available)
            dmu_dl[i] = finite_diff(self.mu[pan_elem], lpanel[pan_elem], stencil)

        # Potential change dmu/dt, second-order differencing after first time step
        if i == 0:
            dmu_dt = self.mu / DEL_T
        if i == 1:
            dmu_dt = (self.mu - self.mu_past[0,:])/DEL_T
        else:
            dmu_dt = (3*self.mu - 4*self.mu_past[0,:] + self.mu_past[1,:])/(2*DEL_T)

        # Unsteady pressure calculation (from Matlab code)
        qpx_tot = dmu_dl*tx + self.sigma*nx
        qpz_tot = dmu_dl*tz + self.sigma*nz

        self.p = -RHO*(qpx_tot**2 + qpz_tot**2)/2. + RHO*dmu_dt + RHO*(qpx_tot*(self.V+self.vx) + qpz_tot*self.vz)
        self.cp = self.p / (0.5*RHO*self.V**2)
        
    def force(self, P, i):
        """Calculates drag and lift forces acting on the body.

        Args:
            RHO (float): Fluid density
            C (float): Body's chord length
            B (float): Body's Span length
            i: Time step number.
        """
        C   = P['C']
        B   = P['B']
        RHO = P['RHO']
        CD_BOD = P['CD_BOD']
        S_W = P['S_W']
        NU = P['NU']
        L_T = P['L_T']
        SW_ADDED_DRAG = P['SW_ADDED_DRAG']
        DRAG_LAW = P['DRAG_LAW']
        (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x, self.AF.z)

        delFx = -self.p * lpanel * B * nx
        delFz = -self.p * lpanel * B * nz
        delF = np.array([delFx, delFz])
        delP = np.sum(-delF * np.array([self.vx.T, self.vz.T]), 1)
        
        force = np.sum(delF,1)
        lift = force[1]
        thrust = -force[0]
        power = np.sum(delP, 0)
        
        if SW_ADDED_DRAG:
            if DRAG_LAW == 'FORM':
                D_add = 0.5 * CD_BOD * RHO * S_W * self.V**2
            elif DRAG_LAW == 'BLASIUS':
                D_add = CD_BOD * RHO * S_W * np.absolute(self.V**1.5 * (NU / L_T)**0.5)
            else:
                print 'ERROR: Invalid drag law "%s"' % DRAG_LAW
                print 'Valid trag laws are:'
                print '    "FORM"'
                print '    "BLASIUS"'
                raise ValueError('Invalid drag law "%s"' % DRAG_LAW)
            
            net_thrust = force[0] - np.sign(self.V) * D_add
            self.Ct_net = -net_thrust / (0.5 * RHO * np.absolute(self.V)**2 * C * B)
        
        self.Cf = np.sqrt(force[0]**2 + force[1]**2) / (0.5 * RHO * np.absolute(self.V)**2 * C * B)
        self.Cl = lift /(0.5 * RHO * np.absolute(self.V)**2 * C * B)
        self.Ct = thrust / (0.5 * RHO * np.absolute(self.V)**2 * C * B)
        self.Cpow = power /  (0.5 * RHO * np.absolute(self.V)**3 * C * B)
        
    def free_swimming(self, P, i):
        """Determines the free-swimming velocity.
        
        Args:
            T (float): Current simulation time
            HAVE (float): Current heave position
            DEL_T (float): Time-step increment
            RHO (float): Fluid density
            C (float): Body's chord length
            B (float): Body's Span length
            M (float): Body's mass
            SW_FREE_SWIM (bool): free swimming (True) or fixed swimming (False)
        """
        T            = P['T'][i]
        HEAVE        = P['HEAVE'][i]
        DEL_T        = P['DEL_T']
        RHO          = P['RHO']
        C            = P['C']
        B            = P['B']
        M            = P['M']
        SW_FREE_SWIM = P['SW_FREE_SWIM']
        
        if SW_FREE_SWIM:
            Fx = -self.Ct_net * 0.5 * RHO * np.absolute(self.V)**2 * C * B
            a_b = Fx / M
            V0_old = np.copy(self.V)
            self.V = a_b * DEL_T + self.V
            # Location of leading edge
            self.AF.x_le = 0.5 * (self.V + V0_old) * DEL_T + self.AF.x_le
            self.AF.z_le = HEAVE
        else:
            # Location of leading edge
            self.AF.x_le = self.V*T
            self.AF.z_le = HEAVE
            
#        np.save('./x_le/%05i.npy' % i, self.V)