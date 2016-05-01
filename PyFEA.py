#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
import scipy.linalg as spla

class PyFEA(object):
    def __init__(self, Solid, SW_SPRING, FRAC_DELT, endTime, E, RHO_S):
        """
        Iniitalizes object related variables needed for other class methods.
        
        Args:
            Solid (object): A solid object created from the solid class
            FRAC_DELT (float): Fraction of the fluid solver time-step to break 
                the structural solver time-step up into.
            endTime (float): Total elapsed time for the structural solver.
            E (float): Young's Modulus of thesolid object.
            RHO_S (float): Solid object's density.
        """
        if SW_SPRING:
            self.deltaT            = FRAC_DELT * endTime
            self.I                 = 0.
            self.kappa_1           = 0.
            self.kappa_2           = 0.
            self.zeta              = 0.
            self.Nf                = 0.
            self.Ni                = 0.
            self.RHO_S             = RHO_S
            self.theta_n           = 0.
            self.thetaDot_n        = 0.
            self.thetaDotDot_n     = 0.
            self.theta_nPlus       = 0.
            self.thetaDot_nPlus    = 0.
            self.thetaDotDot_nPlus = 0.
        else:
            print Solid.Nelements
            self.Nelements = Solid.Nelements
            self.M = np.zeros((3 * (Solid.Nnodes), 3 * (Solid.Nnodes)))
            self.K = np.zeros((3 * (Solid.Nnodes), 3 * (Solid.Nnodes)))
            self.deltaT = FRAC_DELT * endTime
            self.endTime = endTime
            self.E = E
            self.I = np.zeros((Solid.Nelements,1))
            self.A = np.zeros((Solid.Nelements,1))
            self.l = np.zeros((Solid.Nelements,1))
            self.RHO_S = RHO_S
            self.Fload = np.zeros((3*Solid.Nnodes,1))
            self.Fext_n = np.zeros((3*Solid.Nnodes,1))
            self.Fext_nPlus = np.zeros((3*Solid.Nnodes,1))
            self.Fint_n = np.zeros((3*Solid.Nnodes,1))
            self.Fint_nPlus = np.zeros((3*Solid.Nnodes,1))
            
            # Initial Displacements
            temp = 3 * Solid.fixedCounter
            self.U_n = np.zeros((3*Solid.Nnodes,1))
            self.Udot_n = np.zeros((3*Solid.Nnodes,1))
            self.UdotDot_n = np.zeros((3*Solid.Nnodes-temp,1))
            self.R_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
            self.R_n = np.zeros((3*Solid.Nnodes-temp,1))
            
            # Final Displacements
            self.U_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
            self.Udot_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
            self.UdotDot_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
            
            self.initU = np.zeros((3*Solid.Nnodes,1))
            self.initUdot = np.zeros((3*Solid.Nnodes,1))
        
    def elementTangentStiffnessMatrix(self, E, A, I, L, u_bar, theta_b1, theta_b2):
        """
        Calculates the element siffness matrix for bending and axial loads.

        Args:
            E (float): Young's Modulus of thesolid object.
            I (float): Element's area moment of inertia.
            A (float): Element's cross-sectional area.
            l (float): Length of the element.

        Return:
            k_e (float): NumPy 2D array of the element stiffness matrix.
        """
        EA = E * A
        EI = E * I
        
        # Part 1: Calculate matrix/array elements
        N  = EA * ((u_bar / L) + (1./15. * theta_b1**2) - (1./30. * theta_b1 * theta_b2) + (1./15. * theta_b2**2))
        M1 = EA * L *  ((u_bar / L) + (1./15. * theta_b1**2) - (1./30. * theta_b1 * theta_b2) + (1./15. * theta_b2**2)) * ((2./15. * theta_b1) - (1./30. * theta_b2)) + EI / L * (4. * theta_b1 + 2. * theta_b2)
        M2 = EA * L *  ((u_bar / L) + (1./15. * theta_b1**2) - (1./30. * theta_b1 * theta_b2) + (1./15. * theta_b2**2)) * ((2./15. * theta_b2) - (1./30. * theta_b1)) + EI / L * (2. * theta_b1 + 4. * theta_b2)
        KL11 = EA / L
        KL12 = EA * ((2./15 * theta_b1) - (1./30 * theta_b2))
        KL13 = EA * ((2./15. * theta_b2) - (1./30. * theta_b1))
        KL22 = EA * L * ((2./15. * theta_b1) - (1./30. * theta_b2))**2 + 2./15. * EA * L * ((u_bar / L) + (1./15. * theta_b1**2) - (1./30. * theta_b1 * theta_b2) + (1./15. * theta_b2**2)) + 4. * EI / L
        KL23 = EA * L * ((2./15. * theta_b2) - (1./30. * theta_b1)) * ((2./15. * theta_b1) - (1./30. * theta_b2)) - 1./30. * EA * L * ((u_bar / L) + (1./15. * theta_b1**2) - (1./30. * theta_b1 * theta_b2) + (1./15. * theta_b2**2)) + 2. * EI / L
        KL33 = EA * L * ((2./15. * theta_b2) - (1./30. * theta_b1))**2 + 2./15. * EA * L * ((u_bar / L) + (1./15. * theta_b1**2) - (1./30. * theta_b1 * theta_b2) + (1./15. * theta_b2**2)) + 4. * EI / L
        KL21 = np.copy(KL12)
        KL31 = np.copy(KL13)
        KL32 = np.copy(KL23)
        
        # Part 2: Assemble the matrix/array elements
        kl = np.array([[KL11, KL12, KL13],
                       [KL21, KL22, KL23],
                       [KL31, KL32, KL33]])
  
        return (kl, N, M1, M2)
        
    def elementMassMatrix(self, RHO_S, A, l, mType='consistent'):
        """
        Calculates the element mass matrix for bending and axial loads. This can
        return either a 'consistent' or 'lumped' mass matrix

        Args:
            RHO_S (float): Solid object's density.
            A (float): Element's cross-sectional area.
            l (float): Length of the element.
            mType (str): Type of Mass Matrix. must be 'consistent' or 'lumped'.

        Returns:
            m_e (float): NumPy 2D array of the element mass matrix.
            
        Raises:
            ValueError: If 'mType' is not defined as 'consistent' or 'lumped'.
        """
        
        if (mType == 'consistent'):
            C1 = RHO_S * A * l / 420.
            C2 = RHO_S * A * l / 6.
            m_e = np.array(
                           [[2.*C2,          0.,          0.,       1.*C2,         0.,          0.],
                            [    0.,    156.*C1,    22.*l*C1,          0.,     54.*C1,   -13.*l*C1],
                            [    0.,   22.*l*C1,  4.*l**2*C1,          0.,   13.*l*C1, -3.*l**2*C1],
                            [1.*C2,          0.,          0.,       2.*C2,         0.,          0.],
                            [    0.,     54.*C1,    13.*l*C1,          0.,    156.*C1,   -22.*l*C1],
                            [    0.,  -13.*l*C1, -3.*l**2*C1,          0.,  -22.*l*C1,  4.*l**2*C1]]
                          )
        elif (mType == 'lumped'):
            C1 = RHO_S * A * l / 420.
            C2 = RHO_S * A * l / 6.
            m_e = np.array(
                           [[2*C2,         0,          0,          C2,         0,          0],
                            [   0,        C1,          0,           0,         0,          0],
                            [   0,         0,          0,           0,         0,          0],
                            [  C2,         0,          0,        2*C2,         0,          0],
                            [   0,         0,          0,           0,        C1,          0],
                            [   0,         0,          0,           0,         0,          0]]
                          )           
        else:
            #TODO: Figure out how to throw an exception and hault exec.
            # An exception should be thrown and execuition haulted
            print 'ERROR: Invalid mass matrix type "%s"' % mType
            print 'Valid types are:'
            print '    "consistent"'
            print '    "lumped"'
           
        return m_e
        
    def elementTransformation(self, s, c, Ln):
        """
        Calculates the element connectivity matrix and vectors. This is used to
        formulate the global tangent stiffness matricies.
        
        Args:
            s (float): The current sine of alpha.
            c (float): Thecuurent cosine of alpha.
            Ln (float): The element deformed length.
            
        Returns:
            B (float): Element's local to global transformation matrix.
            r (float): Element's local to global transformation vector.
            z (float): Element's local to global transformation vector.
        """
        B = np.array([[-c,      -s, 0.,    c,     s, 0.],
                      [-s/Ln, c/Ln, 1., s/Ln, -c/Ln, 0.],
                      [-s/Ln, c/Ln, 0., s/Ln, -c/Ln, 1.]])
        r = np.array([[-c], [-s], [0.],  [c], [s], [0.]])
        z = np.array([[ s], [-c], [0.], [-s], [c], [0.]])
        
        return (B, r, z)
        
    def elementConnectivityMatrix(self, element, S, C):
        """
        Calculates the element connectivity matrix. This is used to formulate 
        the global mass and stiffness matricies.
        
        Args:
            elememnt (int): The current global element number.
            theta (float): The initial theta displacement.
            
        Returns:
            l_e (float): Element's local to global connectivity matrix.
        """
        element += 1
        l_e = np.zeros((6,3*(self.Nelements+1)))
        temp = np.array(
                        [   [ C,  S,  0,  0,  0,  0],
                            [-S,  C,  0,  0,  0,  0],
                            [ 0,  0,  1,  0,  0,  0],
                            [ 0,  0,  0,  C,  S,  0],
                            [ 0,  0,  0, -S,  C,  0],
                            [ 0,  0,  0,  0,  0,  1]    ]                       
                       )             
        l_e[:,3*element-2-1:5+3*element-2] = np.copy(temp)
        
        return l_e
        
    def globalToLocalDisp(self, L0, Ln, theta1, theta2, alpha):
        """
        Calculates the element displacement quantities from global values.
        
        Args:
            L0 (float): Element undeformed length.
            Ln (float): Element deformed length.
            theta1 (float): Element global theta displacement for local node 1.
            theta2 (float): Element global theta displacement for local node 2.
            alpha (float): Rigid-body rotation undergone by the element.
            
        Returns:
            u_bar (float): Local change in element length.
            theta_b1 (float): Deformational rotation at local node 1.
            theta_b2 (float): Deformational rotation at local node 2.
        """
        u_bar = Ln - L0
        theta_b1 = theta1 - alpha
        theta_b2 = theta2 - alpha
        
        return (u_bar, theta_b1, theta_b2)
        
    def globalMatricies(self, Solid, mType, NRT, alphap, sna0, U):
        
        K = np.copy(self.K)
        M = np.copy(self.M)
        Fint = np.copy(self.Fint_n)
        K.fill(0.)
        M.fill(0.)
        Fint.fill(0.)
        
        for i in xrange(self.Nelements):
            # Get element coordinates
            x1 = np.copy(Solid.nodes_0[i  ,0])
            x2 = np.copy(Solid.nodes_0[i+1,0])
            z1 = np.copy(Solid.nodes_0[i  ,1])
            z2 = np.copy(Solid.nodes_0[i+1,1])
            
            # Get element displacements
            j = 3*i
            u1 = np.copy(U[j  ,0])
            u2 = np.copy(U[j+3,0])
            w1 = np.copy(U[j+1,0])
            w2 = np.copy(U[j+4,0])
            theta1 = np.copy(U[j+2,0])
            theta2 = np.copy(U[j+5,0])
            
            # Get element properties
            L0 = np.copy(self.l[i])
            Ln = np.sqrt((x2 + u2 - x1 - u1)**2 + (z2 + w2 - z1 - w1)**2)
            
            # Get element directionaliies
            c0 = (x2 - x1) / L0
            s0 = (z2 - z1) / L0
            c = (x2 + u2 - x1 - u1) / Ln
            s = (z2 + w2 - z1 - w1) / Ln
            sna = c0 * s - s0 * c
            ca  = c0 * c + s0 * s
            L = np.copy(L0)
            
            # Determine sine(alpha0)
            if (NRT == 1):
                alpha = 0.
            else:
                if (sna0[i] >= 0.):
                    alpha = np.arctan2(sna,ca)
                    if (alpha < 0.):
                        alpha = 2. * np.pi + alpha
                    ad = np.absolute(np.absolute(alpha) - np.absolute(alphap[i]))
                    if (ad > 2.):
                        alpha = np.arctan2(sna,ca)
                        sna0[i] = -1 * sna0[i]
                else:
                    alpha = np.arctan2(sna,ca)
                    if (alpha > 0.):
                        alpha = -2 * np.pi + alpha
                    ad = np.absolute(np.absolute(alpha) - np.absolute(alphap[i]))
                    if (ad > 2.):
                        alpha = np.arctan2(sna,ca)
                        sna0[i] = -1 * sna0[i]    
            alphap[i] = np.copy(alpha)
                
            # Calculate local node displacements
            (u_bar, theta_b1, theta_b2) = self.globalToLocalDisp(L0, Ln, theta1, theta2, alpha)
            
            # Build transformation matrix and vectors
            (B, r, z) = self.elementTransformation(s, c, Ln)
            l_e = self.elementConnectivityMatrix(i, s, c)
            
            # Calculate element modified internal forces forces to avoid membrane locking
            (kl, N, M1, M2) = self.elementTangentStiffnessMatrix(self.E, self.A[i], self.I[i], L, u_bar, theta_b1, theta_b2)
            m_e = self.elementMassMatrix(self.RHO_S, self.A[i], L, mType)
            fint1 = np.array([[N], [M1], [M2]])
            
            # Transform to element global reference frame
            fint = np.dot(np.transpose(B), fint1)
            Ktan1 = np.dot(np.dot(np.transpose(B),kl), B)
            Ktan2 = np.dot(z, np.transpose(z)) * N / Ln
            Ktan3 = (np.dot(r,np.transpose(z)) + np.dot(z,np.transpose(r))) * (M1 + M2) / Ln**2
            Ktan = Ktan1 + Ktan2 + Ktan3
            Mg = np.dot(np.dot(np.transpose(l_e),m_e), l_e)
                           
            # Add element matricies to the global matricies
            j = 3*i
            Fint[j:j+6,:]  = Fint[j:j+6,:]  + np.copy(fint)
            K[j:j+6,j:j+6] = K[j:j+6,j:j+6] + np.copy(Ktan)
            M = M + np.copy(Mg)
        
        return (M, K, Fint, sna0, x1, x2, z1, z2)
        
        
    def steadySolve(self, Body, Solid, nsteps):
        """
        Solves a steady finite element system of equations.
        
        Args: 
            Body (object): A body object created from the swimmer class.
            Solid (object): A solid object created from the solid class.
            outerCorr (int): Current FSI subiteration number.
            mType (str): Type of Mass Matrix. must be 'consistent' or 'lumped'.

        Raises:
            ValueError: If 'method' is not defined as 'HHT', 'NEWMARK', or 'TRAPEZOIDAL'.           
        """
        # Initialize varriables
        alphap = np.zeros(self.Nelements)
        F = np.zeros((3*Solid.Nnodes,1))
        temp = 3 * Solid.fixedCounter
        sna0 = np.zeros(self.Nelements)
        NRT = 0
        
        # Create local variables
        U = np.copy(self.U_n)
        K = np.copy(self.K)
        Fint = np.copy(self.Fint_n)
        Fext = np.copy(self.Fload)
        
        # Define the load increment
        Finc = Fext / nsteps
        
        
        for step in xrange(nsteps):
            # Update the load
            F = F + Finc
            
            for innerCorr in xrange(201):
                # Update equilibrium iteration counter
                NRT = NRT + 1
                K.fill(0.)
                Fint.fill(0.)
                
                for i in xrange(self.Nelements):
                    # Get element coordinates
                    x1 = np.copy(Solid.nodes_0[i  ,0])
                    x2 = np.copy(Solid.nodes_0[i+1,0])
                    z1 = np.copy(Solid.nodes_0[i  ,1])
                    z2 = np.copy(Solid.nodes_0[i+1,1])
                    
                    # Get element displacements
                    j = 3*i
                    u1 = np.copy(U[j  ,0])
                    u2 = np.copy(U[j+3,0])
                    w1 = np.copy(U[j+1,0])
                    w2 = np.copy(U[j+4,0])
                    theta1 = np.copy(U[j+2,0])
                    theta2 = np.copy(U[j+5,0])
                    
                    # Get element properties
                    L0 = np.copy(self.l[i])
                    Ln = np.sqrt((x2 + u2 - x1 - u1)**2 + (z2 + w2 - z1 - w1)**2)
                    
                    # Get element directionaliies
                    c0 = (x2 - x1) / L0
                    s0 = (z2 - z1) / L0
                    c = (x2 + u2 - x1 - u1) / Ln
                    s = (z2 + w2 - z1 - w1) / Ln
                    sna = c0 * s - s0 * c
                    ca  = c0 * c + s0 * s
                    L = np.copy(L0)
                    
                    # Determine sine(alpha0)
                    if (NRT == 1):
                        alpha = 0.
                    else:
                        if (sna0[i] >= 0.):
                            alpha = np.arctan2(sna,ca)
                            if (alpha < 0.):
                                alpha = 2. * np.pi + alpha
                            ad = np.absolute(np.absolute(alpha) - np.absolute(alphap[i]))
                            if (ad > 2.):
                                alpha = np.arctan2(sna,ca)
                                sna0[i] = -1 * sna0[i]
                        else:
                            alpha = np.arctan2(sna,ca)
                            if (alpha > 0.):
                                alpha = -2 * np.pi + alpha
                            ad = np.absolute(np.absolute(alpha) - np.absolute(alphap[i]))
                            if (ad > 2.):
                                alpha = np.arctan2(sna,ca)
                                sna0[i] = -1 * sna0[i]    
                    alphap[i] = np.copy(alpha)
                        
                    # Calculate local node displacements
                    (u_bar, theta_b1, theta_b2) = self.globalToLocalDisp(L0, Ln, theta1, theta2, alpha)
                    
                    # Build transformation matrix and vectors
                    (B, r, z) = self.elementTransformation(s, c, Ln)
                    
                    # Calculate element modified internal forces forces to avoid membrane locking
                    (kl, N, M1, M2) = self.elementTangentStiffnessMatrix(self.E, self.A[i], self.I[i], L, u_bar, theta_b1, theta_b2)
                    fint1 = np.array([[N], [M1], [M2]])
                    
                    # Transform to element global reference frame
                    fint = np.dot(np.transpose(B), fint1)
                    Ktan1 = np.dot(np.dot(np.transpose(B),kl), B)
                    Ktan2 = np.dot(z, np.transpose(z)) * N / Ln
                    Ktan3 = (np.dot(r,np.transpose(z)) + np.dot(z,np.transpose(r))) * (M1 + M2) / Ln**2
                    Ktan = Ktan1 + Ktan2 + Ktan3
                                   
                    # Add element matricies to the global matricies
                    j = 3*i
                    Fint[j:j+6,:]  = Fint[j:j+6,:]  + np.copy(fint)
                    K[j:j+6,j:j+6] = K[j:j+6,j:j+6] + np.copy(Ktan)
                
                # Calculate load residual
                delR = (F - Fint)
                
                # Solve for incremental displacement
                DeltaU = spla.solve(K[temp:,temp:], delR[temp:,:])
                
                # Update nodal displacements
                U[temp:,:] = U[temp:,:] + DeltaU
                
                if (NRT == 1):
                    for i in xrange(self.Nelements):
                        j = 3*i
                        sna0[i] = -((U[j+3,0] - U[j,0]) * (z2 - z1) - (U[j+4,0] - U[j+1,0]) * (x2 - x1))

                # Check if energy is conserved within tolerance
                du = np.linalg.norm(DeltaU, ord=2)
                if du < 1e-6 or innerCorr == 200:
                    if innerCorr == 200:
                        print 'ERROR! Max iterations reached in structural solve'
    #                    raise ValueError('Maximum structural iterations reached')
                    break

        # Store the final displacements
#        self.U_n = np.copy(U)
        self.U_nPlus = np.copy(U[temp:,:])
        
    def dynamicSolve(self, Body, Solid, outerCorr, mType='consistent'):
        """
        Solves a steady finite element system of equations.
        
        Args: 
            Body (object): A body object created from the swimmer class.
            Solid (object): A solid object created from the solid class.
            outerCorr (int): Current FSI subiteration number.
            mType (str): Type of Mass Matrix. must be 'consistent' or 'lumped'.

        Raises:
            ValueError: If 'method' is not defined as 'HHT', 'NEWMARK', or 'TRAPEZOIDAL'.           
        """
        # Initialize varriables
        alphap = np.zeros(self.Nelements)
        F = np.zeros((3*Solid.Nnodes,1))
        temp = 3 * Solid.fixedCounter
        sna0 = np.zeros(self.Nelements)
        NRT = 0
        beta = 0.25
        gamma = 0.5
        
        # Create local variables
        U = np.copy(self.U_n)
        U_n = np.copy(self.U_n)
        Udot_n = np.copy(self.Udot_n)
        UdotDot_n = np.copy(self.UdotDot_n)  
        
        

        dt = self.deltaT
        
        # Set the force acting at the begining and end of the time integration
        if (outerCorr == 1):
            self.Fext_n = np.copy(self.Fext_nPlus)
            self.R_n = np.copy(self.R_nPlus)
        Fext_n = np.copy(self.Fext_n)
        Fext_nPlus = np.copy(self.Fload)
        F = Fext_nPlus
        
        for innerCorr in xrange(1001):
            # Update equilibrium iteration counter
            NRT = NRT + 1

            # Generate global Mass, Tangent Stiffness, and Internal Force Matricies/Vectors
            (M, K, Fint, sna0, x1, x2, z1, z2) = self.globalMatricies(Solid, mType, NRT, alphap, sna0, U)

            # Calculate load residual
            delR = (F - Fint)
            
            # Build linear system of equations
            # Part 1: The left hand coefficient matrix
            A = 1./(beta * dt**2) * M[temp:,temp:] + K[temp:,temp:]
            
            # Part 2: The RHS
            c1 = Udot_n[temp:,:] / beta / dt + (1. / (2. * beta) - 1.) * UdotDot_n
            b = delR[temp:,:] + np.dot(self.M[temp:,temp:], c1)
            
            # Solve for incremental displacement
            DeltaU = spla.solve(A, b)
            
            # Update nodal displacements
            U_nPlus = U[temp:,:] + DeltaU
            Udot_nPlus    = gamma / beta / dt * DeltaU - (gamma / beta - 1.) * Udot_n[temp:,:] - dt * ((gamma / 2. / beta) - 1.) * UdotDot_n
            UdotDot_nPlus = 1. / beta / dt**2 * (DeltaU - dt * Udot_n[temp:,:]) - ((1. / 2. / beta) - 1.) * UdotDot_n
            
            # Update internal force
            Fint_nPlus = Fint[temp:,:] + np.dot(K[temp:,temp:], DeltaU)
            
            # Estimate Error
            Rerr = F[temp:,:] - np.dot(self.M[temp:,temp:], UdotDot_nPlus) - Fint_nPlus
            
            if (NRT == 1):
                for i in xrange(self.Nelements):
                    j = 3*i
                    U2 = np.copy(U)
                    U2.fill(0.)
                    U2[temp:,:] = np.copy(U_nPlus)
                    sna0[i] = -((U2[j+3,0] - U2[j,0]) * (z2 - z1) - (U2[j+4,0] - U2[j+1,0]) * (x2 - x1))

            # Check if energy is conserved within tolerance
            err = np.linalg.norm(Rerr, ord=2)
            Werr = np.linalg.norm(dt * np.dot(np.transpose(Udot_nPlus), Rerr), ord=2)
            if Werr < 1e-6 or innerCorr == 1000:
                if innerCorr == 1000:
                    print '+-----------------------------------------------------------------------------+'
                    print '| WARNING! Maximum iterations reached in structural solve                     |'
                    print '+-----------------------------------------------------------------------------+'
                break
            else:
                U[temp:,:] = np.copy(U_nPlus)
#                F[temp:,:] = F[temp:,:] + Rerr

        # Store the final displacements
        self.U_nPlus = np.copy(U_nPlus)
        self.Udot_nPlus = np.copy(Udot_nPlus)
        self.UdotDot_nPlus = np.copy(UdotDot_nPlus)
        
    def spring_solve(self):
        """Solves passive pitching of the leading edge.
        
        Args:
            Body (obj): A Body object that the fluid acts on
            Solid (obj): A Solid mesh object
            outerCorr (int): The FSI iteration number
        """
        dt = self.deltaT / 1000.
        I = self.I
        kappa_1 = self.kappa_1
        kappa_2 = self.kappa_2
        zeta = self.zeta
        Nf = self.Nf
        Ni = self.Ni
        theta_n           = np.copy(self.theta_n)
        thetaDot_n        = np.copy(self.thetaDot_n)
        thetaDotDot_n     = np.copy(self.thetaDotDot_n)
        theta_nPlus       = np.copy(self.theta_nPlus)
        thetaDot_nPlus    = np.copy(self.thetaDot_nPlus)
        thetaDotDot_nPlus = np.copy(self.thetaDotDot_nPlus)
        
        # Build linear system of equations
        # Part 1: The left hand coefficient matrix
        A = kappa_1 + I * (2. / dt)**2 + zeta * (2. / dt)
        
        for i in xrange(1000):
            # Part 2: The RHS
            b = I * ((2. / dt)**2 * theta_n + (4. / dt) * thetaDot_n + thetaDotDot_n) + zeta * ((2. / dt) * theta_n + thetaDot_n) + Nf + Ni
            
            # Solve for pitching angle theta
            if (kappa_2 == 0.):
                theta_nPlus = b / A
            else:
                # Since the spring is cubic, there is a possible of three 
                # solutions to satisfy the equations; however, the problem 
                # should be posed so that there is only one real valued root. 
                # The following is the solution to that real valued root.
                term1 = (np.sqrt(3.) * np.sqrt(4. * A**3 * kappa_2**3 + 27. * kappa_2**4 * b**2) + 9.*kappa_2**2 * b)**(1./3.)
                term2 = term1 / (2.**(1. / 3.) * 3.**(2. / 3.) * kappa_2)
                if (term1 == 0.):
                    theta_nPlus = b / A
                else:
                    term3 = ((2. / 3.)**(1. / 3.) * A) / term1
                    theta_nPlus = term2 - term3
            
            # Update angular velocity and acceleration
            thetaDotDot_nPlus = (2. / dt)**2 * (theta_nPlus - theta_n) - (4. / dt) * thetaDot_n - thetaDotDot_n
            thetaDot_nPlus = thetaDot_n + 0.5 * (thetaDotDot_n + thetaDotDot_nPlus) * dt
            
            # Update Inital angular position, velocity, and acceleration
            theta_n       = np.copy(theta_nPlus)
            thetaDot_n    = np.copy(thetaDot_nPlus)
            thetaDotDot_n = np.copy(thetaDotDot_nPlus)
            
        
        # Store final values
        self.theta_nPlus = theta_nPlus
        self.thetaDot_nPlus = thetaDot_nPlus
        self.thetaDotDot_nPlus = thetaDotDot_nPlus
