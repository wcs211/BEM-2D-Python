#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np

class PyFEA(object):
    def __init__(self, Solid, FRAC_DELT, endTime, E, RHO_S):
        """
        Iniitalizes object related variables needed for other class methods.
        
        Keyword arguments:

        """
        
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
        
#       Initial Displacements
        temp = 3 * Solid.fixedCounter
        self.U_n = np.zeros((3*Solid.Nnodes,1))
        self.Udot_n = np.zeros((3*Solid.Nnodes,1))
        self.UdotDot_n = np.zeros((3*Solid.Nnodes-temp,1))
        
#       Final Displacements
        self.U_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
        self.Udot_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
        self.UdotDot_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
        
        self.initU = np.zeros((3*Solid.Nnodes,1))
        self.initUdot = np.zeros((3*Solid.Nnodes,1))
        
    def elementStiffnessMatrix(self, E, I, A, l):
        """
        Returns the element siffness matrix for bending and axial loads

        Keyword arguments:
        E -- Young's Modulus
        I -- Area Moment of Inertia
        A -- Cross-sectional Area
        l -- Undeformed Element Length
        """
        C1 = (E * A / l)
        C2 = (E * I / l**3)
        k_e = np.array(
                       [[1.*C1,          0,           0,         -1.*C1,           0,           0],
                        [0,       12*C2,      6*l*C2,           0,      -12*C2,      6*l*C2],
                        [0,      6*l*C2,   4*l**2*C2,           0,     -6*l*C2,   2*l**2*C2],
                        [-1.*C1,         0,           0,          1.*C1,           0,           0],
                        [0,      -12*C2,     -6*l*C2,           0,       12*C2,     -6*l*C2],
                        [0,      6*l*C2,   2*l**2*C2,           0,     -6*l*C2,   4*l**2*C2]]        
                      )
        return k_e
        
    def elementMassMatrix(self, RHO_S, A, l, mType):
        """
        Returns the element mass matrix for bending and axial loads. This can
        return either a 'consistent' or 'lumped' mass matrix

        Keyword arguments:
        RHO_S -- Young's Modulus
        A -- Cross-sectional Area
        l -- Undeformed Element Length
        mType -- Type of Mass Matrix. must be 'consistent' or 'lumped'
        """
        
        if (mType == 'consistent'):
            C1 = RHO_S * A * l / 420
            C2 = RHO_S * A * l / 6
            m_e = np.array(
                           [[2*C2,         0,          0,          1.*C2,         0,          0],
                            [   0,    156*C1,    22*l*C1,           0,     54*C1,   -13*l*C1],
                            [   0,   22*l*C1,  4*l**2*C1,           0,   13*l*C1, -3*l**2*C1],
                            [  1.*C2,         0,          0,        2*C2,         0,          0],
                            [   0,     54*C1,    13*l*C1,           0,    156*C1,   -22*l*C1],
                            [   0,  -13*l*C1, -3*l**2*C1,           0,  -22*l*C1,  4*l**2*C1]]
                          )
        elif (mType == 'lumped'):
            C1 = RHO_S * A * l / 420
            C2 = RHO_S * A * l / 6
            m_e = np.array(
                           [[2*C2,         0,          0,          C2,         0,          0],
                            [   0,        C1,          0,           0,         0,          0],
                            [   0,         0,          0,           0,         0,          0],
                            [  C2,         0,          0,        2*C2,         0,          0],
                            [   0,         0,          0,           0,        C1,          0],
                            [   0,         0,          0,           0,         0,          0]]
                          )           
        else:
            # An exception should be thrown and execuition haulted
            print 'ERROR: Invalid mass matrix type "%s"' % mType
            print 'Valid types are:'
            print '    "Consistent"'
            print '    "lumped"'
            
        return m_e
        
    def elementConnectivityMatrix(self, element, theta):
        """
        Solves a unsteady finite element system of equations.
        
        Keyword arguments:
        element -- The current global element number
        theta -- The initial theta displacement
        """
        element += 1
        l_e = np.zeros((6,3*(self.Nelements+1)))
        C = np.cos(theta)
        S = np.sin(theta)
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
        
    def HHT(self, alpha, beta, gamma, Fext_n, Fext_nPlus, fixedNodes, U_n, Udot_n, UdotDot_n):
        """
        Solves a dynamic system of equations using the Hilber-Hughes-Taylor 
        (HHT) Method. This is a transient, implicit method with numerical
        dissipation in the high frequency domain. This method has second-order
        accuracy.
        
        Keyword arguments:
        alpha -- Integraton constant
        beta -- Integration constant
        gamma -- Integration constant
        Fext_n -- Force exterted at the begining of the time-step
        Fext_nPlus -- Force exterted at the end of the time-step
        fixedNodes -- Count of the number of nodes with a no displacement condition
        """
        
        temp = 3 * fixedNodes
        
        # Form the 'A' matrix
        A = (self.M[temp:,temp:] + beta * self.deltaT**2 * (1 - alpha) * self.K[temp:,temp:]) / (beta * self.deltaT**2)
               
        # Form the 'B' matrix
        B = (1 - alpha) * Fext_nPlus[temp:,:] + 1 / (beta * self.deltaT**2) * \
            np.dot(self.M[temp:,temp:], U_n[temp:,:] + self.deltaT * Udot_n[temp:,:] + self.deltaT**2 * \
            (0.5 - beta) * UdotDot_n) + \
            alpha * Fext_n[temp:,:] - alpha * np.dot(self.K[temp:,temp:], U_n[temp:,:])
            
        # Solve the system to get the displacements
        U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        UdotDot_nPlus = (U_nPlus - (U_n[temp:,:] + self.deltaT * Udot_n[temp:,:] + self.deltaT**2 * (0.5 - beta) * UdotDot_n)) / (beta * self.deltaT**2)
        
        # Solve for the velocities
        Udot_nPlus = (Udot_n[temp:,:] + self.deltaT * (1 - gamma) * UdotDot_n) + gamma * self.deltaT * UdotDot_nPlus
        
        return (U_nPlus, Udot_nPlus, UdotDot_nPlus)
        
    def NEWMARK(self, beta, gamma, Fext_n, Fext_nPlus, fixedNodes):
        """
        Solves a dynamic system of equations using the NEWMARK Method.
        This is a transient, implicit method with numerical dissipation in the 
        high frequency domain. This method has first-order accuracy.
        
        Keyword arguments:
        beta -- Integration constant
        gamma -- Integration constant
        Fext_n -- Force exterted at the begining of the time-step
        Fext_nPlus -- Force exterted at the end of the time-step
        fixedNodes -- Count of the number of nodes with a no displacement condition
        """
        
        # Form the 'A' matrix
        A = (self.M + beta * self.deltaT**2 * self.K) / (beta * self.deltaT**2)
        
        # Form the 'B' matrix
        B = Fext_nPlus + 1 / (beta * self.deltaT**2) * self.M * (self.U_n + \
            self.deltaT * self.Udot_n + \
            self.deltaT**2 * (0.5 - beta) * self.UdotDot_n)
            
        # Solve the system to get the displacements
        self.U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        self.UdotDot_nPlus = (self.U_nPlus - (self.U_n + self.deltaT * self.Udot_n + self.deltaT**2 * (0.5 - beta) * self.UdotDot_n)) / (beta * self.deltaT**2)
        
        # Solve for the velocities
        self.Udot_nPlus = (self.Udot_n + self.deltaT * (1 - gamma) * self.UdotDot_n) + gamma * self.deltaT * self.UdotDot_nPlus
        
    def TRAPEZOIDAL(self, Fext_nPlus, fixedNodes):
        """
        Solves for the system dynamics using the trapezoidal rule.
        
        Keyword arguments:
        Fext_nPlus -- Force exterted at the end of the time-step
        fixedNodes -- Count of the number of nodes with a no displacement condition
        """
        
        # Form the 'A' matrix
        A = (self.K + (2 / self.deltaT)**2 * self.M)
        
        # Form the 'B' matrix
        B = (Fext_nPlus + self.M * ((2 / self.deltaT)**2 * self.U_n + \
            (4 / self.deltaT) * self.Udot_n + self.UdotDot_n))
            
        # Solve the system to get the displacements
        self.U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        self.UdotDot_nPlus = 2 * (self.Udot_nPlus - self.Udot_n) / self.deltaT - self.UdotDot_n
        
        # Solve for the velocities
        self.Udot_nPlus = 2 * (self.U_nPlus - self.U_n) / self.deltaT - self.Udot_n
        
    def solve(self, Body, Solid, outerCorr, t, TSTEP, mType, method, alpha, beta, gamma):
        """
        Solves a unsteady finite element system of equations.
        
        Keyword arguments:
        mType -- Type of Mass Matrix. must be 'consistent' or 'lumped'
        """
        # Determine current pitching angle
#        theta = Body.MP.THETA_MAX * np.sin(2 * np.pi * Body.MP.F * (t + TSTEP) + Body.MP.PHI)
#        theta = 5*np.pi/180*np.tanh(t)
        U_n = np.copy(self.U_n)
        Udot_n = np.copy(self.Udot_n)
        UdotDot_n = np.copy(self.UdotDot_n)        
        
        # Reset mass and stiffness matrix to include all nodes
        self.M = 0
        self.K = 0
        
        # Assemble global mass and stiffness matricies
        for i in xrange(self.Nelements):
            # Clear/initialize old element matrix
            k_e = np.zeros((6, 6))
            m_e = np.zeros((6, 6))
            l_e = np.zeros((6, 2 * (self.Nelements + 1)))
            
            # Determine element stiffness, mass, and connectivity matricies
            k_e = self.elementStiffnessMatrix(self.E, self.I[i], self.A[i], self.l[i])
            m_e = self.elementMassMatrix(self.RHO_S, self.A[i], self.l[i], mType)
            l_e = self.elementConnectivityMatrix(i, -U_n[3*i-1,0])
            
            # Add element matricies to the global matricies
            self.M = self.M + np.dot(np.dot(np.transpose(l_e), m_e), l_e)
            self.K = self.K + np.dot(np.dot(np.transpose(l_e), k_e), l_e)
            
        # Set the zero displacement constraints
        temp = 3 * Solid.fixedCounter
        
        # Solve for the initial acceleration matrix
        Fext_n = np.copy(self.Fload)
        RHS = Fext_n[temp:,:] - np.dot(self.K[temp:, temp:], U_n[temp:])
        UdotDot_n = np.linalg.solve(self.M[temp:,temp:], RHS)
        
        # March through time until the total simulated time has elapsed
        j = np.size(np.arange(self.deltaT,self.endTime+self.deltaT,self.deltaT))
        for i in xrange(j):
            Fext_nPlus = np.copy(Fext_n)
            if (method == 'HHT'):
                (U_nPlus, Udot_nPlus, UdotDot_nPlus) = self.HHT(alpha, beta, gamma, Fext_n, Fext_nPlus, Solid.fixedCounter, U_n, Udot_n, UdotDot_n)
            elif (method == 'NEWMARK'):
                self.NEWMARK(beta, gamma, Fext_n, Fext_nPlus, Solid.fixedCounter)
            elif (method == ' TRAPEZOIDAL'):
                self.TRAPEZOIDAL(Fext_nPlus, Solid.fixedCounter)
            else:
                # Throw exception and hault execuition
                print 'ERROR! Invalid integration scheme "%s".' % method
                print 'Valid schemes are:'
                print '    HHT'
                print '    NEWMARK'
                print '    TRAPEZOIDAL'
            if (i != j):
                U_n[temp:,:] = np.copy(U_nPlus)
                Udot_n[temp:,:] = np.copy(Udot_nPlus)
                UdotDot_n = np.copy(UdotDot_nPlus)

        self.U_nPlus = np.copy(U_nPlus)
        self.Udot_nPlus = np.copy(Udot_nPlus)
        self.UdotDot_nPlus = np.copy(UdotDot_nPlus)