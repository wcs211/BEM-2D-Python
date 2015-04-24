#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np

class PyFEA(object):
    def __init__(self, Nelements, fracDeltaT, endTime, E, I, A, l, rho, Fload, U_n, Udot_n):
        """
        Iniitalizes object related variables needed for other class methods.
        
        Keyword arguments:

        """
        
        self.Nelements = Nelements
        self.M = np.zeros((3 * (self.Nelements + 1)))
        self.K = np.zeros((3 * (self.Nelements + 1)))
        self.deltaT = fracDeltaT * endTime
        self.endTime = endTime
        self.E = E
        self.I = I
        self.A = A
        self.l = l
        self.rho = rho
        self.Fload = Fload
        self.U_n = U_n
        self.Udot_n = Udot_n
        self.UdotDot_n = np.zeros(())
        
    def elementStiffnessMatrix(E, I, A, l):
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
                       [[C1,          0,           0,         -C1,           0,           0],
                        [0,       12*C2,      6*l*C2,           0,      -12*C2,      6*l*C2],
                        [0,      6*l*C2,   4*l**2*C2,           0,     -6*l*C2,   2*l**2*C2],
                        [-C1,         0,           0,          C1,           0,           0],
                        [0,      -12*C2,     -6*l*C2,           0,       12*C2,     -6*l*C2],
                        [0,      6*l*C2,   2*l**2*C2,           0,     -6*l*C2,   4*l**2*C2]]        
                      )
        return k_e
        
    def elementMassMatrix(rho, A, l, mType):
        """
        Returns the element mass matrix for bending and axial loads. This can
        return either a 'consistent' or 'lumped' mass matrix

        Keyword arguments:
        rho -- Young's Modulus
        A -- Cross-sectional Area
        l -- Undeformed Element Length
        mType -- Type of Mass Matrix. must be 'consistent' or 'lumped'
        """
        
        if (mType == 'consistent'):
            C1 = rho * A * l / 420
            C2 = rho * A * l / 6
            m_e = np.array(
                           [[2*C2,         0,          0,          C2,         0,          0],
                            [   0,    156*C1,    22*l*C1,           0,     54*C1,   -13*l*C1],
                            [   0,   22*l*C1,   4*l^2*C1,           0,   13*l*C1,  -3*l^2*C1],
                            [  C2,         0,          0,        2*C2,         0,          0],
                            [   0,     54*C1,    13*l*C1,           0,    156*C1,   -22*l*C1],
                            [   0,  -13*l*C1,  -3*l^2*C1,           0,  -22*l*C1,   4*l^2*C1]]
                          )
        elif (mType == 'lumped'):
            C1 = rho * A * l / 420
            C2 = rho * A * l / 6
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
                       
        l_e[:,3*element-2:5+3*element-2] = temp
        
        return l_e
        
    def HHT(self, alpha, beta, gamma, Fext_n, Fext_nPlus, fixedNodes):
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
        
        # Form the 'A' matrix
        A = np.zeros(self.Nelements - 3 * fixedNodes)
        A = (self.M + beta * self.deltaT**2 * (1 - alpha) * self.K) / (beta * self.deltaT**2)
        
        # Form the 'B' matrix
        B = np.zeros((self.Nelements - 3 * fixedNodes, 1))
        B = (1 - alpha) * Fext_nPlus + 1 / (beta * self.deltaT**2) * self.M * \
            (self.U_n + self.deltaT * self.Udot_n + self.deltaT**2 * (0.5 - beta) * \
            self.UdotDot_n) + alpha * Fext_n - alpha * self.K * self.U_n
            
        # Solve the system to get the displacements
        U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        UdotDot_nPlus = (U_nPlus - (self.U_n + self.deltaT * self.Udot_n + self.deltaT**2 * (0.5 - beta) * self.UdotDot_n)) / (beta * self.deltaT**2)
        
        # Solve for the velocities
        Udot_nPlus = (self.Udot_n + self.deltaT * (1 - gamma) * self.UdotDot_n) + gamma * self.deltaT * UdotDot_nPlus
        
        return U_nPlus, Udot_nPlus, UdotDot_nPlus
        
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
        A = np.zeros(self.Nelements - 3 * fixedNodes)
        A = (self.M + beta * self.deltaT**2 * self.K) / (beta * self.deltaT**2)
        
        # Form the 'B' matrix
        B = np.zeros((self.Nelements - 3 * fixedNodes, 1))
        B = Fext_nPlus + 1 / (beta * self.deltaT**2) * self.M * (self.U_n + \
            self.deltaT * self.Udot_n + \
            self.deltaT**2 * (0.5 - beta) * self.UdotDot_n)
            
        # Solve the system to get the displacements
        U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        UdotDot_nPlus = (U_nPlus - (self.U_n + self.deltaT * self.Udot_n + self.deltaT**2 * (0.5 - beta) * self.UdotDot_n)) / (beta * self.deltaT**2)
        
        # Solve for the velocities
        Udot_nPlus = (self.Udot_n + self.deltaT * (1 - gamma) * self.UdotDot_n) + gamma * self.deltaT * UdotDot_nPlus
        
        return U_nPlus, Udot_nPlus, UdotDot_nPlus
        
    def TRAPEZOIDAL(self, Fext_nPlus, fixedNodes):
        """
        Solves for the system dynamics using the trapezoidal rule.
        
        Keyword arguments:
        Fext_nPlus -- Force exterted at the end of the time-step
        fixedNodes -- Count of the number of nodes with a no displacement condition
        """
        
        # Form the 'A' matrix
        A = np.zeros(self.Nelements - 3 * fixedNodes)
        A = (self.K + (2 / self.deltaT)**2 * self.M)
        
        # Form the 'B' matrix
        B = np.zeros((self.Nelements - 3 * fixedNodes, 1))
        B = (Fext_nPlus + self.M * ((2 / self.deltaT)**2 * self.U_n + \
            (4 / self.deltaT) * self.Udot_n + self.UdotDot_n))
            
        # Solve the system to get the displacements
        U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        UdotDot_nPlus = 2 * (self.Udot_nPlus - self.Udot_n) / self.deltaT - self.UdotDot_n
        
        # Solve for the velocities
        Udot_nPlus = 2 * (self.U_nPlus - self.U_n) / self.deltaT - self.Udot_n
        
        return U_nPlus, Udot_nPlus, UdotDot_nPlus
        
    def solve(self, mType, method, theta, fixedNodes, alpha, beta, gamma):
        """
        Solves a unsteady finite element system of equations.
        
        Keyword arguments:
        mType -- Type of Mass Matrix. must be 'consistent' or 'lumped'
        """
        # Assemble global mass and stiffness matricies
        for i in xrange(self.Nelements):
            # Clear/initialize old element matrix
            k_e = np.zeros((6, 6))
            m_e = np.zeros((6, 6))
            l_e = np.zeros((6, 2 * (self.Nelements + 1)))
            
            # Determine element stiffness, mass, and connectivity matricies
            k_e = self.elementStiffnessMatrix(self.E[i], self.I[i], self.A[i], self.l[i])
            m_e = self.elementMassMatrix(self.rho, self.A[i], self.l[i], mType)
            l_e = self.elementConnectivityMatrix(i, theta)
            
            # Add element matricies to the global matricies
            self.M = self.M + np.transpose(l_e) * m_e * l_e
            self.K = self.K + np.transpose(l_e) * k_e * l_e
            
        # Set the zero displacement constraints
        temp = 3 * fixedNodes
        self.M = np.delete(self.M, np.s_[0:temp,0:temp])
        self.K = np.delete(self.K, np.s_[0:temp,0:temp])
        self.Fload = np.delete(self.Fload, np.s_[0:temp])
        self.U_n = np.delete(self.U_n, np.s_[0:temp])
        self.Udot = np.delete(self.Udot_n, np.s_[0:temp])
        
        # Solve for the initial acceleration matrix
        Fext_n = self.Fload
        self.UdotDot_n = np.linalg.solve(self.M, Fext_n)
        
        # March through time until the total simulated time has elapsed
        for i in xrange(np.size(np.s_[self.deltaT:self.endTime:self.deltaT])):
            Fext_nPlus = Fext_n
            if (method == 'HHT'):
                self.HHT(alpha, beta, gamma, Fext_n, Fext_nPlus, fixedNodes)
            elif (method == 'NEWMARK'):
                self.NEWMARK(beta, gamma, Fext_n, Fext_nPlus, fixedNodes)
            elif (method == ' TRAPEZOIDAL'):
                self.TRAPEZOIDAL(Fext_nPlus, fixedNodes)
            else:
                # Throw exception and hault execuition
                print 'ERROR! Invalid integration scheme "%s".' % method
                print 'Valid schemes are:'
                print '    HHT'
                print '    NEWMARK'
                print '    TRAPEZOIDAL'
        
        
            
        