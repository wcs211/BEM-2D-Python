# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import os
import time
import socket

class printOutput(object):
    'Toolkit for formated terminal output'
    def progTitle(self, version):
        print "/*---------------------------------------------------------------------------*\\"
        print "| .   )\\      Py thon       |                                                 |"
        print "| \\`.-' `-oo                | Written by:  Lehigh Biofluids Group             |"
        print "|  ) _  __,0) B  oundary    | Version:     %s                              |" % version
        print "| /.' )/      E  lement     | Web:         http://www.lehigh.edu/             |"
        print "| '           M  ethod      |                                                 |"
        print "\\*---------------------------------------------------------------------------*/"
        print 'Date     : %s' % time.strftime("%B %d, %Y")
        print 'Time     : %s' % time.strftime("%I:%M:%S %p")
        print 'Host     : %s' % socket.gethostname()
        print 'PID      : %i' % os.getpid()
        print 'Case     : %s' % os.getcwd()
        print ''
        print '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
        
    def calcInput(self,alpha_max,Re,AoA_max,delT):
        print 'Calculated inputs:'
        print '     alpha_max = %f' % alpha_max
        print '     Re = %e'        % Re
        print '     AoA_max = %f'   % AoA_max
        print '     delT = %f'      % delT
        
    def initializeOutput(self,t):
        print '\nInitializing the flow solution for FLOW TIME = %f\n' % t
        
    def timestepHeader(self,i_t,t):
        print '==============================================================================='
        print ' TIME-STEP NUMBER = %i, FLOW TIME = %f' % (i_t-1,t)
        print '-------------------------------------------------------------------------------'
        
    def fsiHeader(self):
        print '|  Iter.  |  Relax.    |  Max Displ |  Max Res  | RMS Res Norm | Max Res Norm |'
        print '+---------+------------+------------+-----------+--------------+--------------+'
        
    def fsiIterOut(self,outerCorr,fsiRelaxationFactor,maxDispl,maxMagFsiResidual,fsiResidualNorm,maxFsiResidualNorm):
        print '| %7i |   %.2E |   %.2E |  %.2E |     %.2E |     %.2E |' % (
        outerCorr,fsiRelaxationFactor,maxDispl,maxMagFsiResidual,
        fsiResidualNorm,maxFsiResidualNorm )
        
    def fsiConverged(self):
        print '| SOLUTION CONVERGED!                                                         |'
        print '+-----------------------------------------------------------------------------+'
        
    def fsiNotConverged(self):
        print '| WARNING! MAX INNER-LOOP ITERATIONS REACHED                                  |'
        print '+-----------------------------------------------------------------------------+'
        
    def solutionOutput(self,D_visc,Cf,Cl,Ct,Cpow,Gamma):
        print '| Solution Information:                                                       |'
        print '|     D_visc   = %13e                                                |' % D_visc
        print '|     Cf       = %13e                                                |' % Cf    
        print '|     Cl       = %13e                                                |' % Cl    
        print '|     Ct       = %13e                                                |' % Ct    
        print '|     Cpow     = %13e                                                |' % Cpow  
        print '|     Gamma    = %13e                                                |' % Gamma 
        
    def solutionAvgOutput(self,Cl_avg,Ct_avg,Tnet_avg,D_avg,Pow_avg,Cpow_avg):
        print '|                                                                             |'
        print '| Solution Average Cycle Information:                                         |'
        print '|     Cl_avg   = %13e                                                |' % Cl_avg  
        print '|     Ct_avg   = %13e                                                |' % Ct_avg  
        print '|     Tnet_avg = %13e                                                |' % Tnet_avg
        print '|     D_avg    = %13e                                                |' % D_avg   
        print '|     Pow_avg  = %13e                                                |' % Pow_avg 
        print '|     Cpow_avg = %13e                                                |' % Cpow_avg
        
    def solutionCompleteOutput(self, perDone):
       print '|     %3i%% Complete                                                           |' % perDone
       print '+-----------------------------------------------------------------------------+\n'