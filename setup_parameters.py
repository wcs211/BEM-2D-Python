# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:36:05 2015

@author: biofluids_1
"""

import numpy as np

# Constants that determine other parameters and don't actually need lookup
RF = np.pi # Reduced frequency
MU = 0.001003

P = PARAMETERS = {

  'COUNTER':        401
, 'DEL_T':          np.pi*0.01/RF
#, 'DEL_T':          0.01
, 'DSTEP':          10**-5
, 'TSTEP':          10**-5

, 'N_BODY':         100
, 'C':              1.
, 'K':              2.-(12.4/180)
, 'EPSILON':        0.075
, 'V0':             -1.0
, 'THETA_MAX':      5.73*np.pi/180
, 'H_C':            0
, 'F':              RF/(2*np.pi)
#, 'F':              0.5
, 'PHI':            0

, 'CE':             0.4
, 'S':              0.1
, 'RHO':            998.2

, 'SW_GEOMETRY':    'VDV'
, 'SW_KUTTA':       1

}

# Constants dependent on declared parameters
P['DELTA_CORE'] = (0.005*P['THETA_MAX']+0.09)*P['C']
P['RE'] = P['RHO']*-P['V0']*P['C']/MU