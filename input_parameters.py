import numpy as np
from functions_general import ramp

# Constants that determine other parameters and don't actually need lookup
MU = 0.001003

P = PARAMETERS = {
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Data I/O                                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  'SW_SAVE_DATA':       False
, 'SW_SV_L_CYCLE':      False
, 'SW_SV_FORCES':       True
, 'SAVE_EVERY':         625
, 'OUTPUT_DIR':         '/home/wcs211/BEM-2D-Python/data'
, 'START_FROM':         'zeroTime'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Geometry Definition                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_GEOMETRY':        'TD'
, 'N_BODY':             2000
, 'C':                  0.10
, 'K':                  2.-(12.4/180)
, 'EPSILON':            0.075
, 'T_MAX':              0.00025
, 'CE':                 0.40
, 'S':                  0.15

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-step and Misc. Parameters                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_STEP':             250
, 'N_CYC':              3
, 'DSTEP':              1e-5
, 'TSTEP':              1e-5
, 'VERBOSITY':          1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'V0':                 -0.05
, 'THETA_MAX':          0.00
, 'HEAVE_MAX':          0.005
, 'F':                  0.7
#, 'PHI':                -3.0*np.pi/8.0
, 'PHI':                -1.0*np.pi/2.
, 'RHO':                1000.
, 'SW_KUTTA':           False

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Skin Friction Solver Constants                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Structural Solver Constants                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'INT_METHOD':         'HHT'
, 'M_TYPE':             'consistent'
, 'ALPHA':              0.03
, 'BETA':               0.25*(1+0.03)**2
, 'GAMMA':              0.5+0.03
, 'N_ELEMENTS_S':       32
, 'MATERIAL':           'Synthesized'
, 'E':                  2.0e9
, 'RHO_S':              1300
, 'FRAC_DELT':          1.0
, 'FLEX_RATIO':         1e-7
, 'T_CONST':            0.0025

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FSI Coupling Constants                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_FSI':             True
, 'N_OUTERCORR_MAX':    200
, 'OUTER_CORR_TOL':     1e-7
, 'FIXED_PT_RELAX':     1e-5
, 'COUPLING_SCHEME':    'Aitken'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Solver Switches                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_ROLLUP':          False
, 'SW_FREE_SWIM':       False
, 'SW_VISC_DRAG':       False
, 'SW_INTERP_MTD':      True
, 'SW_CNST_THK_BM':     True
, 'SW_RAMP':            True
, 'SW_PLOT_FIG':        False
}


##### The Following parameters are based on perviously declared variables #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Multi-Swimmer Parameters                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['N_SWIMMERS']  = 1
#P['X_START'] = [0.25, 0.25, 0.25, 0.0, 0.0]
#P['Z_START'] = [ -0.25, 0.00, 0.25, -0.125, 0.125]

P['X_START']     = [i * 0.0 for i in xrange(P['N_SWIMMERS'])]
P['Z_START']     = [i * 0.4 for i in xrange(P['N_SWIMMERS'])]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-stepParameters                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['DEL_T']   = 1. / P['F'] / P['N_STEP']
P['COUNTER'] = P['N_CYC'] * P['N_STEP'] + 1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Body Motion Parameters                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['T']           = [P['DEL_T'] * i for i in xrange(P['COUNTER'])]
# Zero Angle of Attack Motions
if P['SW_RAMP']:
    slope  = 2. * P['F']
    offset = P['N_CYC'] / 16. / P['F'] - 0.*slope
    RAMP             = [0.5*np.tanh(slope *(P['T'][i] - offset))+0.5 for i in xrange(P['COUNTER'])]
    RAMP_P           = [0.5*np.tanh(slope *((P['T'][i] + P['TSTEP'])-offset))+0.5 for i in xrange(P['COUNTER'])]
    RAMP_M           = [0.5*np.tanh(slope *((P['T'][i] - P['TSTEP'])-offset))+0.5 for i in xrange(P['COUNTER'])]
else:
    RAMP             = [1.0 for i in xrange(P['COUNTER'])]
    RAMP_P           = [1.0 for i in xrange(P['COUNTER'])]
    RAMP_M           = [1.0 for i in xrange(P['COUNTER'])]
#    RAMP             = ramp(np.asarray(P['T']),0.1,0.) + ramp(np.asarray(P['T']),-0.1,10.)
#    RAMP_P           = ramp(np.asarray(P['T']) + P['TSTEP'],0.1,0.) + ramp(np.asarray(P['T']) + P['TSTEP'],-0.1,10.)
#    RAMP_M           = ramp(np.asarray(P['T']) - P['TSTEP'],0.1,0.) + ramp(np.asarray(P['T']) - P['TSTEP'],-0.1,10.)
    
#P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI'])  * RAMP[i] for i in xrange(P['COUNTER'])]
#P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) * RAMP_M[i] for i in xrange(P['COUNTER'])]
#P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) * RAMP_P[i] for i in xrange(P['COUNTER'])]
#H_DOT            = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * (P['T'][i] - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * P['T'][i] + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * (P['T'][i] - offset)))**2 * np.sin(2. * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
#H_DOT_PLUS       = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * ((P['T'][i] + P['TSTEP']) - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * ((P['T'][i] + P['TSTEP']) - offset)))**2 * np.sin(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
#H_DOT_MINUS      = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * ((P['T'][i] - P['TSTEP']) - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * ((P['T'][i] - P['TSTEP']) - offset)))**2 * np.sin(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
#P['THETA']       = [np.arctan(H_DOT[i] / P['V0']) for i in xrange(P['COUNTER'])]
#P['THETA_MINUS'] = [np.arctan(H_DOT_MINUS[i] / P['V0']) for i in xrange(P['COUNTER'])]
#P['THETA_PLUS']  = [np.arctan(H_DOT_PLUS[i] / P['V0']) for i in xrange(P['COUNTER'])]

# Heaving and Pitching Motions
P['THETA']       = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) * RAMP[i] for i in xrange(P['COUNTER'])]
P['THETA_MINUS'] = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) * RAMP_M[i] for i in xrange(P['COUNTER'])]
P['THETA_PLUS']  = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) * RAMP_P[i] for i in xrange(P['COUNTER'])]
P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i]) * RAMP[i] for i in xrange(P['COUNTER'])]
P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP'])) * RAMP_M[i] for i in xrange(P['COUNTER'])]
P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP'])) * RAMP_P[i] for i in xrange(P['COUNTER'])]

P['THETA_DOT']     = [2 * np.pi * P['F'] * P['THETA_MAX'] * np.cos(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) * RAMP[i] for i in xrange(P['COUNTER'])]
P['THETA_DOT_DOT'] = [-4 * np.pi**2 * P['F']**2 * P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) * RAMP[i] for i in xrange(P['COUNTER'])]

# Constants dependent on declared parameters
#P['DELTA_CORE']  = (0.005*P['THETA_MAX']+0.09)*P['C']
P['DELTA_CORE']  = 0.05 * P['C']
#P['DELTA_CORE']  = 2.5e-5 * P['C']
P['RE']          = P['RHO']*-P['V0']*P['C']/MU
