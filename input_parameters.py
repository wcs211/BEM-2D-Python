import numpy as np

# Constants that determine other parameters and don't actually need lookup
MU = 0.001003

P = PARAMETERS = {
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Data I/O                                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  'SW_SAVE_DATA':       False
, 'SW_SV_L_CYCLE':      True
, 'SAVE_EVERY':         1
, 'OUTPUT_DIR':         '/home/wcs211/BEM-2D-Python/data'
, 'START_FROM':         'zeroTime'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Geometry Definition                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_GEOMETRY':        'FP'
, 'N_BODY':             100
, 'C':                  0.100
, 'K':                  2.-(12.4/180)
, 'EPSILON':            0.075
, 'T_MAX':              0.0011
, 'CE':                 0.4
, 'S':                  0.01

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-step and Misc. Parameters                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_STEP':             25
, 'N_CYC':              10
, 'DSTEP':              10**-5
, 'TSTEP':              10**-5
, 'VERBOSITY':          1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'V0':                 -0.05
, 'THETA_MAX':          np.arctan(0.025132741/0.05)
, 'HEAVE_MAX':          0.0025
, 'F':                  0.4
, 'PHI':                0.5*np.pi
, 'RHO':                998.2
, 'SW_KUTTA':           False

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Skin Friction Solver Constants                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Structural Solver Constants                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'INT_METHOD':         'HHT'
, 'M_TYPE':             'consistent'
, 'ALPHA':              0.02
, 'BETA':               0.25*(1+0.02)**2
, 'GAMMA':              0.5+0.02
, 'N_ELEMENTS_S':       100
, 'MATERIAL':           'Polyethylene'
, 'E':                  3.8e9
, 'RHO_S':              935
, 'FRAC_DELT':          1.0
, 'FLEX_RATIO':         0.056410256
, 'T_CONST':            0.95

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FSI Coupling Constants                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_FSI':             True
, 'N_OUTERCORR_MAX':    1000
, 'OUTER_CORR_TOL':     1e-7
, 'FIXED_PT_RELAX':     1e-5
, 'COUPLING_SCHEME':    'Aitken'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Solver Switches                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_ROLLUP':          True
, 'SW_FREE_SWIM':       False
, 'SW_VISC_DRAG':       False
, 'SW_INTERP_MTD':      True
, 'SW_CNST_THK_BM':     True
, 'SW_RAMP':            True
, 'SW_PLOT_FIG':        True
}


##### The Following parameters are based on perviously declared variables #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Multi-Swimmer Parameters                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['N_SWIMMERS']  = 1
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
slope  = P['F']
offset = 0.25 * P['N_CYC'] / P['F']
RAMP             = [0.5*np.tanh(slope *(P['T'][i] - offset))+0.5 for i in xrange(P['COUNTER'])]
RAMP_P           = [0.5*np.tanh(slope *((P['T'][i] + P['TSTEP'])-offset))+0.5 for i in xrange(P['COUNTER'])]
RAMP_M           = [0.5*np.tanh(slope *((P['T'][i] - P['TSTEP'])-offset))+0.5 for i in xrange(P['COUNTER'])]
P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI'])  * RAMP[i] for i in xrange(P['COUNTER'])]
P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) * RAMP_M[i] for i in xrange(P['COUNTER'])]
P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) * RAMP_P[i] for i in xrange(P['COUNTER'])]
H_DOT            = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * (P['T'][i] - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * P['T'][i] + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * (P['T'][i] - offset)))**2 * np.sin(2. * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
H_DOT_PLUS       = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * ((P['T'][i] + P['TSTEP']) - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * ((P['T'][i] + P['TSTEP']) - offset)))**2 * np.sin(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
H_DOT_MINUS      = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * ((P['T'][i] - P['TSTEP']) - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * ((P['T'][i] - P['TSTEP']) - offset)))**2 * np.sin(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
P['THETA']       = [np.arctan(H_DOT[i] / P['V0']) for i in xrange(P['COUNTER'])]
P['THETA_MINUS'] = [np.arctan(H_DOT_MINUS[i] / P['V0']) for i in xrange(P['COUNTER'])]
P['THETA_PLUS']  = [np.arctan(H_DOT_PLUS[i] / P['V0']) for i in xrange(P['COUNTER'])]

#P['THETA']       = [np.tanh(P['T'][i])*5./180.*np.pi for i in xrange(P['COUNTER'])]
#P['THETA_MINUS'] = [np.tanh(P['T'][i])*5./180.*np.pi for i in xrange(P['COUNTER'])]
#P['THETA_PLUS']  = [np.tanh(P['T'][i])*5./180.*np.pi for i in xrange(P['COUNTER'])]
#P['HEAVE']       = [0. for i in xrange(P['COUNTER'])]
#P['HEAVE_MINUS'] = [0. for i in xrange(P['COUNTER'])]
#P['HEAVE_PLUS']  = [0. for i in xrange(P['COUNTER'])]
#
#P['THETA']       = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
#P['THETA_MINUS'] = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
#P['THETA_PLUS']  = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
#P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i]) * np.tanh(3 * P['T'][i])  * RAMP[i] for i in xrange(P['COUNTER'])]
#P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP'])) * np.tanh(3 * (P['T'][i] - P['TSTEP'])) * RAMP[i] for i in xrange(P['COUNTER'])]
#P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP'])) * np.tanh(3 * (P['T'][i] + P['TSTEP'])) * RAMP[i] for i in xrange(P['COUNTER'])]

# Constants dependent on declared parameters
P['DELTA_CORE']  = (0.005*P['THETA_MAX']+0.09)*P['C']
P['RE']          = P['RHO']*-P['V0']*P['C']/MU
