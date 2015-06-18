import numpy as np

# Constants that determine other parameters and don't actually need lookup
RF = 2*np.pi # Reduced frequency
MU = 0.001003

P = PARAMETERS = {
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Geometry Definition                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 'SW_GEOMETRY':    'FP'
, 'N_BODY':         100
, 'C':              0.2
, 'K':              2.-(12.4/180)
, 'EPSILON':        0.075
, 'T_MAX':          0.008
, 'CE':             0.4
, 'S':              0.01

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-step and Misc. Parameters                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'COUNTER':        201
, 'DEL_T':          np.pi*0.1/RF
, 'DSTEP':          10**-5
, 'TSTEP':          10**-5
, 'VERBOSITY':      1      

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'V0':             -1.0
, 'THETA_MAX':      5*np.pi/180
, 'F':              RF/(2*np.pi)
, 'PHI':            0
, 'RHO':            998.2
, 'SW_KUTTA':       1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Skin Friction Solver Constants                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Structural Solver Constants                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'INT_METHOD':         'HHT'
, 'M_TYPE':             'consistent'
, 'ALPHA':              0.2
, 'BETA':               0.25*(1+0.2)**2
, 'GAMMA':              0.5+0.2
, 'N_ELEMENTS_S':       40
, 'MATERIAL':           'Polyethylene'
, 'E':                  3.8e9
, 'RHO_S':              935
, 'FRAC_DELT':          0.1
, 'FLEX_RATIO':         0.3
, 'T_CONST':            0.95

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FSI Coupling Constants                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_FSI':             1
, 'N_OUTERCORR_MAX':    1000
, 'OUTER_CORR_TOL':     1e-5
, 'FIXED_PT_RELAX':     0.00001
, 'COUPLING_SCHEME':    'Aitken'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Solver Switches                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_FREE_SWIM':   0
, 'SW_VISC_DRAG':   0
, 'SW_INTERP_MTD':  1
, 'SW_CNST_THK_BM': 1
, 'SW_RAMP':        1
}
##### The Following parameters are based on perviously declared variables #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Body Motion Parameters                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['T'] = [P['DEL_T'] * i for i in xrange(P['COUNTER'])]
P['THETA'] = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
P['THETA_MINUS'] = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
P['THETA_PLUS'] = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
P['HEAVE'] = [0 for i in xrange(P['COUNTER'])]
P['HEAVE_MINUS'] = [0 for i in xrange(P['COUNTER'])]
P['HEAVE_PLUS'] = [0 for i in xrange(P['COUNTER'])]

# Constants dependent on declared parameters
P['DELTA_CORE'] = (0.005*P['THETA_MAX']+0.09)*P['C']
P['RE'] = P['RHO']*-P['V0']*P['C']/MU
