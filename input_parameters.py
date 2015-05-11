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
, 'VERBOSITY':      1


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_BODY':         100
, 'C':              1.
, 'K':              2.-(12.4/180)
, 'EPSILON':        0.075
, 'V0':             -1.0
, 'THETA_MAX':      5.73*np.pi/180
, 'F':              RF/(2*np.pi)
#, 'F':              0.5
, 'PHI':            0
, 'T_MAX':          0.01

, 'CE':             0.4
, 'S':              0.1
, 'RHO':            998.2

, 'SW_GEOMETRY':    'VDV'
, 'SW_KUTTA':       1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Skin Friction Solver Constants                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Structural Solver Constants                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'INT_METHOD':         'HHT'
, 'M_TYPE':             'consistent'
, 'ALPHA':              0.1
, 'BETA':               0.25*(1+0.1)**2
, 'GAMMA':              0.5+0.1
, 'N_ELEMENTS_S':       100
, 'MATERIAL':           'Aluminum'
, 'E':                  75.0e9
, 'RHO_S':              2710
, 'FRAC_DELT':          0.1
, 'FLEX_RATIO':         0.3
, 'T_CONST':            0.95


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FSI Coupling Constants                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_OUTERCORR_MAX':    50
, 'OUTER_CORR_TOL':     1e-5
, 'FIXED_PT_RELAX':     0.001
, 'COUPLING_SCHEME':    'Aitken'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Solver Switches                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SWITCH_KUTTA':       1
, 'SWITCH_FREE_SWIM':   0
, 'SWITCH_VISC_DRAG':   0
, 'SWITCH_INTERP_MTD':  1
, 'SWITCH_CNST_THK_BM': 0



}

# Constants dependent on declared parameters
P['DELTA_CORE'] = (0.005*P['THETA_MAX']+0.09)*P['C']
P['RE'] = P['RHO']*-P['V0']*P['C']/MU