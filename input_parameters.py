import numpy as np
from functions_general import intermittent, multi_kinematics, accel_multi_kinematics

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
, 'N_BODY':             250
, 'C':                  0.1
, 'B':                  1.0
, 'K':                  2.-(12.4/180)
, 'EPSILON':            0.075
, 'T_MAX':              0.01
, 'CE':                 0.40
, 'S':                  0.15

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-step and Misc. Parameters                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_STEP':             250
, 'N_CYC':              10
, 'DSTEP':              1e-5
, 'TSTEP':              1e-5
, 'VERBOSITY':          1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'V0':                 -0.05
, 'THETA_MAX':          0.00
, 'HEAVE_MAX':          0.0125
, 'F':                  0.3
, 'DC':                 0.5
, 'SIG_WEIGHT':         [0., 1., 0., 0.] # [sine, square, triangle, sawtooth]
, 'PHI':                np.pi/2.
, 'RHO':                1000.
, 'NU':                 1.003e-6
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
, 'N_ELEMENTS_S':       1
, 'MATERIAL':           'Synthesized'
, 'E':                  2.0e9
, 'RHO_S':              1000.
, 'FRAC_DELT':          1.0
, 'FLEX_RATIO':         1e-5
, 'T_CONST':            0.00125
, 'KAPPA':              0.1
, 'ZETA':               0.0

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FSI Coupling Constants                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_FSI':             True
, 'SW_SPRING':          True
, 'N_OUTERCORR_MAX':    200
, 'OUTER_CORR_TOL':     1e-7
, 'FIXED_PT_RELAX':     1e-5
, 'COUPLING_SCHEME':    'Aitken'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Plotting Options                                                            #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'X_TICKS':            [-15., 15., 0.2]
, 'X_FIELD':            [-0.05, 1.75]
, 'Z_FIELD':            [-0.5, 0.5]
, 'X_BODY':             [-0.04, 0.14]
, 'Z_BODY':             [-0.05, 0.05]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Solver Switches                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_MULTI':           True
, 'SW_ROLLUP':          True
, 'SW_FLAT_WAKE':       False
, 'SW_FREE_SWIM':       False
, 'SW_INTERMITTENT':    False
, 'SW_VISC_DRAG':       False
, 'SW_INTERP_MTD':      True
, 'SW_CNST_THK_BM':     True
, 'SW_RAMP':            False
, 'SW_PLOT_FIG':        True
}


##### The Following parameters are based on perviously declared variables #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Virtual Body Properties                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['SW_ADDED_DRAG'] =    False
P['DRAG_LAW']      =    'FORM'
P['S_WP']          =    10.
P['M_STAR']        =    0.5
P['CD_BOD']        =    0.010                                      # Virtual body drag coefficient.
P['BETA_M']        =    1.0                                        # Added mass coefficient.
P['S_P']           =    P['C'] * P['B']                            # Planform area of the propulsor
P['M']             =    P['RHO'] * P['S_P'] * P['C'] / P['M_STAR'] # Virtual body mass, kg.
P['S_WB']          =    (P['S_WP'] - 2.) * P['S_P']                # Surface area of virtual body, m^2.
P['S_W']           =    P['S_WB'] + 2. * P['S_P']                  # Surface area of virtual body, m^2.
P['AR_B']          =    3.5                                        # Ellipsoidal body aspect ratio.
P['L_T']           =    P['C'] + np.sqrt(P['S_WB'] * P['AR_B'])    # Virtual body length, meters.


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
# Zero Angle of Attack Motions
if P['SW_RAMP']:
    slope  = 2. * P['F']
    offset = P['N_CYC'] / 4. / P['F'] - 0.*slope
    P['RAMP']        = [0.5*np.tanh(slope *(P['T'][i] - offset))+0.5 for i in xrange(P['COUNTER'])]
    P['RAMP_P']      = [0.5*np.tanh(slope *((P['T'][i] + P['TSTEP'])-offset))+0.5 for i in xrange(P['COUNTER'])]
    P['RAMP_M']      = [0.5*np.tanh(slope *((P['T'][i] - P['TSTEP'])-offset))+0.5 for i in xrange(P['COUNTER'])]
else:
    P['RAMP']        = [1.0 for i in xrange(P['COUNTER'])]
    P['RAMP_P']      = [1.0 for i in xrange(P['COUNTER'])]
    P['RAMP_M']      = [1.0 for i in xrange(P['COUNTER'])]
    
if P['SW_INTERMITTENT']:
    (THETA, HEAVE, period) = intermittent(P['HEAVE_MAX'], P['THETA_MAX'], P['PHI'], P['DC'], P['F'], P['N_STEP'], P['N_CYC'], 0)
    P['THETA'] =THETA.tolist()
    P['HEAVE'] = HEAVE.tolist()
    
    (THETA, HEAVE, period) = intermittent(P['HEAVE_MAX'], P['THETA_MAX'], P['PHI'], P['DC'], P['F'], P['N_STEP'], P['N_CYC'], -1)
    P['THETA_MINUS'] =THETA.tolist()
    P['HEAVE_MINUS'] = HEAVE.tolist()
    
    (THETA, HEAVE, period) = intermittent(P['HEAVE_MAX'], P['THETA_MAX'], P['PHI'], P['DC'], P['F'], P['N_STEP'], P['N_CYC'], 1)
    P['THETA_PLUS'] =THETA.tolist()
    P['HEAVE_PLUS'] = HEAVE.tolist()
    
elif P['SW_MULTI']:
    (sigTheta, sigThetaMinus, sigThetaPlus) = multi_kinematics(P, 0., scale=[1.0, 1.143727574, 1.693954952, 2.690184034], rate=5)
    P['THETA']       = [P['THETA_MAX'] * sigTheta[i]      * P['RAMP'][i]   for i in xrange(P['COUNTER'])]
    P['THETA_MINUS'] = [P['THETA_MAX'] * sigThetaMinus[i] * P['RAMP_M'][i] for i in xrange(P['COUNTER'])]
    P['THETA_PLUS']  = [P['THETA_MAX'] * sigThetaPlus[i]  * P['RAMP_P'][i] for i in xrange(P['COUNTER'])]
    
    (sigHeave, sigHeaveMinus, sigHeavePlus) = multi_kinematics(P, P['PHI'], scale=[1.0, 1.143727574, 1.693954952, 2.690184034], rate=5)
    P['HEAVE']       = [P['HEAVE_MAX'] * sigHeave[i]      * P['RAMP'][i]   for i in xrange(P['COUNTER'])]
    P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * sigHeaveMinus[i] * P['RAMP_M'][i] for i in xrange(P['COUNTER'])]
    P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * sigHeavePlus[i]  * P['RAMP_P'][i] for i in xrange(P['COUNTER'])]
    
    inertia = accel_multi_kinematics(P, P['PHI'], scale=[1.0, 1.143727574, 1.693954952, 2.690184034], rate=5)[0]
    P['INERTIA']     = [P['HEAVE_MAX'] * inertia[i] * P['RAMP'][i] for i in xrange(P['COUNTER'])]
    
else:
    # Zero attack angle motions
    P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI'])  * P['RAMP'][i] for i in xrange(P['COUNTER'])]
    P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) * P['RAMP_M'][i] for i in xrange(P['COUNTER'])]
    P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) * P['RAMP_P'][i] for i in xrange(P['COUNTER'])]
    H_DOT            = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * (P['T'][i] - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * P['T'][i] + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * (P['T'][i] - offset)))**2 * np.sin(2. * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
    H_DOT_PLUS       = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * ((P['T'][i] + P['TSTEP']) - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * ((P['T'][i] + P['TSTEP']) - offset)))**2 * np.sin(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
    H_DOT_MINUS      = [2. * np.pi * P['F'] * P['HEAVE_MAX'] * (0.5 * np.tanh(slope * ((P['T'][i] - P['TSTEP']) - offset)) + 0.5) * np.cos(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) + P['HEAVE_MAX'] / (np.cosh(slope * ((P['T'][i] - P['TSTEP']) - offset)))**2 * np.sin(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
    P['THETA']       = [np.arctan(H_DOT[i] / P['V0']) for i in xrange(P['COUNTER'])]
    P['THETA_MINUS'] = [np.arctan(H_DOT_MINUS[i] / P['V0']) for i in xrange(P['COUNTER'])]
    P['THETA_PLUS']  = [np.arctan(H_DOT_PLUS[i] / P['V0']) for i in xrange(P['COUNTER'])]

# Constants dependent on declared parameters
#P['DELTA_CORE']  = (0.005*P['THETA_MAX']+0.09)*P['C']
P['DELTA_CORE']  = 0.05 * P['C']
#P['DELTA_CORE']  = 2.5e-5 * P['C']
P['RE']          = P['RHO']*-P['V0']*P['C']/MU
