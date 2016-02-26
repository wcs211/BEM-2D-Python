import os
import sys
import numpy as np
from scipy import array
from scipy.interpolate import PchipInterpolator

def panel_vectors(x, z):
    """
    Calculates the normal and tangential unit vectors for 2D body panel.
    
    Args:
        x (float: 1 x N): x panel endpoint coordinates
        z (float: 1 x N): z panel endpoint coordinates
        
    Returns:
        tx (float: 1 x N-1): x unit vector component of the panel's tangent vector
        tz (float: 1 x N-1): z unit vector component of the panel's tangent vector
        nx (float: 1 x N-1): x unit vector component of the panel's normal vector
        nz (float: 1 x N-1): z unit vector component of the panel's normal vector
        lpanel (float: 1 x N-1): length of each panel
    """
    lpanel = np.sqrt((x[1:]-x[:-1])**2 + (z[1:]-z[:-1])**2)
    tx = (x[1:]-x[:-1])/lpanel
    tz = (z[1:]-z[:-1])/lpanel
    nx = -tz
    nz = tx
    return (tx,tz,nx,nz,lpanel)

def point_vectors(xdp, xdm, zdp, zdm):
    """
    Calculates the normal and tangential unit vectors for each midpoint or 
    collocation point.
    
    Args:
        xdp (float): 
        xdm (float):
        zdp (float):
        zdm (float):
        
    Returns:
        tx (float):
        tz (float):
        nx (float):
        nz (float):
    """
    tx = (xdp-xdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    tz = (zdp-zdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    nx = -tz
    nz = tx
    return(tx,tz,nx,nz)

def archive(array, axis=0):
    """
    Shifts array values along an axis (row-wise by default).

    Used for arrays that keep past values for differencing with respect to time.

    Args:
        array: The array that will be manipulated.
        axis: The axis to shift values along (0==row-wise, 1==column-wise).
    """
    if len(np.shape(array)) == 1:
        array[1:] = array[:-1]
    elif len(np.shape(array)) == 2:
        if axis == 0:
            array[1:,:] = array[:-1,:]
        else:
            array[:,1:] = array[:,:-1]

def transformation(xt,zt,xi,zi):
    """
    Transforms points into a panel's reference frame.
    
    Velocity and velocity potential equations are defined in panel coordinates 
    so a transformation should be done. Each row of xp1/xp2/zp is a target, and 
    each column is an influence.
    
    Args:
        xt (float): x coordinate of target points
        zt (float): z coordinate of target points
        xi (float): x coordinate of influences
        zi (float): z coordinate of influences
        
    Returns:
        xp1 (float):
        xp2 (float):
        zp (float):
    """
    # NI is N influences, NT is N targets
    NT = np.size(xt)
    NI = np.size(xi)-1

    (tx,tz,nx,nz) = panel_vectors(xi,zi)[:-1]

    # Intermediary variables to reduce number of tile/repeat operations
    # From normalvectors: tx==nz, tz==-nx
    x_tile = np.repeat(xt[:,np.newaxis],NI,1) - np.repeat(xi[:-1,np.newaxis].T,NT,0)
    z_tile = np.repeat(zt[:,np.newaxis],NI,1) - np.repeat(zi[:-1,np.newaxis].T,NT,0)
    tx_tile = np.repeat(tx[:,np.newaxis].T,NT,0)
    tz_tile = np.repeat(tz[:,np.newaxis].T,NT,0)

    # Transforming left side collocation points from global to local coordinates
    xp1 = x_tile*tx_tile + z_tile*tz_tile
    zp = x_tile*(-tz_tile) + z_tile*tx_tile

    # Transforming right side panel points into local coordinate system
    dummy = (xi[1:]-xi[:-1])*tx + (zi[1:]-zi[:-1])*tz
    xp2 = xp1 - np.repeat(dummy[:,np.newaxis].T,NT,0)

    return(xp1,xp2,zp)

def absoluteToBody(Body, Solid, P, i):
    """
    Transforms absolute reference frame to body reference frame. Needed in 
    FSI simulations after solid displacements are determined.
    
    Args:
        Body (object):
        Solid (object):
        P (dict):
        i (int):
    """
    THETA = P['THETA'][i]
    
    Body.BF.x = ((Body.AF.x - Body.AF.x_le) * np.cos(-1*THETA) - (Body.AF.z - Body.AF.z_le) * np.sin(-1*THETA))
    Body.BF.z = ((Body.AF.z - Body.AF.z_le) * np.cos(-1*THETA) + (Body.AF.x - Body.AF.x_le) * np.sin(-1*THETA))
    Body.BF.x_col = ((Body.BF.x[1:] + Body.BF.x[:-1])/2)
    Body.BF.z_col = ((Body.BF.z[1:] + Body.BF.z[:-1])/2)

    Solid.nodesNew[:,0] = (Solid.nodes[:,0] - Body.AF.x_le) * np.cos(-1*THETA) - (Solid.nodes[:,1] - Body.AF.z_le) * np.sin(-1*THETA)
    Solid.nodesNew[:,1] = (Solid.nodes[:,1] - Body.AF.z_le) * np.cos(-1*THETA) + (Solid.nodes[:,0] - Body.AF.x_le) * np.sin(-1*THETA)

def ramp(t, slope, startTime):
    """
    This function can generate a ramp signal based on the following inputs:

    Args:
        t: array of time samples
        slope: slope of the ramp signal
        startTime: location where the ramp turns on
    """
    # Get the number of samples in the output signal
    N = t.size

    # Initialize the ramp signal
    r = np.zeros(N)

    # Find the index where the ramp turns on
    if (np.median(np.diff(t)) > 0):
        startInd = np.min((t>=startTime).nonzero())
        popInd =np.arange(startInd,N)
    elif (np.median(np.diff(t)) < 0):
        # Time-reversed ramp
        startTime = -1. * startTime
        startInd = np.max((t>=startTime).nonzero())
        popInd = np.arange(startInd)
        slope = -1. * slope

    # For indicies greater than the start time, compute the
    # proper signal value based on slope
    r[popInd] = slope * (t[popInd] + startTime) - 2 * startTime * slope

    return (r)

def geom_setup(P, PC, Swimmer, solid=None, FSI=None, PyFEA=None):
    """
    Creates new objects need to start a new simulation.
    
    Args:
        P (dict):
        PC (object):
        solid (object, optional):
        FSI (object, optional):
        PyFEA (object, optional):
    
    Returns:
        SwiL (list):
        GeoL (list):
        MotL (list):
        Swimmers (list):
        SolidL (list):
        FSIL (list):
        PyFEAL (list):
    """
    # Initialize lists of objects 
    SwiL     = [None for x in xrange(P['N_SWIMMERS'])]
    GeoL     = [None for x in xrange(P['N_SWIMMERS'])]
    MotL     = [None for x in xrange(P['N_SWIMMERS'])]
    Swimmers = [None for x in xrange(P['N_SWIMMERS'])]
    SolidL   = [None for x in xrange(P['N_SWIMMERS'])]
    FSIL     = [None for x in xrange(P['N_SWIMMERS'])]
    PyFEAL   = [None for x in xrange(P['N_SWIMMERS'])]

    for i in xrange(P['N_SWIMMERS']):
        # Add the Swimmer's parameter class to the list
        SwiL[i] = PC.SwimmerParameters(P['CE'], P['DELTA_CORE'], P['SW_KUTTA'])
        
        # Determine which geometry parameter list to create and add it to the 
        # list. If the geometry parameter does not exist, raise a value error.
        if (P['SW_GEOMETRY'] == 'FP'):
            GeoL[i] = PC.GeoFPParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
        elif (P['SW_GEOMETRY'] == 'TD'):
            GeoL[i] = PC.GeoTDParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
        elif (P['SW_GEOMETRY'] == 'VDV'):
            GeoL[i] = PC.GeoVDVParameters(P['N_BODY'], P['S'], P['C'], P['K'], P['EPSILON'])
        else:
            print "ERROR! Invalid geometry type. Valid geometry types are:'"
            print "    'FP'"
            print "    'TD'"
            print "    'VDV'"
            raise ValueError('ERROR! Invalid geometry type.')

        # Add the Swimmer's motion parameters to the motion parameter list.
        MotL[i] = PC.MotionParameters(P['X_START'][i], P['Z_START'][i], P['V0'], P['THETA_MAX'], P['F'], P['PHI'])

        # Create a Swimmer object and add it to the Swimmer object list.
        Swimmers[i] = Swimmer(SwiL[i], GeoL[i], MotL[i], P['COUNTER']-1)

        # Create more objects if this is an FSI simulation.
        if P['SW_FSI']:
            # Create a solid mesh object and add it to the solid mesh list.
            SolidL[i] = solid(Swimmers[i].Body, P['N_ELEMENTS_S'], P['T_MAX'])
            
            # Create an FSI object and add it to the FSI list.
            FSIL[i] = FSI(Swimmers[i].Body, SolidL[i])
            
            # Create an FEA object and add it to the list.
            PyFEAL[i] = PyFEA(SolidL[i], P['SW_SPRING'], P['FRAC_DELT'], P['DEL_T'], P['E'], P['RHO_S'])

            # Initialize the mesh for the solid object.
            SolidL[i].initMesh()
            
            # Create the appropriate solid properties base on the input 
            # geometry. If it does not exist, raise a value error and inform 
            # the user of valid choices.
            if (P['SW_GEOMETRY'] == 'FP'):
                SolidL[i].initThinPlate(P['T_MAX'],P['C'],P['SW_CNST_THK_BM'],P['T_CONST'],P['FLEX_RATIO'])
            elif (P['SW_GEOMETRY'] == 'TD'):
                SolidL[i].initTearDrop(P['T_MAX'],P['C'],P['SW_CNST_THK_BM'],P['T_CONST'],P['FLEX_RATIO'])
            else:
                print "ERROR! Invalid geometry type. Valid geometry types are:'"
                print "    'FP'"
                print "    'TD'"
                raise ValueError('ERROR! Invalid geometry type.')

    return (SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL)

def simulation_startup(P, DIO, PC, Swimmer, solid=None, FSI=None, PyFEA=None):
    """
    Handles starting a simulation. This will either start a new simulation or 
    start a simulation from a previous save.
    
    Args:
        P (dict):
        DIO (object):
        Swimmer (object):
        solid (object, optional):
        FSI (object, optional):
        PyFEA (object, optional):
        
    Returns:
        START_COUNTER (int):
        COUNTER (int):
        SwiL (list):
        GeoL (list):
        MotL (list):
        Swimmers (list):
        SolidL (list):
        FSIL (list):
        PyFEAL (list):
        
    Raises:
        ValueError: If given an unsupported 'START_FROM' keyword.
        ValueError: If input parameters are inconsistent with starting data.
    """
    
    # Check if the specified data output directory exists. If it does not, then
    # create a new output directory with that same name and for the simulation 
    # to start from zeroTime.
    if (os.path.exists(P['OUTPUT_DIR']) == False or os.listdir(P['OUTPUT_DIR']) == []):
        P['START_FROM'] = 'zeroTime'

    # If the keyword for starting the simulation is latestTime, look in the 
    # specified data output directory for the latest save state.
    if (P['START_FROM'] == 'latestTime'):
        startTime = 0.
        for file in os.listdir(''.join((P['OUTPUT_DIR'], '/'))):
            startTime = max(float(file), startTime)

        (sP, i, FLOWTIME, SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL) = DIO.read_data(''.join((P['OUTPUT_DIR'], '/', '%.8f' % startTime)))
        if not (sP['DEL_T'] == P['DEL_T']) and (sP['N_SWIMMERS'] == P['N_SWIMMERS']) and (sP['N_BODY'] == P['N_BODY']):
            raise ValueError('ERROR! Inconsistent input parameters with starting data file.')

        if (Swimmers[0].Wake.x.shape[0] < P['COUNTER']):
            for Swim in Swimmers:
                Swim.Wake.x.resize(P['COUNTER'])
                Swim.Wake.z.resize(P['COUNTER'])
                Swim.Wake.mu.resize(P['COUNTER']-1)
                Swim.Wake.gamma.resize(P['COUNTER'])

        START_COUNTER = i + 1
        COUNTER = P['COUNTER']

    # If the keyword for starting the simulation is firstTime, look in the 
    # specified data output directory for the earliest save state.
    elif (P['START_FROM'] == 'firstTime'):
        startTime = sys.float_info.max
        for file in os.listdir(''.join((P['OUTPUT_DIR'], '/'))):
            startTime = max(float(file), startTime)

        (sP, i, FLOWTIME, SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL) = DIO.read_data(''.join((P['OUTPUT_DIR'], '/', '%.8f' % startTime)))
        if not (sP['DEL_T'] == P['DEL_T']) and (sP['N_SWIMMERS'] == P['N_SWIMMERS']) and (sP['N_BODY'] == P['N_BODY']):
            raise ValueError('ERROR! Inconsistent input parameters with starting data file.')

        if (Swimmers[0].Wake.x.shape[0] < P['COUNTER']):
            for Swim in Swimmers:
                Swim.Wake.x.resize(P['COUNTER'])
                Swim.Wake.z.resize(P['COUNTER'])
                Swim.Wake.mu.resize(P['COUNTER']-1)
                Swim.Wake.gamma.resize(P['COUNTER'])

        START_COUNTER = i + 1
        COUNTER = P['COUNTER']

    # If the keyword for starting the simulation is zeroTime, call the geometry
    # setup function to create new objects for a fresh simulation.
    elif (P['START_FROM'] == 'zeroTime'):
        startTime = '0.00000000'
        (SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL) = geom_setup(P, PC, Swimmer, solid, FSI, PyFEA)

        START_COUNTER = 0
        COUNTER = P['COUNTER']

    # If none of those keywords existed, raise a value error and tell the user
    # acceptable values.
    else:
        print 'ERROR! Invalid START_FROM. Valid values are:'
        print '    latestTime'
        print '    firstTime'
        print '    zeroTime'
        raise ValueError('ERROR! Invalid START_FROM.')

    return (START_COUNTER, COUNTER, SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL)

def extrap1d(interpolator):
    """
    A wrapper for numpy scipy interpolators that allows for linear extrapolation
    beyond the provided interpolation set values.
    
    Args:
        interpolator (struct):
        
    Returns:
        ufunclike (struct):
    """
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike
    
def intermittent_ref(HEAVE_MAX, THETA_MAX, phi, DC, mag):
    """
    Reference intermittent signal function. This is used to interpolate the
    signal for a smooth signal in the intermittent function.
    
    Args:
        HEAVE_MAX (float):
        THETA_MAX (float):
        phi (float):
        DC (float):
        mag (float):
        
    Returns:
        t_T (float):
        y_pitch (float):
        y_pitch_plus (float):
        y_pitch_minus (float):
        y_heave (float):
        y_heave_plus (float):
        y_heave_minus (float):
    """
    # Constants
    N = 5000      # Defines the number of points for the active and passive portions of the reference signal
    Tstep = 1e-5
    # Defining a smooth intermittent signal
    # Modifying the slope of the envelope signal as DC --> 1
    if (DC > 0.9):
        a = mag / (1 - DC)
    else:
        a = 10 * mag
    
    # Defining the non-dimenional time for the active and passive portions of
    # the cycle. t_T_active is the non-dimensional time over the active portion, 
    # t is normalized by the active period. t_T_passive is the non-dimensional 
    # time over the passive portion, t is normalized by the active period.
    t_T_active = np.linspace(0.0, 1.0, N).T
    if (DC < 1.0):
        t_T_passive = np.linspace(1.0, 1. / DC, N).T
#        t_T_passive = np.copy(t_T_passive[1:]) # Removing overlapping point
    else:
        t_T_passive = np.empty((N,1)) # For DC = 1 the passive time vector must be a null vector 
    
    # Creating the amplitude envelope and modifying the active portion of the
    # signal
    amp_env = -np.tanh(a * t_T_active) * np.tanh(a * (t_T_active - 1.0))
    y_pitch_sin       = THETA_MAX * np.sin(2. * np.pi * t_T_active           + DC * phi)
    y_pitch_sin_plus  = THETA_MAX * np.sin(2. * np.pi * (t_T_active + Tstep) + DC * phi)
    y_pitch_sin_minus = THETA_MAX * np.sin(2. * np.pi * (t_T_active - Tstep) + DC * phi)
    
    y_heave_sin       = HEAVE_MAX * np.sin(2. * np.pi * t_T_active          )
    y_heave_sin_plus  = HEAVE_MAX * np.sin(2. * np.pi * (t_T_active + Tstep))
    y_heave_sin_minus = HEAVE_MAX * np.sin(2. * np.pi * (t_T_active - Tstep))
    
    sig_pitch_mod       = amp_env * y_pitch_sin
    sig_pitch_mod_plus  = amp_env * y_pitch_sin_plus
    sig_pitch_mod_minus = amp_env * y_pitch_sin_minus
    
    sig_heave_mod       = amp_env * y_heave_sin
    sig_heave_mod_plus  = amp_env * y_heave_sin_plus
    sig_heave_mod_minus = amp_env * y_heave_sin_minus
    
    # Assembling an intermittent reference signal for one full cycle
    y_pitch       = np.hstack((sig_pitch_mod,       np.zeros(t_T_passive.size-1).T))
    y_pitch_plus  = np.hstack((sig_pitch_mod_plus,  np.zeros(t_T_passive.size-1).T))
    y_pitch_minus = np.hstack((sig_pitch_mod_minus, np.zeros(t_T_passive.size-1).T))
    
    y_heave       = np.hstack((sig_heave_mod,       np.zeros(t_T_passive.size-1).T))
    y_heave_plus  = np.hstack((sig_heave_mod_plus,  np.zeros(t_T_passive.size-1).T))
    y_heave_minus = np.hstack((sig_heave_mod_minus, np.zeros(t_T_passive.size-1).T))
    
    t_T = np.hstack((t_T_active, t_T_passive[1:]))
    
    return(t_T, y_pitch, y_pitch_plus, y_pitch_minus, y_heave, y_heave_plus, y_heave_minus)
    
def intermittent(HEAVE_MAX, THETA_MAX, phi, DC, f, N_STEP, N_CYC, s):
    """
    Creates an intermittent kinematic signal.
    
    Args:
        HEAAVE_MAX (float):
        THETA_MAX (float):
        phi (float):
        DC (float):
        f (float):
        N_STEP (float):
        N_CYC (float):
        s (float):
        
    Returns:
        angle_pitch (float):
        phase_heave (float):
        period (float):
    """
    # Constants
    mag = 3 # defines the slope of the hyperbolic tangent smoothing function inside of intermittent_ref
    
    # Defining an intermittent signal
    (t_T_ref, y_pitch_ref, y_pitch_ref_plus, y_pitch_ref_minus, y_heave_ref, 
     y_heave_ref_plus, y_heave_ref_minus) = intermittent_ref(HEAVE_MAX, THETA_MAX, phi, DC, mag)
    
    # Defining the time step spacing, delT/T_active, that is the time step
    # normalized by the period of the active motion, which is T_active = 1/f
    delT_T = 1.0 / DC / N_STEP / f
    
    # Creating time signal
    t_T_single = np.arange(0., 1. / DC / f, delT_T).T

    period = t_T_single[-1]
    t_actual = t_T_ref / f
    
    # Sampling the pitch reference signal at t_T    
    y_pitch_single       = PchipInterpolator(t_actual, y_pitch_ref      )(t_T_single)
    y_pitch_single_plus  = PchipInterpolator(t_actual, y_pitch_ref_plus )(t_T_single)
    y_pitch_single_minus = PchipInterpolator(t_actual, y_pitch_ref_minus)(t_T_single)
    
    # Sampling the heave reference signal at t_T
    y_heave_single       = PchipInterpolator(t_actual, y_heave_ref      )(t_T_single)
    y_heave_single_plus  = PchipInterpolator(t_actual, y_heave_ref_plus )(t_T_single)
    y_heave_single_minus = PchipInterpolator(t_actual, y_heave_ref_minus)(t_T_single)
    
    # Copying the signal for Ncyc cycles
    if (s == 0):
        angle_pitch = np.hstack((y_pitch_single, np.tile(y_pitch_single, (N_CYC-1)))) # Full pitch signal for intermittent swimming
    elif (s == 1):
        angle_pitch = np.hstack((y_pitch_single_plus, np.tile(y_pitch_single_plus, (N_CYC-1))))
    elif (s == -1):
        angle_pitch = np.hstack((y_pitch_single_minus, np.tile(y_pitch_single_minus, (N_CYC-1))))
    
    if (s==0):
        phase_heave = np.hstack((y_heave_single, np.tile(y_heave_single, (N_CYC-1)))) # Full heave signal for intermittent swimming
    elif (s==1):
        phase_heave = np.hstack((y_heave_single_plus, np.tile(y_heave_single_plus, (N_CYC-1))))
    elif (s==-1):
        phase_heave = np.hstack((y_heave_single_minus, np.tile(y_heave_single_minus, (N_CYC-1))))
    
    return(angle_pitch, phase_heave, period)

def squarewave_ref(a):
    """
    This function creates an squarewave reference signal that needs to be
    sampled at specific times, scaled, and shifted to construct the actual  
    signal to be used in the computations.
    
    Args:
        a (float): Defines the slope of the hyperbolic tangent function.
        
    Returns:
        t_T (float):
        y (float):
    """
    # Constants
    N = 1e5 # Defines the number of points for the half-period reference signal
    
    # Defining an smooth squarewave signal
    # Defining the non-dimensional time for the cycle
    t_T_h = np.linspace(0., N, 0.5).T # non-dimensional time, t is normalized by the period.
    
    # Creating the squarewave signal
    y_h1 = -np.tanh(a * t_T_h) * np.tanh(a * (t_T_h - 0.5))
    y_h2 = -(np.flipud(y_h1)) # Mirror and invert first half-wave to create second half
    
    # Assembling full wave
    y_full = np.vstack((y_h1, y_h2[1:])) # Combine first and second half to create full wave
    
    # Final reference signal
    y = (1. / np.max(y_full)) * y_full # Normalize the y-values to -1 <= y <= 1
    t_T = np.vstack((t_T_h, t_T_h[1:] + 0.5)) # Create a full set of x values (0 to 1)
    
    return(t_T, y)

def squarewave(f, A, a, N_STEP, N_CYC):
    """
    Args:
        f (float):
        A (float):
        a (float):
        N_STEP (int):
        N_CYC (int):
        
    Returns:
        t (float):
        y (float):
        t_T_ref (float):
        y_ref (float):
    """
    (t_T_ref, y_ref) = squarewave_ref(a)
    
    # Defining the time step spacing, that is the time step
    # normalized by the period of the motion, which is T = 1/f
    delT_T = 1. / N_STEP
    
    # Creating time signal for a single cycle
    t_T_single = np.linspace(0., 1., delT_T).T
    
    # Sampling the reference signal at t_T
    y_single = PchipInterpolator(t_T_ref, y_ref)(t_T_single)
    
    # Copying the signal for Ncyc cycles
    y = A * np.hstack((y_single, np.tile(y_single[1:], (N_CYC-1,1)))) # Full signal vector for squarewave signal ***************             
    
    # Calculate a shift to the time signal for each extra cycle beyond the first
    shift  = t_T_single[-1]
    shift_vec = np.tile(shift*np.arange(1,N_CYC), (np.size(t_T_single[1:]), 1))
    (m, n) = np.shape(shift_vec)
    shift_vec = shift_vec.reshape(m*n, 1,order='F').copy()
    
    t = 1. / f * np.vstack((t_T_single, shift_vec + np.tile(t_T_single[1:], (N_CYC-1, 1)))) # Full time vector for squarewave signal ***************

    return(t, y, t_T_ref, y_ref)

def multi_kinematics(P, PHI=0., scale=None, rate=50):
    """
    Creates a sine, square, triangle, or sawtooth kinematic signal. The signal 
    forms can be weighted for any superposition of the signals.
    
    Args:
        P (dict):
        PHI (float, optional):
        scale (float, optional):
        rate (float, optional):
    
    Returns:
        signal (float):
        signalMinus (float):
        signalPlus (float):
    """
    delta = 1. / rate
    if (scale == None):
        x = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        y = [2.1790849673, 1.5669799009, 1.3790259591, 1.2876797596, 1.2332922318, 1.1968681428, 1.1709601874, 1.1512761466, 1.136077986, 1.1238710886]
        scaleSig = PchipInterpolator(x, y)
        scale = [scaleSig(rate), scaleSig(rate), scaleSig(rate), scaleSig(rate)]
    
    # Sine Wave Input
    sineFrac      = [scale[0] * np.sin(2 * np.pi * P['F'] * P['T'][i]                + PHI) for i in xrange(P['COUNTER'])]
    sineFracMinus = [scale[0] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + PHI) for i in xrange(P['COUNTER'])]
    sineFracPlus  = [scale[0] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + PHI) for i in xrange(P['COUNTER'])]
    
    # Square Wave Input
    a = 9.117
    squareFrac      = [np.tanh(a * P['T'][i]               ) for i in xrange(P['COUNTER'])]
    squareFracMinus = [np.tanh(a * (P['T'][i] - P['TSTEP'])) for i in xrange(P['COUNTER'])]
    squareFracPlus  = [np.tanh(a * (P['T'][i] + P['TSTEP'])) for i in xrange(P['COUNTER'])]
    
    for j in [i+1 for i in xrange(P['N_CYC'])]:
        squareFrac      = [squareFrac[i]      * np.tanh(a * (P['T'][i]                - (j-0.5)/P['F'])) * np.tanh(a * (P['T'][i]                - j / P['F'])) for i in xrange(P['COUNTER'])]
        squareFracMinus = [squareFracMinus[i] * np.tanh(a * ((P['T'][i] - P['TSTEP']) - (j-0.5)/P['F'])) * np.tanh(a * ((P['T'][i] - P['TSTEP']) - j / P['F'])) for i in xrange(P['COUNTER'])]
        squareFracPlus  = [squareFracPlus[i]  * np.tanh(a * ((P['T'][i] + P['TSTEP']) - (j-0.5)/P['F'])) * np.tanh(a * ((P['T'][i] + P['TSTEP']) - j / P['F'])) for i in xrange(P['COUNTER'])]
    
#    squareFrac      = [scale[1] * 2. * np.arctan(np.sin(2. * np.pi * P['F'] * P['T'][i]                + PHI) / delta) / np.pi for i in xrange(P['COUNTER'])]
#    squareFracMinus = [scale[1] * 2. * np.arctan(np.sin(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + PHI) / delta) / np.pi for i in xrange(P['COUNTER'])]
#    squareFracPlus  = [scale[1] * 2. * np.arctan(np.sin(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + PHI) / delta) / np.pi for i in xrange(P['COUNTER'])]
    
    # Triangle wave
    triangleFrac      = [scale[2] * (1. - 2. * np.arccos((1. - delta) * np.sin(2. * np.pi * P['F'] * P['T'][i]                + PHI)) / np.pi) for i in xrange(P['COUNTER'])]
    triangleFracMinus = [scale[2] * (1. - 2. * np.arccos((1. - delta) * np.sin(2. * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + PHI)) / np.pi) for i in xrange(P['COUNTER'])]
    triangleFracPlus  = [scale[2] * (1. - 2. * np.arccos((1. - delta) * np.sin(2. * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + PHI)) / np.pi) for i in xrange(P['COUNTER'])]
    
    # Saw-tooth Wave
    xt  = [0.25 * (2. * P['F'] * P['T'][i]                - 1.) for i in xrange(P['COUNTER'])]
    xtm = [0.25 * (2. * P['F'] * (P['T'][i] - P['TSTEP']) - 1.) for i in xrange(P['COUNTER'])]
    xtp = [0.25 * (2. * P['F'] * (P['T'][i] + P['TSTEP']) - 1.) for i in xrange(P['COUNTER'])]
    
    xs  = [0.5 * P['F'] * P['T'][i]                for i in xrange(P['COUNTER'])]
    xsm = [0.5 * P['F'] * (P['T'][i] - P['TSTEP']) for i in xrange(P['COUNTER'])]
    xsp = [0.5 * P['F'] * (P['T'][i] + P['TSTEP']) for i in xrange(P['COUNTER'])]
    
    trg  = [1. - 2. * np.arccos((1. - delta) * np.sin(2. * np.pi * xt[i] )) / np.pi for i in xrange(P['COUNTER'])]
    trgM = [1. - 2. * np.arccos((1. - delta) * np.sin(2. * np.pi * xtm[i])) / np.pi for i in xrange(P['COUNTER'])]
    trgP = [1. - 2. * np.arccos((1. - delta) * np.sin(2. * np.pi * xtp[i])) / np.pi for i in xrange(P['COUNTER'])]
    
    sqr  = [2. * np.arctan(np.sin(2. * np.pi * xs[i] ) / delta) / np.pi for i in xrange(P['COUNTER'])]
    sqrM = [2. * np.arctan(np.sin(2. * np.pi * xsm[i]) / delta) / np.pi for i in xrange(P['COUNTER'])]
    sqrP = [2. * np.arctan(np.sin(2. * np.pi * xsp[i]) / delta) / np.pi for i in xrange(P['COUNTER'])]
    
    sawFrac      = [scale[3] * trg[i]  * sqr[i]  for i in xrange(P['COUNTER'])]
    sawFracMinus = [scale[3] * trgM[i] * sqrM[i] for i in xrange(P['COUNTER'])]
    sawFracPlus  = [scale[3] * trgP[i] * sqrP[i] for i in xrange(P['COUNTER'])]

    # Form composite signal
    signal      = [P['SIG_WEIGHT'][0] * sineFrac[i]      + P['SIG_WEIGHT'][1] * squareFrac[i]      + P['SIG_WEIGHT'][2] * triangleFrac[i]      + P['SIG_WEIGHT'][3] * sawFrac[i]      for i in xrange(P['COUNTER'])]
    signalMinus = [P['SIG_WEIGHT'][0] * sineFracMinus[i] + P['SIG_WEIGHT'][1] * squareFracMinus[i] + P['SIG_WEIGHT'][2] * triangleFracMinus[i] + P['SIG_WEIGHT'][3] * sawFracMinus[i] for i in xrange(P['COUNTER'])]
    signalPlus  = [P['SIG_WEIGHT'][0] * sineFracPlus[i]  + P['SIG_WEIGHT'][1] * squareFracPlus[i]  + P['SIG_WEIGHT'][2] * triangleFracPlus[i]  + P['SIG_WEIGHT'][3] * sawFracPlus[i]  for i in xrange(P['COUNTER'])]
    
    return(signal, signalMinus, signalPlus)

def vel_multi_kinematics(P, sig):
    """
    Creates an velcoity signal based on an input kinematic signal using 
    fourth-order finite difference stencils (forward, backward, and central).
    
    Args:
        P (dict):
        sig (list, str):
    
    Returns:
        accell (float):
        accell_minus (float):
        accell_plus (float):
    
    """
    DEL_T = P['DEL_T']
    signal       = sig[0]
    signal_minus = sig[1]
    signal_plus  = sig[2]
    vel       = [0. for i in xrange(P['COUNTER'])]
    vel_minus = [0. for i in xrange(P['COUNTER'])]
    vel_plus  = [0. for i in xrange(P['COUNTER'])]
    
    for i in [0, 1]:
        vel[i]       = (-25. / 12.) * signal[i]       + 4. * signal[i+1]       - 3. * signal[i+2] + (4. / 3.) * signal[i+3]       - (1. / 4.) * signal[i+4]
        vel_minus[i] = (-25. / 12.) * signal_minus[i] + 4. * signal_minus[i+1] - 3. * signal[i+2] + (4. / 3.) * signal_minus[i+3] - (1. / 4.) * signal_minus[i+4]
        vel_plus[i]  = (-25. / 12.) * signal_plus[i]  + 4. * signal_plus[i+1]  - 3. * signal[i+2] + (4. / 3.) * signal_plus[i+3]  - (1. / 4.) * signal_plus[i+4]

    for i in xrange(P['COUNTER']-4):
        vel[i+2]       = (1. / 12.) * signal[i]       - (2. / 3.) * signal[i+1]       + (2. / 3.) * signal[i+3]       - (1. / 12.) * signal[i+4]
        vel_minus[i+2] = (1. / 12.) * signal_minus[i] - (2. / 3.) * signal_minus[i+1] + (2. / 3.) * signal_minus[i+3] - (1. / 12.) * signal_minus[i+4]
        vel_plus[i+2]  = (1. / 12.) * signal_plus[i]  - (2. / 3.) * signal_plus[i+1]  + (2. / 3.) * signal_plus[i+3]  - (1. / 12.) * signal_plus[i+4]
        
    for i in [-1, -2]:
        vel[i]       = (25. / 12.) * signal[i]       - 4. * signal[i-1]       + 3. * signal[i-2]       - (4. / 3.) * signal[i-3]       + (1. / 4.) * signal[i-4]
        vel_minus[i] = (25. / 12.) * signal_minus[i] - 4. * signal_minus[i-1] + 3. * signal_minus[i-2] - (4. / 3.) * signal_minus[i-3] + (1. / 4.) * signal_minus[i-4]
        vel_plus[i]  = (25. / 12.) * signal_plus[i]  - 4. * signal_plus[i-1]  + 3. * signal_plus[i-2]  - (4. / 3.) * signal_plus[i-3]  + (1. / 4.) * signal_plus[i-4]

    vel       = [vel[i]       / DEL_T for i in xrange(P['COUNTER'])]
    vel_minus = [vel_minus[i] / DEL_T for i in xrange(P['COUNTER'])]
    vel_plus  = [vel_plus[i]  / DEL_T for i in xrange(P['COUNTER'])]

    return(vel, vel_minus, vel_plus)

def accel_multi_kinematics(P, sig):
    """
    Creates an accelleration signal based on an input kinematic signal using 
    fourth-order finite difference stencils (forward, backward, and central).
    
    Args:
        P (dict):
        sig (list, str):
    
    Returns:
        accell (float):
        accell_minus (float):
        accell_plus (float):
    
    """
    DEL_T = P['DEL_T']
    signal       = sig[0]
    signal_minus = sig[1]
    signal_plus  = sig[2]
    accell       = [0. for i in xrange(P['COUNTER'])]
    accell_minus = [0. for i in xrange(P['COUNTER'])]
    accell_plus  = [0. for i in xrange(P['COUNTER'])]
    
    for i in [0, 1]:
        accell[i]       = (15. / 4.) * signal[i]       - (77. / 6.) * signal[i+1]       + (107. / 6.) * signal[i+2] - 13. * signal[i+3]       + (61. / 12.) * signal[i+4]       - (5. / 6.) * signal[i+5]
        accell_minus[i] = (15. / 4.) * signal_minus[i] - (77. / 6.) * signal_minus[i+1] + (107. / 6.) * signal[i+2] - 13. * signal_minus[i+3] + (61. / 12.) * signal_minus[i+4] - (5. / 6.) * signal_minus[i+5]
        accell_plus[i]  = (15. / 4.) * signal_plus[i]  - (77. / 6.) * signal_plus[i+1]  + (107. / 6.) * signal[i+2] - 13. * signal_plus[i+3]  + (61. / 12.) * signal_plus[i+4]  - (5. / 6.) * signal_plus[i+5]

    for i in xrange(P['COUNTER']-4):
        accell[i+2]       = (-1. / 12.) * signal[i]       + (4. / 3.) * signal[i+1]       - (5. / 2.) * signal[i+2]       + (4. / 3.) * signal[i+3]       - (1. / 12.) * signal[i+4]
        accell_minus[i+2] = (-1. / 12.) * signal_minus[i] + (4. / 3.) * signal_minus[i+1] - (5. / 2.) * signal_minus[i+2] + (4. / 3.) * signal_minus[i+3] - (1. / 12.) * signal_minus[i+4]
        accell_plus[i+2]  = (-1. / 12.) * signal_plus[i]  + (4. / 3.) * signal_plus[i+1]  - (5. / 2.) * signal_plus[i+2]  + (4. / 3.) * signal_plus[i+3]  - (1. / 12.) * signal_plus[i+4]
        
    for i in [-1, -2]:
        accell[i]       = (15. / 4.) * signal[i]       - (77. / 6.) * signal[i-1]       + (107. / 6.) * signal[i-2]       - 13. * signal[i-3]       + (61. / 12.) * signal[i-4]       - (5. / 6.) * signal[i-5]
        accell_minus[i] = (15. / 4.) * signal_minus[i] - (77. / 6.) * signal_minus[i-1] + (107. / 6.) * signal_minus[i-2] - 13. * signal_minus[i-3] + (61. / 12.) * signal_minus[i-4] - (5. / 6.) * signal_minus[i-5]
        accell_plus[i]  = (15. / 4.) * signal_plus[i]  - (77. / 6.) * signal_plus[i-1]  + (107. / 6.) * signal_plus[i-2]  - 13. * signal_plus[i-3]  + (61. / 12.) * signal_plus[i-4]  - (5. / 6.) * signal_plus[i-5]

    accell       = [accell[i]       / DEL_T**2 for i in xrange(P['COUNTER'])]
    accell_minus = [accell_minus[i] / DEL_T**2 for i in xrange(P['COUNTER'])]
    accell_plus  = [accell_plus[i]  / DEL_T**2 for i in xrange(P['COUNTER'])]

    return(accell, accell_minus, accell_plus)