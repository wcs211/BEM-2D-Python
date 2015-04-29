import numpy as np
       
    # x,z components of each panel's tangential and normal vectors
def panel_vectors(x,z):
    lpanel = np.sqrt((x[1:]-x[:-1])**2 + (z[1:]-z[:-1])**2)
    tx = (x[1:]-x[:-1])/lpanel
    tz = (z[1:]-z[:-1])/lpanel
    nx = -tz
    nz = tx
    return (tx,tz,nx,nz,lpanel)
    
    # x,z components of each midpoint's/collocation point's tangential and normal vectors
def point_vectors(xdp,xdm,zdp,zdm):
    tx = (xdp-xdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    tz = (zdp-zdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    nx = -tz
    nz = tx
    return(tx,tz,nx,nz)
    
def archive(array, axis=0):
    """
    Shifts array values along an axis (row-wise by default).
    Used for arrays that keep past values for differencing with respect to time.
    """
    
    if axis == 0:
        array[1:,:] = array[:-1,:]
    elif axis == 1:
        array[:,1:] = array[:,:-1]