import numpy as np
       
    #x,z components of each panel's tangential and normal vectors
def PanelVectors(x,z):
    lpanel=np.sqrt((x[1:]-x[:-1])**2 + (z[1:]-z[:-1])**2)
    tx=(x[1:]-x[:-1])/lpanel
    tz=(z[1:]-z[:-1])/lpanel
    nx=-tz
    nz=tx
    return (tx,tz,nx,nz,lpanel)
    
    #x,z components of each collacation point's tangential and normal vectors
def PointVectors(xdp,xdm,zdp,zdm):
    txc = (xdp-xdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    tzc = (zdp-zdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    nxc = -tzc
    nzc = txc
    return(txc,tzc,nxc,nzc)