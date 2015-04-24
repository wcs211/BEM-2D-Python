import numpy as np
    
# VandeVooren airfoil type
def van_de_vooren(Body):
# Uses Body.(C, EPSILON, K, N)
# Gets Body.(xb, zb, xb_col, zb_col)
# Others: A, phi, r1, r2, phi1, phi2
    
    A = Body.C*((1+Body.EPSILON)**(Body.K-1))*(2**(-Body.K))
    
    phi = np.linspace(0,np.pi,Body.N/2+1)
    
    r1 = np.sqrt((A*np.cos(phi)-A)**2+(A**2)*np.sin(phi)**2)
    r2 = np.sqrt((A*np.cos(phi)-Body.EPSILON*A)**2+(A**2)*np.sin(phi)**2)
    
    phi1 = np.arctan2((A*np.sin(phi)) , (A*np.cos(phi)-A))
    phi2 = np.arctan2(A*np.sin(phi) ,(A*np.cos(phi)-Body.EPSILON*A))
    
    xb = ((r1**Body.K)/(r2**(Body.K-1)))*(np.cos(Body.K*phi1)*np.cos((Body.K-1)*phi2) + np.sin(Body.K*phi1)*np.sin((Body.K-1)*phi2))
    zb_top = ((r1**Body.K)/(r2**(Body.K-1)))*(np.sin(Body.K*phi1)*np.cos((Body.K-1)*phi2) - np.cos(Body.K*phi1)*np.sin((Body.K-1)*phi2))
    zb_bot = -((r1**Body.K)/(r2**(Body.K-1)))*(np.sin(Body.K*phi1)*np.cos((Body.K-1)*phi2) - np.cos(Body.K*phi1)*np.sin((Body.K-1)*phi2))
    
    xb = xb-xb[-1] # Carrying the leading edge to the origin
    xb[0] = Body.C
    
    zb_top[0] = 0
    zb_bot[0] = 0
    zb_bot[-1] = 0
    
    # Merge top and bottom surfaces together
    Body.xb = np.hstack((xb , xb[-2::-1]))
    Body.zb = np.hstack((zb_bot , zb_top[-2::-1]))
    
    Body.xb_col = ((Body.xb[1:] + Body.xb[:-1])/2)
    Body.zb_col = ((Body.zb[1:] + Body.zb[:-1])/2)

#Flat plate geometry
def flat_plate(Body):
    
    Body.xf = np.linspace(1,0,Body.N/2+1)
    Body.zf = np.zeros(Body.N/2+1)