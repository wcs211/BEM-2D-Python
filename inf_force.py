import numpy as np
from general_functions import panel_vectors

# Velocity and velocity potential equations are defined in panel coordinates so a transformation should be done
# Each row of xpleft/xpright/zp is an influence, and each column is a target
# NI is N influences, NT is N targets
# xi/zi is x/z of influences, xt/zt is x/z of target points
def transformation(xt,zt,xi,zi):
# Returns xp1, xp2, zp
# Others: NT, NI, tx, tz, nx, nz, dummy, x_tile, z_tile, tx_tile, tz_tile
    
    NT = np.size(xt)
    NI = np.size(xi)-1
    
    (tx,tz,nx,nz) = panel_vectors(xi,zi)[:-1]
    
    # Intermediary variables to reduce number of tile/repeat operations
    # From normalvectors: tx==nz, tz==-nx
    x_tile = np.repeat(xt[:,np.newaxis].T,NI,0) - np.repeat(xi[:-1,np.newaxis],NT,1)
    z_tile = np.repeat(zt[:,np.newaxis].T,NI,0) - np.repeat(zi[:-1,np.newaxis],NT,1)
    tx_tile = np.repeat(tx[:,np.newaxis],NT,1)
    tz_tile = np.repeat(tz[:,np.newaxis],NT,1)
    
    # Transforming left side collocation points from global to local coordinates
    xp1 = x_tile*tx_tile + z_tile*tz_tile
    zp = x_tile*(-tz_tile) + z_tile*tx_tile
    
    # Transforming right side panel points into local coordinate system
    dummy = (xi[1:]-xi[:-1])*tx + (zi[1:]-zi[:-1])*tz
    xp2 = xp1 - np.repeat(dummy[:,np.newaxis],NT,1)
    
    return(xp1,xp2,zp)

# This method constructs the influence coefficent matrices
def influence_matrices(Body,Edge,Wake,i):
# Uses transformation()
# Uses Body.(N, x, z, x_col, z_col, V0, vx, vz)
# Uses Edge.(x, z)
# Uses Wake.(x, z)
# Gets Body.(phi_s, sigma, phi_db, phi_dw)
# Others: xp1, xp2, zp, nx, nz, tx, tz, xp1_e, xp2_e, zp_e
    
    if i > 0:
        # Tangential and normal body panels and panel length calculations
        (xp1,xp2,zp) = transformation(Body.AFC.x_col,Body.AFC.z_col,Body.AFC.x,Body.AFC.z)
    
        # Body source singularities influencing the body
        # Transpose so that row elements represent the effect on the (row number)th panel
        Body.phi_s = np.transpose((xp1 * np.log(xp1**2 + zp**2) - xp2 * np.log(xp2**2 + zp**2) \
                                  + 2*zp*(np.arctan2(zp,xp2) - np.arctan2(zp,xp1)))/(4*np.pi))
        
        # Body source strength calculations
        (nx,nz) = panel_vectors(Body.AFC.x,Body.AFC.z)[2:4]
        Body.sigma = nx*(Body.V0 + Body.vx) + nz*Body.vz   # normal vector pointing outward(overall sigma pointing outward)
        
        # Body doublet singularities influencing body itself
        # Transpose similar to phi_s
        Body.phi_db = np.transpose(-(np.arctan2(zp,xp2)\
                                  - np.arctan2(zp,xp1))/(2*np.pi))
        
        # Edge doublet influencing the body
        (xp1,xp2,zp) = transformation(Body.AFC.x_col,Body.AFC.z_col,Edge.x,Edge.z)
        
        if i > 1: # No wake panels until i==2
            # Wake doublets influencing the body
            (xp1_w,xp2_w,zp_w) = transformation(Body.AFC.x_col,Body.AFC.z_col,Wake.x[:i],Wake.z[:i])
            # Join edge and wake doublet influences into a single matrix
            xp1 = np.insert(xp1, Edge.N, xp1_w, axis=0)
            xp2 = np.insert(xp2, Edge.N, xp2_w, axis=0)
            zp = np.insert(zp, Edge.N, zp_w, axis=0)
    
        Body.phi_dw = np.transpose(-(np.arctan2(zp,xp2) \
                                   -np.arctan2(zp,xp1))/(2*np.pi))
        
def kutta(Body,Edge,Wake,RHO,i,DEL_T,SWITCH_KUTTA):
# Uses Body.(N, phi_dw, phi_s, sigma, phi_db), Wake.mu, pressure()
# Gets Edge.(mu, gamma), Body.(phi_dinv, mu, mu_past, gamma)
# Others: n_iter, c, mu_guess, delta_cp, rhs, slope
    
    if i == 0:
        pass
    else:
        n_iter = 0
        mu_guess = np.empty(2) # [0] is current guess, [1] is previous
        delta_cp = np.empty(2) # [0] is current delta_cp, [1] is previous
        
        while True:
            n_iter += 1
            
            if n_iter == 1:
                # Begin with explicit Kutta condition as first guess
                # Construct the augmented body matrix by combining body and trailing edge panel arrays
                c = np.zeros((Body.N,Body.N))
                c[:,0] = -Body.phi_dw[:,0] # c contains the explicit Kutta condition (first and last columns)
                c[:,-1] = Body.phi_dw[:,0]
                Body.phi_dinv = np.linalg.inv(Body.phi_db + c)
                # Get rhs
                if i == 1:
                    rhs = -np.dot(Body.phi_s,Body.sigma)
                else:
                    rhs = -np.dot(Body.phi_s,Body.sigma) - np.dot(Body.phi_dw[:,1:],Wake.mu[:i-1])
                # Solve for body doublet strengths using explicit Kutta
                Body.mu = np.dot(Body.phi_dinv,rhs)
                # First mu_guess (from explicit Kutta)
                mu_guess[0] = Body.mu[-1]-Body.mu[0]
                
            else:
                if n_iter == 2: # Make a second initial guess
                    # Update phi_dinv so it no longer includes explicit Kutta condition
                    Body.phi_dinv = np.linalg.inv(Body.phi_db)
                    
                    mu_guess[1] = mu_guess[0]
                    delta_cp[1] = delta_cp[0]
                    mu_guess[0] = 0.8*mu_guess[1] # Multiply first (explicit) guess by arbitrary constant to get second guess
                    
                else: # Newton method to get delta_cp == 0
                    # Get slope, which is change in delta_cp divided by change in mu_guess
                    slope = (delta_cp[0]-delta_cp[1])/(mu_guess[0]-mu_guess[1])
                    mu_guess[1] = mu_guess[0]
                    delta_cp[1] = delta_cp[0]
                    mu_guess[0] = mu_guess[1] - delta_cp[0]/slope
                
                # Form right-hand side including mu_guess as an influence
                if i == 1:
                    rhs = -np.dot(Body.phi_s,Body.sigma) - np.squeeze(np.dot(Body.phi_dw,mu_guess[0]))
                else:              
                    rhs = -np.dot(Body.phi_s,Body.sigma) - np.dot(Body.phi_dw,np.insert(Wake.mu[:i-1],0,mu_guess[0]))
                
                Body.mu = np.dot(Body.phi_dinv,rhs)
                
            
            pressure(Body,RHO,i,DEL_T)
            if SWITCH_KUTTA == 0:
                break
            delta_cp[0] = np.absolute(Body.cp[-1]-Body.cp[0])
            if delta_cp[0] < 0.0001:
                break
                    
        # mu_past used in differencing for pressure
        Body.mu_past[1,:] = Body.mu_past[0,:]
        Body.mu_past[0,:] = Body.mu
        
        Edge.mu = mu_guess[0]
        Edge.gamma[0] = -Edge.mu
        Edge.gamma[1] = Edge.mu
        
        # Get gamma of body panels for use in wake rollup
        Body.gamma[0] = -Body.mu[0]
        Body.gamma[1:-1] = Body.mu[:-1]-Body.mu[1:]
        Body.gamma[-1] = Body.mu[-1]

def pressure(Body,RHO,i,DEL_T):
# Uses Body.(N, x, z, V0, vx, vz, sigma, mu, mu_past), panel_vectors()
# Gets Body.cp and Body.p
# Others: dmu_dl, dmu_dt, qpx_tot, qpz_tot
    
    if i > 0:
        
        (tx,tz,nx,nz,lpanel) = panel_vectors(Body.AFC.x,Body.AFC.z)
        
        # Tangential panel velocity dmu/dl, first-order differencing
        dmu_dl = np.empty(Body.N)
        dmu_dl[0] = (Body.mu[0]-Body.mu[1]) / (lpanel[0]/2 + lpanel[1]/2)
        dmu_dl[1:-1] = (Body.mu[:-2]-Body.mu[2:]) / (lpanel[:-2]/2 + lpanel[1:-1] + lpanel[2:]/2)
        dmu_dl[-1] = (Body.mu[-2]-Body.mu[-1]) / (lpanel[-2]/2 + lpanel[-1]/2)
        
        # Potential change dmu/dt, second-order differencing after first time step
        if i == 1:
            dmu_dt = (Body.mu - Body.mu_past[0,:])/DEL_T
        else:
            dmu_dt = (3*Body.mu - 4*Body.mu_past[0,:] + Body.mu_past[1,:])/(2*DEL_T)
        
        # Unsteady pressure calculation (from Matlab code)
        qpx_tot = dmu_dl*tx + Body.sigma*nx
        qpz_tot = dmu_dl*tz + Body.sigma*nz
        
        Body.p = -RHO*(qpx_tot**2 + qpz_tot**2)/2. + RHO*dmu_dt + RHO*(qpx_tot*(Body.V0+Body.vx) + qpz_tot*Body.vz)
        Body.cp = Body.p / (0.5*RHO*Body.V0**2)
       
# Calculation of drag and lift forces affecting overall airfoil      
def force(Body,i):
# Uses Body.(x, z, p, N), panel_vectors()
# Gets Body.(drag, lift)
# Others: tx, tz, nx, nz, lpanel
    
    if i == 0:
        pass

    else:
        (tx,tz,nx,nz,lpanel) = panel_vectors(Body.x,Body.z)
                              
        Body.drag[i-1] = np.dot(Body.p[i-1,:]*lpanel, np.reshape(tx,(Body.N,1)))\
                      + np.dot(Body.p[i-1,:]*lpanel, np.reshape(-tz,(Body.N,1)))

        Body.lift[i-1] = np.dot(Body.p[i-1,:]*lpanel, np.reshape(-nz,(Body.N,1)))\
                      + np.dot(Body.p[i-1,:]*lpanel, np.reshape(nx,(Body.N,1)))