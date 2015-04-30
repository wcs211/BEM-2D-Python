import numpy as np
from general_functions import panel_vectors, transformation

# This method constructs the influence coefficent matrices
def influence_matrices(Body, Edge, Wake, i):
# Uses transformation()
# Uses Body.(N, x, z, x_col, z_col, V0, vx, vz)
# Uses Edge.(x, z)
# Uses Wake.(x, z)
# Gets Body.(phi_s, sigma, phi_db, phi_dw)
# Others: xp1, xp2, zp, nx, nz, tx, tz, xp1_e, xp2_e, zp_e
    
    if i > 0:
        # Tangential and normal body panels and panel length calculations
        (xp1,xp2,zp) = transformation(Body.AF.x_col,Body.AF.z_col,Body.AF.x,Body.AF.z)
    
        # Body source singularities influencing the body
        # Transpose so that row elements represent the effect on the (row number)th panel
        Body.phi_s = np.transpose((xp1 * np.log(xp1**2 + zp**2) - xp2 * np.log(xp2**2 + zp**2) \
                                  + 2*zp*(np.arctan2(zp,xp2) - np.arctan2(zp,xp1)))/(4*np.pi))
        
        # Body source strength calculations
        (nx,nz) = panel_vectors(Body.AF.x,Body.AF.z)[2:4]
        Body.sigma = nx*(Body.V0 + Body.vx) + nz*Body.vz   # normal vector pointing outward(overall sigma pointing outward)
        
        # Body doublet singularities influencing body itself
        # Transpose similar to phi_s
        Body.phi_db = np.transpose(-(np.arctan2(zp,xp2)\
                                  - np.arctan2(zp,xp1))/(2*np.pi))
        
        # Edge doublet influencing the body
        (xp1,xp2,zp) = transformation(Body.AF.x_col,Body.AF.z_col,Edge.x,Edge.z)
        
        if i > 1: # No wake panels until i==2
            # Wake doublets influencing the body
            (xp1_w,xp2_w,zp_w) = transformation(Body.AF.x_col,Body.AF.z_col,Wake.x[:i],Wake.z[:i])
            # Join edge and wake doublet influences into a single matrix
            xp1 = np.insert(xp1, Edge.N, xp1_w, axis=0)
            xp2 = np.insert(xp2, Edge.N, xp2_w, axis=0)
            zp = np.insert(zp, Edge.N, zp_w, axis=0)
    
        Body.phi_dw = np.transpose(-(np.arctan2(zp,xp2) \
                                   -np.arctan2(zp,xp1))/(2*np.pi))
        
def kutta(Body, Edge, Wake, RHO, SWITCH_KUTTA, DEL_T, i):
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
                
            
            Body.pressure(RHO, DEL_T, i)
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