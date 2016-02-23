#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
from functions_general import panel_vectors, transformation
from functions_influence import inf_doubletpanel, quilt
from pyfmmlib import fmm_part
import matplotlib.pyplot as plt
    
def influence_matrices(Swimmers, i):
    """Constructs the influence coefficient matrices.

    Args:
        Swimmers: List of Swimmer objects being simulated.
        i: Time step number.

    Returns:
        sigma_all: Array containing all Swimmers' body source strengths.
        a_bodydoublet: Body doublet panels' influence matrix.
        b_bodysource: Body source panels' influence matrix.
        b_edgedoublet: Edge panels' influence matrix.
        b_wakedoublet: Wake panels' influence matrix.
        a_explicit: Augments the a_bodydoublet matrix when doing explicit Kutta.
    """
    ep = 2.2204460492503131e-16
    n_b = 0
    n_e = 0
    n_w = 0
    for Swim in Swimmers:
        Swim.i_b = n_b
        Swim.i_e = n_e
        Swim.i_w = n_w
        n_b += Swim.Body.N
        n_e += 1
        n_w += i

    sigma_all = np.empty(n_b)
    mu_w_all = np.empty(n_w)
    for Swim in Swimmers:
        (r0, rn) = (Swim.i_b, Swim.i_b+Swim.Body.N)
        sigma_all[r0:rn] = Swim.Body.sigma[:]
        (r0, rn) = (Swim.i_w, Swim.i_w+i)
        mu_w_all[r0:rn] = Swim.Wake.mu[:i]

    (xp1, xp2, zp) = quilt(Swimmers, 'Body', n_b, n_b, i)
    # Calculate 'Mask' - if a source or doublet gets too close then ignore its influence
    mask = np.greater_equal(np.absolute(zp),ep).astype(int)
    
    # Body doublet singularities influencing bodies themselves (the A matrix)
    a_bodydoublet = inf_doubletpanel(xp1, xp2, zp, mask)

    (xp1, xp2, zp) = quilt(Swimmers, 'Edge', n_b, n_e, i)
    # Calculate 'Mask' - if a source or doublet gets too close then ignore its influence
    mask = np.greater_equal(np.absolute(zp),ep).astype(int)
    
    # Edge doublet singularities influencing the bodies (part of RHS)
    b_edgedoublet = inf_doubletpanel(xp1, xp2, zp, mask)

    a_explicit = np.zeros((n_b,n_b))

    return(sigma_all, mu_w_all, a_bodydoublet, a_explicit, b_edgedoublet)

def solve_phi(Swimmers, RHO, DEL_T, i, outerCorr=0):
    """Solves the boundary integral equation using a Kutta condition and the
    Fast Multipole Method.

    Args:
        Swimmers: List of Swimmer objects being simulated.
        RHO: Fluid density.
        DEL_T: Time step length.
        i: Time step number.
    """
    for Swim in Swimmers:  
        if (outerCorr <= 1):
            # mu_past used in differencing for pressure
            Swim.Body.mu_past[1:4,:] = Swim.Body.mu_past[0:3,:]
            Swim.Body.mu_past[0,:] = Swim.Body.mu

    (sigma_all, mu_w_all, a_b, a_e, b_e) = influence_matrices(Swimmers, i)
    
    for SwimI in Swimmers:
        a_e[:,SwimI.i_b] = -b_e[:, SwimI.i_e]
        a_e[:,SwimI.i_b+SwimI.Body.N-1] = b_e[:, SwimI.i_e]
    a = a_b + a_e

    # Prepare FMM inputs
    b_places = np.vstack((Swimmers[0].Body.AF.x_mid[0,:], Swimmers[0].Body.AF.z_mid[0,:]))
    target = np.vstack((Swimmers[0].Body.AF.x_col, Swimmers[0].Body.AF.z_col))
    lpanel = panel_vectors(Swimmers[0].Body.AF.x, Swimmers[0].Body.AF.z)[-1]
    if (i > 0):
        w_places = np.vstack((0.5*(Swimmers[0].Wake.x[1:i+1] + Swimmers[0].Wake.x[:i]), 0.5*(Swimmers[0].Wake.z[1:i+1] + Swimmers[0].Wake.z[:i])))
        (nx, nz, lwpanel) = panel_vectors(Swimmers[0].Wake.x[:i+1], Swimmers[0].Wake.z[:i+1])[2:]
        n_w = np.vstack((nx,nz))

    # Get right-hand side
    b_body = np.real(fmm_part("P", iprec=5, kernel=0, sources=b_places.T, target=target.T, mop_charge=sigma_all * lpanel)) / 2. / np.pi 
    if i == 0:
        b = -b_body
    else:
        b_wake = np.real(fmm_part("P", iprec=5, kernel=0, sources=w_places.T, target=target.T, dipvec=n_w.T , dip_charge=mu_w_all * lwpanel)) / 2. / np.pi
        b = -b_body - b_wake
    
    # Solve for bodies' doublet strengths using explicit Kutta
    mu_b_all = np.linalg.solve(a, b)
    
    for Swim in Swimmers:
        Swim.Body.pressure(RHO, DEL_T, i)     

        Swim.mu_guess = np.empty(2) # [0] is current guess, [1] is previous
        Swim.delta_p = np.empty(2) # [0] is current delta_p, [1] is previous
        Swim.Body.mu[:] = mu_b_all[Swim.i_b:Swim.i_b+Swim.Body.N]
        Swim.mu_guess[0] = Swim.Body.mu[-1]-Swim.Body.mu[0]
  
        Swim.Edge.mu = Swim.mu_guess[0]
        Swim.Edge.gamma[0] = -Swim.Edge.mu
        Swim.Edge.gamma[1] = Swim.Edge.mu

        # Get gamma of body panels for use in wake rollup
        Swim.Body.gamma[0] = -Swim.Body.mu[0]
        Swim.Body.gamma[1:-1] = Swim.Body.mu[:-1]-Swim.Body.mu[1:]
        Swim.Body.gamma[-1] = Swim.Body.mu[-1]

def wake_rollup(Swimmers, DEL_T, i, P):
    """Performs wake rollup on the swimmers' wake panels.

    Args:
        Swimmers: List of Swimmer objects being simulated.
        DEL_T: Time step length.
        i: Time step number.
    """
    if (P['SW_ROLLUP']):
        # Wake panels initialize when i==1
        if i == 0:
            pass
    
        else:
            NT = i # Number of targets (wake panel points that are rolling up)
            for SwimT in Swimmers:
                SwimT.Wake.vx = np.zeros(NT)
                SwimT.Wake.vz = np.zeros(NT)
                
                wake_x_midT = 0.5*(SwimT.Wake.x[1:i+1] + SwimT.Wake.x[:i])
                wake_z_midT = 0.5*(SwimT.Wake.z[1:i+1] + SwimT.Wake.z[:i])

                target = np.vstack((wake_x_midT, wake_z_midT))   
                
                for SwimI in Swimmers:
                    wake_x_midI = 0.5*(SwimI.Wake.x[1:i+1] + SwimI.Wake.x[:i])
                    wake_z_midI = 0.5*(SwimI.Wake.z[1:i+1] + SwimI.Wake.z[:i])
                    edge_x_mid  = 0.5*(SwimI.Edge.x[1]     + SwimI.Edge.x[0] )
                    edge_z_mid  = 0.5*(SwimI.Edge.z[1]     + SwimI.Edge.z[0] )
                    
                    ed_places = np.vstack((edge_x_mid              , edge_z_mid              ))
                    bd_places = np.vstack((SwimT.Body.AF.x_col     , SwimT.Body.AF.z_col     ))
                    bs_places = np.vstack((SwimT.Body.AF.x_mid[0,:], SwimT.Body.AF.z_mid[0,:]))
                    wd_places = np.vstack((wake_x_midI             , wake_z_midI             ))   
                    
                    (nx_b, nz_b, lbpanel) = panel_vectors(SwimI.Body.AF.x   , SwimI.Body.AF.z   )[2:]
                    (nx_e, nz_e, lepanel) = panel_vectors(SwimI.Edge.x      , SwimI.Edge.z      )[2:]
                    (nx_w, nz_w, lwpanel) = panel_vectors(SwimI.Wake.x[:i+1], SwimI.Wake.z[:i+1])[2:]
                    
                    n_b = np.vstack((nx_b,nz_b))
                    n_e = np.vstack((nx_e,nz_e))
                    n_w = np.vstack((nx_w,nz_w))
                    
                    # Body source influence on wake velocity
                    u_bs = np.real(fmm_part("G", iprec=5, kernel=0, sources=bs_places.T, target=target.T, mop_charge=SwimI.Body.sigma * lbpanel)) / np.pi
                    SwimT.Wake.vx +=  2. * u_bs[:,0]
                    SwimT.Wake.vz += -2. * u_bs[:,1]
                    
                    # Body doublet influence on wake velocity
                    u_bd = np.real(fmm_part("G", iprec=5, kernel=0, sources=bd_places.T, target=target.T, dipvec=n_b.T , dip_charge=SwimI.Body.mu * lbpanel)) / np.pi
                    SwimT.Wake.vx +=  2. * u_bd[:,0]
                    SwimT.Wake.vz += -2. * u_bd[:,1]
                    
                    # TE panel influence on wake velocity
                    u_te = np.real(fmm_part("G", iprec=5, kernel=0, sources=ed_places.T, target=target.T, dipvec=n_e.T , dip_charge=SwimI.Edge.mu * lepanel)) / np.pi
                    SwimT.Wake.vx +=  2. * u_te[:,0]
                    SwimT.Wake.vz += -2. * u_te[:,1]
                    
                    # Wake double influence on wake velocity
                    u_w = np.real(fmm_part("G", iprec=5, kernel=0, sources=wd_places.T, target=target.T, dipvec=n_w.T , dip_charge=SwimI.Wake.mu[1:i+1] * lwpanel)) / np.pi
                    u_w = np.nan_to_num(u_w)
                    SwimT.Wake.vx +=  2. * u_w[:,0]
                    SwimT.Wake.vz += -2. * u_w[:,1]
    
            for Swim in Swimmers:
                # Modify wake with the total induced velocity
                Swim.Wake.x[1:i+1] += Swim.Wake.vx*DEL_T
                Swim.Wake.z[1:i+1] += Swim.Wake.vz*DEL_T
