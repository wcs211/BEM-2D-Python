#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
from functions_influence import influence_matrices


def solve_phi(Swimmers, RHO, DEL_T, i, outerCorr):
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
    
    (sigma_all, mu_w_all, a_b, a_e, b_b, b_e, b_w) = influence_matrices(Swimmers, i)
        
    if i == 0:
        b = -np.dot(b_b, sigma_all)
    else:
        b = -np.dot(b_b, sigma_all) - np.dot(b_w, mu_w_all)
            
    places = np.zeros((2,Swimmers[0].Body.N))
    places[0,:] = Swimmers[0].Body.AF.x_mid[0,:]
    places[1,:] = Swimmers[0].Body.AF.z_mid[0,:]
    
    (mu_b_all, info) = mygmres(places.T, b)
    
    for Swim in Swimmers:
        Swim.Body.pressure(RHO, DEL_T, i) 
    
    for Swim in Swimmers:
        Swim.mu_guess = np.empty(2) # [0] is current guess, [1] is previous
        Swim.delta_p = np.empty(2) # [0] is current delta_p, [1] is previous
        Swim.Body.mu[:] = mu_b_all[Swim.i_b:Swim.i_b+Swim.Body.N]
        Swim.mu_guess[0] = Swim.Body.mu[-1]-Swim.Body.mu[0]
    
    for Swim in Swimmers:  
        Swim.Edge.mu = Swim.mu_guess[0]
        Swim.Edge.gamma[0] = -Swim.Edge.mu
        Swim.Edge.gamma[1] = Swim.Edge.mu

        # Get gamma of body panels for use in wake rollup
        Swim.Body.gamma[0] = -Swim.Body.mu[0]
        Swim.Body.gamma[1:-1] = Swim.Body.mu[:-1]-Swim.Body.mu[1:]
        Swim.Body.gamma[-1] = Swim.Body.mu[-1]
