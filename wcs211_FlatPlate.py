#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-3D
A 3D boundary element method code

"""
import numpy as np
import matplotlib.pyplot as plt
from TDinput import PARAMETERS as IP
from mpl_toolkits.mplot3d import Axes3D 

# Flat plate geometry
# Stepping through each spanwise position to calculate the positions of the
# fluke neutral plane at the given time step.

# Geometry parameters
Nx = IP['Nx']
Ny = IP['Ny']
h = IP['h']
w = IP['w']
D = 0.2*h

# Defining mid-chord line
MC = 0.5 * h * np.ones(Ny)

# Defining Chord length function
C = h*np.ones(Ny)

# Defining Leading and Trailing Edge Line

start = 0
stop  = np.pi
step  = np.copy(Ny)


   
TE = MC[:,0] + C/2
LE = MC[:,0] - C/2

theta = np.linspace(start,stop,step)
xb = (C*np.cos(theta).T + C)/(2.)

#print xb
start = np.pi
stop  = 0
step  = np.copy(Ny)
theta = np.linspace(start,stop,step)
xt = (C*np.cos(theta).T + C)/(2.)