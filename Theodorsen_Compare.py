#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
from input_parameters import PARAMETERS as P
from scipy.special import j0, j1, y0, y1
import matplotlib.pyplot as plt

tau = np.linspace(0, 5, P['N_STEP']*P['N_CYC']+1)

b = 0.5 * P['C']
a = -1.0
w = 2. * np.pi * P['F']
t = np.asarray(P['T'])
RHO = P['RHO']
U = np.absolute(P['V0'])
THETA_MAX = P['THETA_MAX']
HEAVE_MAX = P ['HEAVE_MAX']
PHI = P['PHI']

k = w * b / np.absolute(P['V0'])
k2 = np.pi * P['F'] * P['C'] / np.absolute(P['V0'])
St = 2. * P['F'] * P['HEAVE_MAX'] / np.absolute(P['V0'])

F = (j1(k)*(j1(k)+y0(k)) + y1(k)*(y1(k)-j0(k))) / ((j1(k)+y0(k))**2 + (y1(k)-j0(k))**2)
G = -(y1(k)*y0(k) + j1(k)*j0(k)) / ((j1(k)+y0(k))**2 + (y1(k)-j0(k))**2)

L = -RHO * b**2 * (U * np.pi * THETA_MAX * w * np.cos(w * t + PHI) - np.pi * HEAVE_MAX * w**2 * np.sin(w * t) + np.pi * b * a * THETA_MAX * w**2 * np.sin(w * t + PHI)) - 2. * np.pi * RHO * U * b * F * (U * THETA_MAX * np.sin(w * t + PHI) + HEAVE_MAX * w * np.cos(w * t) + b * (0.5 - a) * THETA_MAX * w * np.cos(w * t + PHI)) - 2. * np.pi * RHO * U * b * G * (U * THETA_MAX * np.cos(w * t + PHI) - HEAVE_MAX * w * np.sin(w * t) - b * (0.5 - a) * THETA_MAX * w * np.sin(w * t + PHI))

Cl = np.real(L) / (0.5 * P['RHO'] * np.absolute(P['V0'])**2 * P['C'])

execfile("rigid_bem2d.py")
expCsv = np.genfromtxt('forces.csv', delimiter=',')

#P['SW_KUTTA'] = True
#
#execfile("rigid_bem2d_kutta.py")
#impCsv = np.genfromtxt('forces.csv', delimiter=',')

#Determine Error
indx = (P['N_CYC'] - 1) * P['N_STEP']
diff = np.absolute((expCsv[indx:,2] - Cl[indx:]))
print 'Maximum Error = ', np.max(diff)

# Plot the results
figure = plt.figure(1)
figure.add_subplot(1, 1, 1, axisbg='1') # Change background color here
plt.tick_params(labelsize=28)
plt.plot(np.asanyarray(P['T'])/(1. / P['F']) - (P['N_CYC']-1), expCsv[:,2],'-o')
plt.plot(np.asanyarray(P['T'])/(1. / P['F']) - (P['N_CYC']-1), Cl)

#plt.plot(np.asanyarray(P['T']), impCsv[:,2],'x')
plt.legend(('BEM','Theodorsen'), fontsize=28)
#plt.legend(('Theodorsen','Explicit Kutta','Implicit Kutta'))
plt.xlabel(r'$\tau = t/T $', fontsize=28)
plt.ylabel(r'$ C_l $', fontsize=28)
xmin = (P['N_CYC'] - 1) * P['N_STEP'] * P['DEL_T'] /(1. / P['F'])  - (P['N_CYC']-1)
xmax = P['N_CYC'] * P['N_STEP'] * P['DEL_T'] /(1. / P['F']) - (P['N_CYC']-1)
ymin = np.min(Cl[indx:]) + 0.125*np.min(Cl[indx:])
ymax = -np.min(Cl[indx:]) - 0.125*np.min(Cl[indx:])
plt.axis([xmin, xmax, ymin, ymax])
plt.show()

