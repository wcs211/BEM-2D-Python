#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
from input_parameters import PARAMETERS as P

importCsv = np.genfromtxt('forces.csv', delimiter=',')
dataLength = importCsv.shape[0]
avgData = np.zeros((P['N_CYC'], 6))
norm = P['HEAVE_MAX']

for i in xrange(P['N_CYC']):
    avgData[i,0] = i+1
    avgData[i,1] = np.trapz(importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,1], x=importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,0]*P['DEL_T']) / (P['T'][(i+1)*P['N_STEP']] - P['T'][i*P['N_STEP']])
    avgData[i,2] = np.trapz(importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,2], x=importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,0]*P['DEL_T']) / (P['T'][(i+1)*P['N_STEP']] - P['T'][i*P['N_STEP']])
    avgData[i,3] = np.trapz(importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,3], x=importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,0]*P['DEL_T']) / (P['T'][(i+1)*P['N_STEP']] - P['T'][i*P['N_STEP']])
    avgData[i,4] = np.trapz(importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,4], x=importCsv[i*P['N_STEP']:(i+1)*P['N_STEP']+1,0]*P['DEL_T']) / (P['T'][(i+1)*P['N_STEP']] - P['T'][i*P['N_STEP']])

    for i_t in xrange(dataLength):
        if (i_t >= i * P['N_STEP'] and i_t <= (i+1) * P['N_STEP']):
            if (np.abs(importCsv[i_t, 5]) > avgData[i,5]):
                avgData[i,5] = np.abs(importCsv[i_t, 5])

np.savetxt("avg_Performance.csv", avgData, delimiter=",", header="Cycle No. [-], Cf_avg [-], Cl_avg [-], Ct_avg [-], Cpow_avg [-], TE_A_max [m]")