#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import os
import numpy as np
import pickle
from data_IO_class import DataIO
#from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from functions_influence import solve_phi, wake_rollup
from terminal_output import print_output as po
import functions_graphics as graph
from SolidClass import solid
from PyFEA import PyFEA
from FSIClass import FSI
from functions_general import archive, absoluteToBody, simulation_startup

# Input Directory
INPUT_DIR = '/home/wcs211/Python_Post-processing/data/'

i = len([name for name in os.listdir(INPUT_DIR) if os.path.isfile(os.path.join(INPUT_DIR, name))])
IF = np.empty(i)
j = 0
for file in os.listdir(INPUT_DIR):
    IF[j] = float(file)
    j = j + 1

INPUT_FILES = sorted(IF)
csvData = np.zeros((np.shape(INPUT_FILES)[0], 6))

for i_t in xrange(np.shape(INPUT_FILES)[0]):
    INPT_FL = "%s%.8f" % (INPUT_DIR, INPUT_FILES[i_t])
    
    with open(INPT_FL, 'rb') as f:
        P, i, FLOWTIME, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA = pickle.load(f)
    
    csvData[i_t,0] = i
    csvData[i_t,1] = P['T'][i]
    csvData[i_t,2] = Swimmers[0].Body.Cf
    csvData[i_t,3] = Swimmers[0].Body.Cl
    csvData[i_t,4] = Swimmers[0].Body.Ct
    csvData[i_t,5] = Swimmers[0].Body.Cpow
    
np.savetxt("performance.csv", csvData, delimiter=",", header="i [-], Flow Time [s], Cf [-], Cl [-], Ct [-], Cpow [-]")