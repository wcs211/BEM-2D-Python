#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import os
import pickle
import numpy as np

class DataIO(object):
    def __init__(self, P):
        """
        Handles data input and output for the simualtion.
        
        Args:
            P (dict): Input Parameters
        """
        self.OUTPUT_DIR = P['OUTPUT_DIR']
        
        # Determine if the output directory exists. If not, create the directory.
        if not os.path.exists(self.OUTPUT_DIR):
            os.makedirs(self.OUTPUT_DIR)
            
    def write_data(self, P, i, DEL_T, SwiP, GeoP, MotP, Swimmers, solid=None, FSI=None, PyFEA=None):
        """
        Writes the simulation state to a binary file named with the current flow 
        time. Also writes a csv containing force information for each time-step.
        
        Args:
            P (dict): Input Parameters
            i (int): Time-step Counter
            DEL_T (float): Time-step
            SwiP (list): Swiimer Parameters
            GeoP (list): Geometry Parameters
            MotP (list): Motion Parameters
            Swimmers (list): Swimmer body, edge panel, and wake
            solid (Optional [list]): Solid Parameters
            FSI (Optional [list]): Fluid-structure Interaction Parameters
            PyFEA (Optional [list]): Finite Element Instance
        """
        # Determine if a csv file containg force information should be written
        if P['SW_SV_FORCES']:
            Swimmers[0].Body.forceData = np.append(Swimmers[0].Body.forceData, np.array([[i, Swimmers[0].Body.Cf, Swimmers[0].Body.Cl, Swimmers[0].Body.Ct, Swimmers[0].Body.Cpow, Swimmers[0].Body.AF.z[0], Swimmers[0].Body.V/Swimmers[0].Body.V0]]), axis=0)
            np.savetxt("forces.csv", Swimmers[0].Body.forceData, delimiter=",", header="i [-], Cf [-], Cl [-], Ct [-], Cpow [-], TE_A [-], V/V0 [-]")
#            np.savetxt("forces.csv", Swimmers[0].Body.forceData, delimiter=",")   
        
        # Determine if a save-state file should be written
        # Check if the save-state should be written for the last cycle only
        if (i >= (P['N_CYC']-1) * P['N_STEP'] and P['SW_SV_L_CYCLE'] == True):
            # Check if a save-state should be written to disk
            if (np.fmod(i,P['SAVE_EVERY']) == 0 and P['SW_SAVE_DATA'] == True):
                # Determine the output file name
                outfile = ''.join((self.OUTPUT_DIR, "/%.8f" % (i * DEL_T)))
                # Write the save-state to disk
                with open(outfile, 'wb') as f:
                    pickle.dump([P, i, i*DEL_T, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA], f, 1)
        else:
            # Check if a save-state should be written to disk
            if (np.fmod(i,P['SAVE_EVERY']) == 0 and P['SW_SAVE_DATA'] == True and P['SW_SV_L_CYCLE'] == False):
                # Determine the output file name
                outfile = ''.join((self.OUTPUT_DIR, "/%.8f" % (i * DEL_T)))
                # Write the save-state to disk
                with open(outfile, 'wb') as f:
                    pickle.dump([P, i, i*DEL_T, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA], f, 1)
            
    def read_data(self, INPUT_FILE):
        """
        Reads a simulation state file.
        
        Args:
            INPUT_FILE (str): File name to read
        """
        # Read the save-state file
        with open(INPUT_FILE, 'rb') as f:
            P, i, FLOWTIME, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA = pickle.load(f)
            
        return (P, i, FLOWTIME, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA)