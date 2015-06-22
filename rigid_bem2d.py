#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""

import time
import numpy as np
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from functions_influence import quilt, wake_rollup
from terminal_output import print_output as po
import functions_graphics as graph
from functions_general import archive, geom_setup

po().prog_title('1.0.0')
start_time = time.time()

COUNTER = P['COUNTER']
DEL_T = P['DEL_T']
DSTEP = P['DSTEP']
TSTEP = P['TSTEP']
T = [DEL_T*i for i in xrange(COUNTER)]
RHO = P['RHO']
RE = P['RE']

#SPLIT = 0.4
#SwiP = PC.SwimmerParameters(P['CE'], P['DELTA_CORE'], P['SW_KUTTA'])
##    GeoP = PC.GeoVDVParameters(P['N_BODY'], P['S'], P['C'], P['K'], P['EPSILON'])
#GeoP = PC.GeoFPParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
##    GeoP = PC.GeoTDParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
#MotP1 = PC.MotionParameters(0., 0*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
#MotP2 = PC.MotionParameters(0., 1*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
#MotP3 = PC.MotionParameters(0., 2*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
#MotP4 = PC.MotionParameters(0., 3*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
#MotP5 = PC.MotionParameters(0., 4*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
#
#S1 = Swimmer(SwiP, GeoP, MotP1, COUNTER-1)
#S2 = Swimmer(SwiP, GeoP, MotP2, COUNTER-1)
#S3 = Swimmer(SwiP, GeoP, MotP3, COUNTER-1)
#S4 = Swimmer(SwiP, GeoP, MotP4, COUNTER-1)
#S5 = Swimmer(SwiP, GeoP, MotP5, COUNTER-1)
#Swimmers = [S1]

(SwiP, GeoP, MotP, Swimmers) = geom_setup(P, PC, Swimmer)[0:4]

po().calc_input(MotP[0].THETA_MAX/np.pi*180.,RE,MotP[0].THETA_MAX/np.pi*180.,DEL_T)

# Data points per cycle == 1/(F*DEL_T)
for i in xrange(COUNTER):
    if i == 0:
        po().initialize_output(T[i])

        for Swim in Swimmers:
                Swim.Body.panel_positions(DSTEP, T[i], P['THETA'][i])
                Swim.Body.surface_kinematics(DSTEP, TSTEP, P['THETA_MINUS'][i], P['THETA_PLUS'][i], DEL_T, T[i], i)
                Swim.edge_shed(DEL_T, i)
                Swim.wake_shed(DEL_T, i)
        quilt(Swimmers, RHO, DEL_T, i)
        wake_rollup(Swimmers, DEL_T, i)
#        archive(S1.Body.AF.x_mid)
#        archive(S1.Body.AF.z_mid)
#        graph.body_plot(S1.Edge, S1.Body)  
        for Swim in Swimmers:
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.body_plot(Swimmers[0].Edge, Swimmers[0].Body)  
    else:
        if np.fmod(i,P['VERBOSITY']) == 0:
            po().timestep_header(i,T[i])

        for Swim in Swimmers:
            Swim.Body.panel_positions(DSTEP, T[i], P['THETA'][i])
            Swim.Body.surface_kinematics(DSTEP, TSTEP, P['THETA_MINUS'][i], P['THETA_PLUS'][i], DEL_T, T[i], i)
            Swim.edge_shed(DEL_T, i)
            Swim.wake_shed(DEL_T, i)
        quilt(Swimmers, RHO, DEL_T, i)

        if np.fmod(i,P['VERBOSITY']) == 0:
            po().solution_output(0,0,0,0,0,0)
            po().solution_complete_output(i/float(COUNTER-1)*100.)
        wake_rollup(Swimmers, DEL_T, i)
#        archive(S1.Body.AF.x_mid)
#        archive(S1.Body.AF.z_mid)
#        graph.body_plot(S1.Edge, S1.Body)
#        archive(Swimmers[0].Body.AF.x_mid)
#        archive(Swimmers[0].Body.AF.z_mid)
        for Swim in Swimmers:
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.body_plot(Swimmers[0].Edge, Swimmers[0].Body) 


total_time = time.time()-start_time
print "Simulation time:", np.round(total_time, 3), "seconds"

graph.body_wake_plot(Swimmers)
#    graph.cp_plot(S1.Body)