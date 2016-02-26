#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import time
import numpy as np
from data_IO_class import DataIO
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
if P['SW_FMM']:
    from functions_fmm import solve_phi, wake_rollup
else:
    from functions_influence import solve_phi, wake_rollup
from terminal_output import print_output as po
import functions_graphics as graph
from functions_general import archive, simulation_startup

#po().prog_title('1.0.0')
DIO = DataIO(P)
start_time = time.time()

# Defining local variables to minimize dictionary lookups
DEL_T     = P['DEL_T']
DSTEP     = P['DSTEP']
TSTEP     = P['TSTEP']
T         = P['T']
RHO       = P['RHO']
RE        = P['RE']
VERBOSITY = P['VERBOSITY']

(START_COUNTER, COUNTER, SwiL, GeoL, MotL, Swimmers) = simulation_startup(P, DIO, PC, Swimmer)[0:6]

po().calc_input(MotL[0].THETA_MAX/np.pi*180.,RE,MotL[0].THETA_MAX/np.pi*180.,DEL_T)
po().initialize_output(T[START_COUNTER])
outerCorr = 1

for i in xrange(START_COUNTER, COUNTER):
    if i == 0:
        for Swim in Swimmers:
                Swim.Body.panel_positions(P, i)
                Swim.Body.surface_kinematics(P, i)
                Swim.edge_shed(DEL_T, i)
                Swim.wake_shed(DEL_T, i)
        solve_phi(Swimmers, P, i, outerCorr)
        for Swim in Swimmers:
            Swim.Body.force(P, i)
            Swim.Body.free_swimming(P, i)
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.plot_n_go(Swimmers, i, P)
        DIO.write_data(P, i, DEL_T, SwiL, GeoL, MotL, Swimmers)

    else:
        if np.fmod(i, VERBOSITY) == 0:
            po().timestep_header(i,T[i])

        for Swim in Swimmers:
            Swim.Body.panel_positions(P, i)
            Swim.Body.surface_kinematics(P, i)
            Swim.edge_shed(DEL_T, i)
            Swim.wake_shed(DEL_T, i)
        solve_phi(Swimmers, P, i, outerCorr)
        wake_rollup(Swimmers, DEL_T, i, P)
        for Swim in Swimmers:
            Swim.Body.force(P, i)
            Swim.Body.free_swimming(P, i)
            if np.fmod(i, VERBOSITY) == 0:
                po().solution_output(Swim.Body.Cf, Swim.Body. Cl,Swim.Body.Ct,Swim.Body.Cpow)
                po().solution_complete_output(i/float(COUNTER-1)*100.)
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.plot_n_go(Swimmers, i, P)
        DIO.write_data(P, i, DEL_T, SwiL, GeoL, MotL, Swimmers)

total_time = time.time()-start_time
print "Simulation time:", np.round(total_time, 3), "seconds"