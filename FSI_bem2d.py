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
from SolidClass import solid
from PyFEA import PyFEA
from FSIClass import FSI
from functions_general import archive, absoluteToBody, simulation_startup

# Turn on SIGFPE handling
np.seterr(all='raise')

DIO = DataIO(P)
start_time = time.time()

# Defining local variables to minimize dictionary lookups
DEL_T           = P['DEL_T']
DSTEP           = P['DSTEP']
TSTEP           = P['TSTEP']
T               = P['T']
RHO             = P['RHO']
RE              = P['RE']
VERBOSITY       = P['VERBOSITY']
COUPLING_SCHEME = P['COUPLING_SCHEME']
M_TYPE          = P['M_TYPE']
OUTER_CORR_TOL  = P['OUTER_CORR_TOL']
N_OUTERCORR_MAX = P['N_OUTERCORR_MAX']

(START_COUNTER, COUNTER, SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL) = simulation_startup(P, DIO, PC, Swimmer, solid, FSI, PyFEA)

po().calc_input(MotL[0].THETA_MAX/np.pi*180.,RE,MotL[0].THETA_MAX/np.pi*180.,DEL_T)
po().initialize_output(T[START_COUNTER])
outerCorr = 2

for i in xrange(START_COUNTER, COUNTER):
    if i == 0:
        for Swim in Swimmers:
            Swim.Body.free_swimming(P, i)
            Swim.Body.panel_positions(P, i)
            Swim.Body.surface_kinematics(P, i)
            Swim.edge_shed(DEL_T, i)
            Swim.wake_shed(DEL_T, i)
        solve_phi(Swimmers, P, i, outerCorr)
        for Swim in Swimmers:        
            Swim.Body.force(P, i)
            
        wake_rollup(Swimmers, DEL_T, i, P)
        archive(Swimmers[0].Body.AF.x_mid)
        archive(Swimmers[0].Body.AF.z_mid)
        SolidL[0].updateSolid(P['THETA'][i])
        graph.plot_n_go(Swimmers, i, P)
        DIO.write_data(P, i, DEL_T, SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL)
    else:
        if np.fmod(i, VERBOSITY) == 0:
            po().timestep_header(i,T[i])
            po().fsi_header()

        FSIL[0].readFsiControls(P)       
        FSIL[0].__init__(Swimmers[0].Body, SolidL[0])
        outerCorr = 0
        while True:
            outerCorr += 1
            FSIL[0].setInterfaceDisplacemet(outerCorr, COUPLING_SCHEME)
            for Swim in Swimmers:
                if (outerCorr == 1):
                    Swim.Body.free_swimming(P, i)
                    Swim.Body.panel_positions(P, i)
                else:
                    Swim.Body.fsi_panel_positions(FSIL[0], P, i)
                    
                Swim.Body.surface_kinematics(P, i)
                Swim.edge_shed(DEL_T, i)
                if (outerCorr == 1):
                    Swim.wake_shed(DEL_T, i)
                      
            solve_phi(Swimmers, P, i, outerCorr)
            
            for Swim in Swimmers:        
                Swim.Body.force(P, i)
                

            #TODO: Replace '0' with viscous drag component when available
            FSIL[0].setInterfaceForce(SolidL[0], Swimmers[0].Body, PyFEAL[0], 0., P, outerCorr, i)
            PyFEAL[0].dynamicSolve(Swimmers[0].Body, SolidL[0], outerCorr, M_TYPE)
            FSIL[0].getDisplacements(SolidL[0], Swimmers[0].Body, PyFEAL[0], P, i)
            FSIL[0].calcFSIResidual(SolidL[0], outerCorr)

            if np.fmod(i, VERBOSITY) == 0:
                po().fsi_iter_out(outerCorr,FSIL[0].fsiRelaxationFactor,FSIL[0].maxDU,FSIL[0].maxMagFsiResidual,FSIL[0].fsiResidualNorm,FSIL[0].maxFsiResidualNorm)

            if (FSIL[0].fsiResidualNorm <= OUTER_CORR_TOL or outerCorr >= N_OUTERCORR_MAX):
                if (FSIL[0].fsiResidualNorm <= OUTER_CORR_TOL):
                    po().fsi_converged()
                else:
                    po().fsi_not_converged()
                if np.fmod(i, VERBOSITY) == 0:
                    po().solution_output(Swimmers[0].Body.Cf, Swimmers[0].Body. Cl,Swimmers[0].Body.Ct,Swimmers[0].Body.Cpow)
                    po().solution_complete_output(i/float(COUNTER-1)*100.)
                wake_rollup(Swimmers, DEL_T, i, P)
                absoluteToBody(Swimmers[0].Body, SolidL[0], P, i)
                archive(Swimmers[0].Body.AF.x_mid)
                archive(Swimmers[0].Body.AF.z_mid)
                graph.plot_n_go(Swimmers, i, P)
                DIO.write_data(P, i, DEL_T, SwiL, GeoL, MotL, Swimmers, SolidL, FSIL, PyFEAL)
                break

total_time = time.time()-start_time
print "Simulation time:", np.round(total_time, 3), "seconds"