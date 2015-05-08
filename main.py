import time
import numpy as np
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from functions_influence import quilt, wake_rollup
from terminal_output import print_output as po
import functions_graphics as graph

def main():

    po().prog_title('1.0.0')
    start_time = time.time()
    
    COUNTER = P['COUNTER']
    DEL_T = P['DEL_T']
    DSTEP = P['DSTEP']
    TSTEP = P['TSTEP']
    T = [DEL_T*i for i in xrange(COUNTER)]
    RHO = P['RHO']
    RE = P['RE']

    SwiP = PC.SwimmerParameters(P['CE'], P['DELTA_CORE'], P['SW_KUTTA'])
    GeoP = PC.GeoVDVParameters(P['N_BODY'], P['S'], P['C'], P['K'], P['EPSILON'])
    MotP1 = PC.MotionParameters(0., 0.3, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
    MotP2 = PC.MotionParameters(0., -0.3, P['V0'], P['THETA_MAX'], P['F'], P['PHI']+np.pi)
    
    Swimmer1 = Swimmer(SwiP, GeoP, MotP1, COUNTER-2)
    Swimmer2 = Swimmer(SwiP, GeoP, MotP2, COUNTER-2)
    Swimmers = [Swimmer1, Swimmer2]
    
#    FSI1 = FSI(Npanels,Nelements)
#    PyFEA1 = PyFEA(Nelements, fracDeltaT, endTime, E, I, A, l, rho, Fload, U_n, Udot_n)
#    Solid1 = solid(Nnodes,xp_0,zp_0,tmax)
    
    po().calc_input(MotP1.THETA_MAX/np.pi*180.,RE,MotP1.THETA_MAX/np.pi*180.,DEL_T)
    
    # Data points per cycle == 1/(F*DEL_T)
    for i in xrange(COUNTER):
        if i == 0:
            po().initialize_output(T[i])
    
        else: #i > 0:
        
#            readFsiControls(fixedPtRelax, nOuterCorrMax)
#            FSI1.__init__(Npanels,Nelements)
            
            while True:
                
#                FSI1.setInterfaceDisplacemet(displ, relaxationFactor, 
#                                             residual, outerCorr, couplingScheme)
                for Swim in Swimmers:
                    Swim.Body.panel_positions(DSTEP, T[i])
                    Swim.Body.surface_kinematics(DSTEP, TSTEP, DEL_T, T[i], i)
                    Swim.edge_shed(DEL_T, i)
                    Swim.wake_shed(DEL_T, i)
                    
                quilt(Swimmers, RHO, DEL_T, i)
                
#                FSI1.setInterfaceForce(outerCorr, nodes, nodesNew, theta, heave, 
#                                       x_b, z_b, xp, zp, xc, zc, P_b, ViscDrag, vn, delFs, 
#                                       interpMtd, meanline_c0, tBeamStruct, fixedCounter, c,
#                                       U_nPlus, Udot_nPlus, i_t)
#                PyFEA1.solve(mType, method, theta, fixedNodes, alpha, beta, gamma)
#                (DU, tempNodes) = FSI1.getDisplacements(theta, heavePos, xp, zp, nodalDelxp, nodalDelzp, tBeam, nodes_0, U_nPlus, interpMtd, meanline_p0, fixedNodes, flexionRatio)
#                FSI1.calcFSIResidual(DU, nodes, tempNodes, outerCorr)
#                
#                if (FSI1.fsiResidualNorm <= FSI1.outerCorrTolerance or FSI1.outerCorr >= FSI1.nOuterCorr):
                wake_rollup(Swimmers, DEL_T, i)
                
                #force(Body1,i)
                #po().solution_output(d_visc,cf,cl,ct,cpow,gamma)
                 
                if np.fmod(i,10) == 0:
                    po().timestep_header(i+1, T[i])
                    po().solution_output(0,0,0,0,0,0)
                    po().solution_complete_output(i/float(COUNTER-1)*100.)
                
                break # TODO: remove this once FSI is ready
    
    total_time = time.time()-start_time
    print "Simulation time:", np.round(total_time, 3), "seconds"

    graph.body_wake_plot(Swimmers)
#    graph.cp_plot(Swimmer1.Body)

if __name__ == '__main__':
    main()