import time
import numpy as np
from terminal_output import print_output as po
from setup_parameters import PARAMETERS as P
from body_class import Body, GeoParametersVDV, MotionParameters, SwimmerParameters
from edge_class import Edge
from wake_class import Wake
from general_functions import archive
from kinematics import neutral_axis, panel_positions, surface_kinematics, edge_shed, wake_shed, wake_rollup
from inf_force import influence_matrices, kutta
import graphics as graph

def main():

    po().prog_title('1.0.0')
    start_time = time.time()
    
    # Recurring simulation constants
    COUNTER = P['COUNTER']
    DEL_T = P['DEL_T']
    DSTEP = P['DSTEP']
    TSTEP = P['TSTEP']
    t = [DEL_T*i for i in xrange(COUNTER)]
    
    PGeo = GeoParametersVDV(P['N_BODY'], P['C'], P['K'], P['EPSILON'])
    PMotion = MotionParameters(P['V0'], P['THETA_MAX'], P['H_C'], P['F'], P['PHI'])
    PSwimmer = SwimmerParameters(P['CE'], P['S'], P['SW_GEOMETRY'], P['SW_KUTTA'])
    
    Body1 = Body.from_van_de_vooren(PGeo, PMotion)
    Edge1 = Edge(P['V0'],P['CE'],COUNTER)
    Wake1 = Wake(P['V0'],COUNTER-2)
#    FSI1 = FSI(Npanels,Nelements)
#    PyFEA1 = PyFEA(Nelements, fracDeltaT, endTime, E, I, A, l, rho, Fload, U_n, Udot_n)
#    Solid1 = solid(Nnodes,xp_0,zp_0,tmax)
    
    po().calc_input(P['THETA_MAX']/np.pi*180.,P['RE'],P['THETA_MAX']/np.pi*180.,DEL_T)
    
    # Data points per cycle == 1/(F*DEL_T)
    for i in xrange(COUNTER):
        if i == 0:
            po().initialize_output(t[i])
    
        else: #i > 0:
        
#            readFsiControls(fixedPtRelax, nOuterCorrMax)
#            FSI1.__init__(Npanels,Nelements)
            
            while True:
                
#                FSI1.setInterfaceDisplacemet(displ, relaxationFactor, 
#                                             residual, outerCorr, couplingScheme)
            
                (Body1.x_neut,Body1.z_neut) = neutral_axis(Body1.PMotion,Body1.BFC.x_col[:Body1.N/2],DSTEP,TSTEP,t[i])
                archive(Body1.x_mid)
                archive(Body1.z_mid)
                (Body1.AFC, Body1.x_mid[0,:], Body1.z_mid[0,:]) = panel_positions(Body1.BFC,Body1.PMotion,P['S'],DSTEP,t[i])
                surface_kinematics(Body1,DSTEP,TSTEP,t[i],i,DEL_T)
                edge_shed(Body1,Edge1,i,DEL_T)
                wake_shed(Edge1,Wake1,i,DEL_T)
                
                influence_matrices(Body1,Edge1,Wake1,i)
                kutta(Body1,Edge1,Wake1,P['RHO'],i,DEL_T,P['SW_KUTTA'])
                
#                FSI1.setInterfaceForce(outerCorr, nodes, nodesNew, theta, heave, 
#                                       x_b, z_b, xp, zp, xc, zc, P_b, ViscDrag, vn, delFs, 
#                                       interpMtd, meanline_c0, tBeamStruct, fixedCounter, c,
#                                       U_nPlus, Udot_nPlus, i_t)
#                PyFEA1.solve(mType, method, theta, fixedNodes, alpha, beta, gamma)
#                (DU, tempNodes) = FSI1.getDisplacements(theta, heavePos, xp, zp, nodalDelxp, nodalDelzp, tBeam, nodes_0, U_nPlus, interpMtd, meanline_p0, fixedNodes, flexionRatio)
#                FSI1.calcFSIResidual(DU, nodes, tempNodes, outerCorr)
#                
#                if (FSI1.fsiResidualNorm <= FSI1.outerCorrTolerance or FSI1.outerCorr >= FSI1.nOuterCorr):
                wake_rollup(Body1,Edge1,Wake1,P['DELTA_CORE'],i,DEL_T)
                
                #force(Body1,i)
                #po().solution_output(d_visc,cf,cl,ct,cpow,gamma)
                 
                if np.fmod(i,10) == 0:
                    po().timestep_header(i+1,t[i])
                    po().solution_output(0,0,0,0,0,0)
                    po().solution_complete_output(i/float(COUNTER-1)*100.)
                
                break # TODO: remove this once FSI is ready
    
    total_time = time.time()-start_time
    print "Simulation time:", np.round(total_time, 3), "seconds"

    graph.body_wake_plot(Body1,Edge1,Wake1)
    #graph.cp_plot(Body1)

if __name__ == '__main__':
    main()