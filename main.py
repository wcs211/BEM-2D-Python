import time
import numpy as np
from terminal_output import print_output as po
from setup_parameters import PARAMETERS as P
from body_class import Body
from edge_class import Edge
from wake_class import Wake
from kinematics import neutral_axis, panel_positions, surface_kinematics, edge_shed, wake_shed, wake_rollup
from inf_force import influence_matrices, kutta
import graphics as graph

def main():

    po().prog_title('1.0.0')
    start_time = time.time()
    
    # Recurring parameters
    COUNTER = P['COUNTER']
    DEL_T = P['DEL_T']
    DSTEP = P['DSTEP']
    TSTEP = P['TSTEP']
    
    Body1 = Body.from_van_de_vooren([P[key] for key in ['N_BODY','C','K','EPSILON','V0','THETA_MAX','H_C','F','PHI']])
    Edge1 = Edge(Body1.V0,P['CE'],COUNTER)
    Wake1 = Wake(Body1,COUNTER-2)

    t = DEL_T*np.arange(0,COUNTER)
    
    po().calc_input(P['THETA_MAX']/np.pi*180.,P['RE'],P['THETA_MAX']/np.pi*180.,DEL_T)
    
    # Data points per cycle == 1/(F*DEL_T)
    for i in xrange(COUNTER):
        if i == 0:
            po().initialize_output(t[i])
    
        else: #i > 0:
            
            (Body1.x_neut,Body1.z_neut) = neutral_axis(Body1,Body1.xb_col[:Body1.N/2],DSTEP,TSTEP,t[i])
            panel_positions(Body1,P['S'],DSTEP,t[i])
            surface_kinematics(Body1,DSTEP,TSTEP,t[i],i,DEL_T)
            edge_shed(Body1,Edge1,i,DEL_T)
            wake_shed(Edge1,Wake1,i,DEL_T)
            
            influence_matrices(Body1,Edge1,Wake1,i)
            kutta(Body1,Edge1,Wake1,P['RHO'],i,DEL_T,P['SWITCH_KUTTA'])
            wake_rollup(Body1,Edge1,Wake1,P['DELTA_CORE'],i,DEL_T)
            
            #force(Body1,i)
            #po().solution_output(d_visc,cf,cl,ct,cpow,gamma)
             
            if np.fmod(i,10) == 0:
                po().timestep_header(i+1,t[i])
                po().solution_output(0,0,0,0,0,0)
                po().solution_complete_output(i/float(COUNTER-1)*100.)
    
    total_time = time.time()-start_time
    print "Simulation time:", np.round(total_time, 3), "seconds"

    graph.body_wake_plot(Body1,Edge1,Wake1)
    #graph.cp_plot(Body1)

if __name__ == '__main__':
    main()