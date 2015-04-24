import time
import numpy as np
from terminal_output import print_output as po
from body_class import Body
from edge_class import Edge
from wake_class import Wake
from geometries import van_de_vooren
from kinematics import neutral_axis, panel_positions, surface_kinematics, edge_shed, wake_shed, wake_rollup
from inf_force import influence_matrices, kutta
import graphics as graph

def main():

    po().prog_title('1.0.0')
    start_time = time.time()
    
    CE = 0.4
    COUNTER = 401
    
    RF = np.pi # Reduced frequency
    SWITCH_KUTTA = 1 # 0 for explicit, 1 for unsteady
    
    N_BODY = 100
    V0 = -1.
    EPSILON = 0.075
    K = 2.-(12.4/180)
    C = 1.
    THETA_MAX = 5.73*np.pi/180
    H_C = 0
    #F = 0.5
    F = RF/(2*np.pi)
    PHI = 0
    DELTA_CORE = (0.005*THETA_MAX+0.09)*C
    Body1 = Body(N_BODY,V0,EPSILON,K,C,THETA_MAX,H_C,F,PHI,COUNTER)
    van_de_vooren(Body1)
    
    Edge1 = Edge(Body1.V0,CE,COUNTER)
    Wake1 = Wake(Body1,COUNTER-2)
    
    DSTEP = 10**-5
    TSTEP = 10**-5
    S = 0.1
    #DEL_T = 0.01
    DEL_T = (0.01/RF)*np.pi
    RHO = 998.2
    MU = 0.001003
    RE = RHO*-V0*C/MU

    t = DEL_T*np.arange(0,COUNTER)
    
    po().calc_input(THETA_MAX/np.pi*180.,RE,THETA_MAX/np.pi*180.,DEL_T)
    
    # Data points per cycle == 1/(F*DEL_T)
    for i in np.arange(0,COUNTER):
        if i == 0:
            po().initialize_output(t[i])
    
        else: #i > 0:
            
            (Body1.x_neut,Body1.z_neut) = neutral_axis(Body1,Body1.xb_col[:Body1.N/2],DSTEP,TSTEP,t[i])
            panel_positions(Body1,S,DSTEP,t[i])
            surface_kinematics(Body1,DSTEP,TSTEP,t[i],i,DEL_T)
            edge_shed(Body1,Edge1,i,DEL_T)
            wake_shed(Edge1,Wake1,i,DEL_T)
            
            influence_matrices(Body1,Edge1,Wake1,i)
            kutta(Body1,Edge1,Wake1,RHO,i,DEL_T,SWITCH_KUTTA)
            wake_rollup(Body1,Edge1,Wake1,DELTA_CORE,i,DEL_T)
            
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