import time
import numpy as np
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from functions_influence import quilt, wake_rollup
from terminal_output import print_output as po
import functions_graphics as graph
from functions_general import archive

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
    MotP = PC.MotionParameters(0., 0., P['V0'], P['THETA_MAX'], P['F'], P['PHI'])

    S1 = Swimmer(SwiP, GeoP, MotP, COUNTER-1)
    Swimmers = [S1]

    po().calc_input(MotP.THETA_MAX/np.pi*180.,RE,MotP.THETA_MAX/np.pi*180.,DEL_T)

    # Data points per cycle == 1/(F*DEL_T)
    for i in xrange(COUNTER):
        if i == 0:
            po().initialize_output(T[i])

        for Swim in Swimmers:
            Swim.Body.panel_positions(DSTEP, T[i])
            Swim.Body.surface_kinematics(DSTEP, TSTEP, DEL_T, T[i], i)
            Swim.edge_shed(DEL_T, i)
            Swim.wake_shed(DEL_T, i)
        quilt(Swimmers, RHO, DEL_T, i)
        wake_rollup(Swimmers, DEL_T, i)
        archive(S1.Body.AF.x_mid)
        archive(S1.Body.AF.z_mid)

        if np.fmod(i,10) == 0:
            po().timestep_header(i+1, T[i])
            po().solution_output(0,0,0,0,0,0)
            po().solution_complete_output(i/float(COUNTER-1)*100.)

    total_time = time.time()-start_time
    print "Simulation time:", np.round(total_time, 3), "seconds"

    graph.body_wake_plot(Swimmers)
    graph.cp_plot(S1.Body)

if __name__ == '__main__':
    main()