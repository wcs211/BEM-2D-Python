import time
import numpy as np
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from functions_influence import quilt, wake_rollup, u_psi
from terminal_output import print_output as po
import functions_graphics as graph
from SolidClass import solid
from PyFEA import PyFEA
from FSIClass import FSI
from functions_general import panel_vectors, archive, absoluteToBody

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

    SPLIT = 0.4
    SwiP = PC.SwimmerParameters(P['CE'], P['DELTA_CORE'], P['SW_KUTTA'], P['SW_WAKE'])
    GeoP = PC.GeoVDVParameters(P['N_BODY'], P['S'], P['C'], P['K'], P['EPSILON'])
#    GeoP = PC.GeoFPParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
#    GeoP = PC.GeoTDParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
    MotP1 = PC.MotionParameters(0., 0*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
    MotP2 = PC.MotionParameters(0., 1*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
    MotP3 = PC.MotionParameters(0., 2*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
    MotP4 = PC.MotionParameters(0., 3*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])
    MotP5 = PC.MotionParameters(0., 4*SPLIT, P['V0'], P['THETA_MAX'], P['F'], P['PHI'])

    S1 = Swimmer(SwiP, GeoP, MotP1, COUNTER)
    S2 = Swimmer(SwiP, GeoP, MotP2, COUNTER)
    S3 = Swimmer(SwiP, GeoP, MotP3, COUNTER)
    S4 = Swimmer(SwiP, GeoP, MotP4, COUNTER)
    S5 = Swimmer(SwiP, GeoP, MotP5, COUNTER)
    Swimmers = [S1]

#    Solid1 = solid(S1.Body, P['N_ELEMENTS_S'], 0.001)
#    FSI1 = FSI(S1.Body, Solid1)
#    PyFEA1 = PyFEA(Solid1, P['FRAC_DELT'], P['DEL_T'], P['E'], P['RHO_S'])

    po().calc_input(MotP1.THETA_MAX/np.pi*180.,RE,MotP1.THETA_MAX/np.pi*180.,DEL_T)

#    Solid1.initMesh()
#    Solid1.initThinPlate(P['T_MAX'],P['C'],P['SWITCH_CNST_THK_BM'],P['T_CONST'],P['FLEX_RATIO'])
#    Solid1.initTearDrop(P['T_MAX'],P['C'],P['SWITCH_CNST_THK_BM'],P['T_CONST'],P['FLEX_RATIO'])

    # Data points per cycle == 1/(F*DEL_T)
    for i in xrange(COUNTER):
        if P['SW_FSI'] == 1:
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
                theta = S1.Body.MP.THETA_MAX * np.sin(2 * np.pi * S1.Body.MP.F * (T[i] + TSTEP) + S1.Body.MP.PHI)
                Solid1.nodes[:,0] = (Solid1.nodesNew[:,0] - Solid1.nodesNew[0,0])*np.cos(theta)
                Solid1.nodes[:,1] = (Solid1.nodesNew[:,0] - Solid1.nodesNew[0,0])*np.sin(theta)
                graph.plot_n_go(S1.Body)

            else:
                if np.fmod(i,P['VERBOSITY']) == 0:
                    po().timestep_header(i,T[i])
                    po().fsi_header()

                FSI1.readFsiControls(P['FIXED_PT_RELAX'], P['N_OUTERCORR_MAX'])
                FSI1.__init__(S1.Body, Solid1)
                outerCorr = 0
                while True:
                    outerCorr += 1
                    FSI1.setInterfaceDisplacemet(outerCorr, P['COUPLING_SCHEME'])
                    for Swim in Swimmers:
                        if (outerCorr == 1):
                            Swim.Body.panel_positions(DSTEP, T[i])
                        else:
                            Swim.Body.AF.x += (FSI1.fluidNodeDispl[:,0] - FSI1.fluidNodeDisplOld[:,0])
                            Swim.Body.AF.z += (FSI1.fluidNodeDispl[:,1] - FSI1.fluidNodeDisplOld[:,1])
                            Swim.Body.AF.x_mid[0,:] = (Swim.Body.AF.x[:-1]+Swim.Body.AF.x[1:])/2
                            Swim.Body.AF.z_mid[0,:] = (Swim.Body.AF.z[:-1]+Swim.Body.AF.z[1:])/2
                            theta = Swim.Body.MP.THETA_MAX * np.sin(2 * np.pi * Swim.Body.MP.F * (T[i] + TSTEP) + Swim.Body.MP.PHI)
                            Swim.Body.BF.x = (Swim.Body.AF.x - Swim.Body.AF.x_le) * np.cos(-1*theta) - (Swim.Body.AF.z - Swim.Body.AF.z_le) * np.sin(-1*theta)
                            Swim.Body.BF.z = (Swim.Body.AF.z - Swim.Body.AF.z_le) * np.cos(-1*theta) + (Swim.Body.AF.x - Swim.Body.AF.x_le) * np.sin(-1*theta)
                            Swim.Body.BF.x_col = ((Swim.Body.BF.x[1:] + Swim.Body.BF.x[:-1])/2)
                            Swim.Body.BF.z_col = ((Swim.Body.BF.z[1:] + Swim.Body.BF.z[:-1])/2)
                            Swim.Body.AF.x_col = Swim.Body.AF.x_mid[0,:] - Swim.Body.S*panel_vectors(Swim.Body.AF.x, Swim.Body.AF.z)[2]*np.absolute(Swim.Body.BF.z_col)
                            Swim.Body.AF.z_col = Swim.Body.AF.z_mid[0,:] - Swim.Body.S*panel_vectors(Swim.Body.AF.x, Swim.Body.AF.z)[3]*np.absolute(Swim.Body.BF.z_col)

                        Swim.Body.surface_kinematics(DSTEP, TSTEP, DEL_T, T[i], i)
                        if (outerCorr == 1):
                            Swim.edge_shed(DEL_T, i)
                            Swim.wake_shed(DEL_T, i)

                    quilt(Swimmers, RHO, DEL_T, i)

                    #TODO: Replace '0' with viscous drag component when available
                    FSI1.setInterfaceForce(Solid1, S1.Body, PyFEA1, T[i], TSTEP, outerCorr,
                                  P['SWITCH_VISC_DRAG'], 0, P['SWITCH_INTERP_MTD'], P['C'], i)
                    PyFEA1.solve(S1.Body, Solid1, outerCorr, T[i], TSTEP, P['M_TYPE'], P['INT_METHOD'], P['ALPHA'], P['BETA'], P['GAMMA'])
                    FSI1.getDisplacements(Solid1, S1.Body, PyFEA1, T[i], TSTEP, P['SWITCH_INTERP_MTD'], P['FLEX_RATIO'])
                    FSI1.calcFSIResidual(Solid1, outerCorr)

                    if np.fmod(i,P['VERBOSITY']) == 0:
                        po().fsi_iter_out(outerCorr,FSI1.fsiRelaxationFactor,FSI1.maxDU,FSI1.maxMagFsiResidual,FSI1.fsiResidualNorm,FSI1.maxFsiResidualNorm)

                    if (FSI1.fsiResidualNorm <= P['OUTER_CORR_TOL'] or outerCorr >= P['N_OUTERCORR_MAX']):
                        if (FSI1.fsiResidualNorm <= P['OUTER_CORR_TOL']):
                            po().fsi_converged()
                        else:
                            po().fsi_not_converged()
                        if np.fmod(i,P['VERBOSITY']) == 0:
                            po().solution_output(0,0,0,0,0,0)
                            po().solution_complete_output(i/float(COUNTER-1)*100.)
                        wake_rollup(Swimmers, DEL_T, i)
                        absoluteToBody(S1.Body, Solid1, T[i], TSTEP)
                        archive(S1.Body.AF.x_mid)
                        archive(S1.Body.AF.z_mid)
                        graph.plot_n_go(S1.Body)
                        break
        else:
            if i == 0:
                po().initialize_output(T[i])

                for Swim in Swimmers:
                        Swim.Body.panel_positions(DSTEP, T[i])
                        Swim.Body.surface_kinematics(DSTEP, TSTEP, DEL_T, T[i], i)
                        Swim.Body.calc_sigma(Swim.SW_WAKE, Swim.Body.u_psi)
                        Swim.edge_shed(Swim.SW_WAKE, DEL_T, i)
                        Swim.wake_shed(Swim.SW_WAKE, DEL_T, i)
                quilt(Swimmers, RHO, DEL_T, i)
                archive(S1.Body.AF.x_mid)
                archive(S1.Body.AF.z_mid)
#                graph.plot_n_go(S1.Body)
            else:
                if np.fmod(i,P['VERBOSITY']) == 0:
                    po().timestep_header(i,T[i])

                for Swim in Swimmers:
                    Swim.Body.panel_positions(DSTEP, T[i])
                    Swim.Body.surface_kinematics(DSTEP, TSTEP, DEL_T, T[i], i)
                    Swim.Body.calc_sigma(Swim.SW_WAKE, Swim.Body.u_psi)
                    Swim.edge_shed(Swim.SW_WAKE, DEL_T, i)
                    Swim.wake_shed(Swim.SW_WAKE, DEL_T, i)
                quilt(Swimmers, RHO, DEL_T, i)

                if np.fmod(i,P['VERBOSITY']) == 0:
                    po().solution_output(0,0,0,0,0,0)
                    po().solution_complete_output(i/float(COUNTER-1)*100.)
                wake_rollup(Swimmers, DEL_T, i)
                u_psi(Swimmers, i)
                archive(S1.Body.AF.x_mid)
                archive(S1.Body.AF.z_mid)
#                graph.plot_n_go(S1.Body)


    total_time = time.time()-start_time
    print "Simulation time:", np.round(total_time, 3), "seconds"

    graph.body_wake_plot(Swimmers)
#    graph.cp_plot(S1.Body)

if __name__ == '__main__':
    main()