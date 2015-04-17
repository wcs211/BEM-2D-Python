import numpy as np
from Terminal_Output import printOutput as PO
from BodyClass import Body
from EdgeClass import Edge
from WakeClass import Wake
from Geometries import VanDeVooren
from Kinematics import NeutralAxis, SurfaceKinematics, EdgeShed, WakeShed, WakeRollup, CollocationPoints
from InfForce import InfluenceMatrices, Kutta
import Postprocess as post
import time # for timing simulation

PO().progTitle('1.0.0')
start_time=time.time()

Ce=0.4
counter=151

RF = 0.5 # reduced frequency
switch_Kutta = 1 # 0 for explicit, 1 for unsteady

N_Body=100
V0=-1.
epsilon=0.075
k=2.-(12.4/180)
c=1.
theta_max=5.0*np.pi/180
h_c=0
#f=0.5
f = RF/2/np.pi
phi=0
delta_core=(0.005*theta_max+0.09)*c
Body1=Body(N_Body,V0,epsilon,k,c,theta_max,h_c,f,phi,counter)
VanDeVooren(Body1)

Edge1=Edge(Body1.V0,Ce,counter)
Wake1=Wake(Body1,counter-2)

dstep=10**-5
tstep=10**-5
S=0.1
#delt=0.01
delt = (0.01/RF)*np.pi
rho=998.2
mu=0.001003
Re=rho*-V0*c/mu

t=delt*np.arange(0,counter)

PO().calcInput(theta_max/np.pi*180.,Re,theta_max/np.pi*180.,delt)

# data points per cycle == 1/(f*delt)
for i in np.arange(0,counter):
    if i==0:
        PO().initializeOutput(t[i])

    else: #i>0:
        
        (Body1.x_neut,Body1.z_neut)=NeutralAxis(Body1,Body1.x_col,dstep,tstep,t[i])
        SurfaceKinematics(Body1,dstep,tstep,t[i],i)
        CollocationPoints(Body1,S,i)
        EdgeShed(Body1,Edge1,i,delt)
        WakeShed(Edge1,Wake1,i,delt)
        
        InfluenceMatrices(Body1,Edge1,Wake1,dstep,tstep,S,i,counter)
        Kutta(Body1,Edge1,Wake1,rho,i,delt,switch_Kutta)
        WakeRollup(Body1,Edge1,Wake1,delta_core,i,delt)
        
        #Force(Body1,i)
        #PO().solutionOutput(D_visc,Cf,Cl,Ct,Cpow,Gamma)
         
#        if np.fmod(i,10)==0:
#            PO().timestepHeader(i+1,t[i])
#            PO().solutionOutput(0,0,0,0,0,0)
#            PO().solutionCompleteOutput(i/float(counter-1)*100.)

total_time=time.time()-start_time
print "Simulation time:", np.round(total_time, 3), "seconds"

#post.BodyWakePlot(Body1,Edge1,Wake1)
#post.CpPlot(Body1)
#plt.gca().invert_yaxis()
#post.DragVsPeriod(Body1,rho,t)