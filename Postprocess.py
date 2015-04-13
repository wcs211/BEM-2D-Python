import matplotlib.pyplot as plt
import numpy as np
import copy

# global figure variable
# this is to make sure each plot is drawn in a new window, no matter which plotting methods are used
N_fig=1

def BodyWakePlot(Body,Edge,Wake):
    
    global N_fig
    figure = plt.figure(N_fig)
    plt.clf()
    figure.add_subplot(1, 1, 1, axisbg='1') # change background color here
    plt.gca().set_aspect('equal')
    maxpercentile=95 # for truncating outliers
    
    # make color map based on vorticity
    color = copy.deepcopy(Wake.gamma[1:-1])
    # take a look at positive and negative circulations separately
    if np.min(color)<0: # check if negative circulations exist (in case of short simulations)
        # truncate any negative outliers
        color[color<np.percentile(color[color<0],100-maxpercentile)]=np.percentile(color[color<0],100-maxpercentile)
        # normalize negative circulations to [-1,0)
        color[color<0]=-color[color<0]/np.min(color)
    if np.max(color)>0: # check if positive circulations exist (in case of short simulations)
        # truncate any positive outliers
        color[color>np.percentile(color[color>0],maxpercentile)]=np.percentile(color[color>0],maxpercentile)
        # normalize positive circulations to (0,1]
        color[color>0]=color[color>0]/np.max(color)
    
    # scatter plot of wake points with red-white-blue colormap, as well as body outline and edge panel segment
    plt.scatter(Wake.X[1:-1],Wake.Z[1:-1], s=30, c=color, edgecolors='none', cmap=plt.get_cmap('bwr_r'))
    plt.plot(Body.X,Body.Z,'k')
    plt.plot(Edge.X,Edge.Z,'g')
    
    N_fig+=1
    
def CpPlot(Body):
    
    global N_fig
    figure = plt.figure(N_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # change background color here
    plt.gca().set_aspect('equal')
    plt.gca().invert_yaxis()
    
    plt.plot(Body.X_mid[:Body.N/2],Body.Cp[:Body.N/2],'g')
    plt.plot(Body.X_mid[Body.N/2:],Body.Cp[Body.N/2:],'b')
    plt.plot(Body.X,-Body.Z,'k')
    
    N_fig+=1
    
def DragVsPeriod(Body,rho,t):
    
    global N_fig
    figure = plt.figure(N_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # change background color here
    plt.xlabel('tau')
    plt.ylabel('Coefficent of drag')
    
    plt.plot(t[4:]*Body.f,-Body.D[3:]/(0.5*rho*Body.V0**2),'b')
    
    N_fig+=1
    
def LiftVsPeriod(Body,rho,t):
    
    global N_fig
    figure = plt.figure(N_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # change background color here
    plt.xlabel('tau')
    plt.ylabel('Coefficent of lift')
    
    plt.plot(t[4:]*Body.f,-Body.L[3:]/(0.5*rho*Body.V0**2),'g')
    
    N_fig+=1
    
def PressureDistrib(Body,rho):
    
    global N_fig
    figure = plt.figure(N_fig)
    figure.add_subplot(1, 1, 1, axisbg='1') # change background color here
    plt.xlabel('x/c')
    plt.ylabel('Pressure')
    plt.plot(np.hstack((Body.xcol,Body.xcol[::-1]))/Body.c,-Body.Pb[-1,:]/(0.5*rho*Body.V0**2),'r')
    
    N_fig+=1