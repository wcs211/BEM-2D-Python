import numpy as np
from general_functions import panel_vectors, transformation
    
def edge_shed(Body,Edge,i,DEL_T):
# Uses Body.(x, z, x_neut, z_neut, V0) and Edge.CE
# Gets Edge.(x, z)
# Others: none
    
    if i == 0:
        pass
    
    # The trailing edge panel whose length can be varied to optimize the solution
    else:
        Edge.x[0] = Body.AF.x[0]
        Edge.z[0] = Body.AF.z[0]
        Edge.x[1] = Body.AF.x[0] + Edge.CE*panel_vectors(Body.AF.x_neut,Body.AF.z_neut)[0][0]*Body.V0*DEL_T
        Edge.z[1] = Body.AF.z[0] + Edge.CE*panel_vectors(Body.AF.x_neut,Body.AF.z_neut)[1][0]*Body.V0*DEL_T
    
# Wake position calculations
def wake_shed(Edge,Wake,i,DEL_T):
# Uses Edge.(x, z, mu) and Wake.(x, z, V0, mu)
# Gets Wake.(x, z, mu, gamma)
# Others: none
    
    # First timestep is before t==0 (geometry initialization timestep)
    if i == 0:
        pass
    
    # Initialize wake coordinates when i==1
    elif i == 1:
        
        Wake.x[0] = Edge.x[-1]
        Wake.z[0] = Edge.z[-1]
        
        Wake.x[1:] = Wake.x[0] + np.arange(1,np.size(Wake.x))*(-Wake.V0)*DEL_T
        Wake.z[1:] = Wake.z[0]
    
    else:
        Wake.x[1:] = Wake.x[:-1]
        Wake.z[1:] = Wake.z[:-1]
        Wake.mu[1:] = Wake.mu[:-1]
        
        Wake.x[0] = Edge.x[-1]
        Wake.z[0] = Edge.z[-1]
        Wake.mu[0] = Edge.mu
        
        Wake.gamma[0] = -Wake.mu[0]
        Wake.gamma[1:-1] = Wake.mu[:-1]-Wake.mu[1:]
        Wake.gamma[-1] = Wake.mu[-1]

def wake_rollup(Body,Edge,Wake,DELTA_CORE,i,DEL_T):
# Uses Transformation
# Uses Body.(x, z, gamma, N, sigma)
# Uses Edge.(x, z, gamma, N)
# Uses Wake.(x, z, gamma)
# Gets Wake.(x, z)
# Others: NT, NI, xp1, xp2, zp, nx, nz, beta, dummy1, dummy2, dummy3, dummy4, xp, zp, r_b, r_e, r_w, vx, vz
    
    # Wake panels (not including the trailing edge panel) initialize when i==2
    if i < 2:
        pass
    
    else:
        
        NT = i-1 # Number of targets (wake panel points that are rolling up)
        
        vx = np.zeros(NT)
        vz = np.zeros(NT)
        
        # Coordinate transformation for body panels influencing wake
        (xp1,xp2,zp) = transformation(Wake.x[1:i],Wake.z[1:i],Body.AF.x,Body.AF.z)
        
        # Angle of normal vector with respect to global z-axis
        (nx,nz) = panel_vectors(Body.AF.x,Body.AF.z)[2:4]
        beta = np.arctan2(-nx,nz)
        
        # Katz-Plotkin eqns 10.20 and 10.21 for body source influence
        dummy1 = np.transpose(np.log((xp1**2+zp**2)/(xp2**2+zp**2))/(4*np.pi))
        dummy2 = np.transpose((np.arctan2(zp,xp2)-np.arctan2(zp,xp1))/(2*np.pi))
        
        # Rotate back to global coordinates
        dummy3 = dummy1*np.cos(beta) - dummy2*np.sin(beta)
        dummy4 = dummy1*np.sin(beta) + dummy2*np.cos(beta)
        
        # Finish eqns 10.20 and 10.21 for induced velocity by multiplying with Body.sigma
        vx = np.dot(dummy3,Body.sigma)
        vz = np.dot(dummy4,Body.sigma)
        
        # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
        NI = Body.N+1
        xp = np.repeat(Wake.x[1:i,np.newaxis].T,NI,0) - np.repeat(Body.AF.x[:,np.newaxis],NT,1)
        zp = np.repeat(Wake.z[1:i,np.newaxis].T,NI,0) - np.repeat(Body.AF.z[:,np.newaxis],NT,1)
        
        # Find distance r_b between each influence/target
        r_b = np.sqrt(xp**2+zp**2)
        
        # Katz-Plotkin eqns 10.9 and 10.10 for body doublet (represented as point vortices) influence
        dummy1 = np.transpose(zp/(2*np.pi*(r_b**2+DELTA_CORE**2)))
        dummy2 = np.transpose(-xp/(2*np.pi*(r_b**2+DELTA_CORE**2)))
        
        # Finish eqns 10.9 and 10.10 by multiplying with Body.gamma, add to induced velocity
        vx += np.dot(dummy1,Body.gamma)
        vz += np.dot(dummy2,Body.gamma)
        
        # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
        NI = Edge.N+1
        xp = np.repeat(Wake.x[1:i,np.newaxis].T,NI,0) - np.repeat(Edge.x[:,np.newaxis],NT,1)
        zp = np.repeat(Wake.z[1:i,np.newaxis].T,NI,0) - np.repeat(Edge.z[:,np.newaxis],NT,1)
        
        # Find distance r_e between each influence/target
        r_e = np.sqrt(xp**2+zp**2)
        
        # Katz-Plotkin eqns 10.9 and 10.10 for edge (as point vortices) influence
        dummy1 = np.transpose(zp/(2*np.pi*(r_e**2+DELTA_CORE**2)))
        dummy2 = np.transpose(-xp/(2*np.pi*(r_e**2+DELTA_CORE**2)))
        
        # Finish eqns 10.9 and 10.10 by multiplying with Edge.gamma, add to induced velocity
        vx += np.dot(dummy1,Edge.gamma)
        vz += np.dot(dummy2,Edge.gamma)
        
        # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
        NI = i
        xp = np.repeat(Wake.x[1:i,np.newaxis].T,NI,0) - np.repeat(Wake.x[:i,np.newaxis],NT,1)
        zp = np.repeat(Wake.z[1:i,np.newaxis].T,NI,0) - np.repeat(Wake.z[:i,np.newaxis],NT,1)
        
        # Find distance r_w between each influence/target
        r_w = np.sqrt(xp**2+zp**2)
        
        # Katz-Plotkin eqns 10.9 and 10.10 for wake (as point vortices) influence
        dummy1 = np.transpose(zp/(2*np.pi*(r_w**2+DELTA_CORE**2)))
        dummy2 = np.transpose(-xp/(2*np.pi*(r_w**2+DELTA_CORE**2)))
        
        # Finish eqns 10.9 and 10.10 by multiplying with Wake.gamma, add to induced velocity
        vx += np.dot(dummy1,Wake.gamma[:i])
        vz += np.dot(dummy2,Wake.gamma[:i])
        
        # Modify wake with the total induced velocity
        Wake.x[1:i] += vx*DEL_T
        Wake.z[1:i] += vz*DEL_T