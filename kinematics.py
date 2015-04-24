import numpy as np
from normal_vectors import panel_vectors, point_vectors
from inf_force import transformation

# Neutral axis is the axis which coincides with the chord line and divides the symmetric airfoil into two
def neutral_axis(Body,x,DSTEP,TSTEP,t):
# Uses Body.(THETA_MAX, F, PHI)
# Returns x_neut, z_neut (just pitching for now)
    
    # Pitching motion of the airfoil(x and z points of the neutral axis due to pitching motion)
    x_neut_p = (x+DSTEP)*np.cos(Body.THETA_MAX*np.sin(2*np.pi*Body.F*(t+TSTEP) + Body.PHI))
    z_neut_p = (x+DSTEP)*np.sin(Body.THETA_MAX*np.sin(2*np.pi*Body.F*(t+TSTEP) + Body.PHI))

    # Heaving motion of the airfoil(x and z points of the neutral axis due to the heaving motion)
    # xneuth = x+DSTEP
    # zneuth = Body.H_C*Body.C*np.sin(2*np.pi*Body.F*(t+TSTEP))*np.ones(len(x))
    
    # Burst and coast behaviour of the airfoil
    ####################

    # Overall neutral axis position
    return(x_neut_p,z_neut_p)
    
# Gets the panel endpoint positions and collocation point positions of the body
def panel_positions(Body,S,DSTEP,t):
# Uses neutral_axis()
# Uses Body.(xb, zb, zb_col, V0)
# Gets Body.(x, z, x_col, z_col, x_mid, z_mid)
# Others: x_neut, z_neut, xdp_s, zdp_s, xdm_s, zdm_s
    
    # Body surface panel endpoint calculations
    (x_neut,z_neut) = neutral_axis(Body,Body.xb,0,0,t)
    
    # Infinitesimal differences on the neutral axis to calculate the tangential and normal vectors
    (xdp_s,zdp_s) = neutral_axis(Body,Body.xb,DSTEP,0,t)
    (xdm_s,zdm_s) = neutral_axis(Body,Body.xb,-DSTEP,0,t)
    
    # Displaced airfoil's surface points for time t0
    Body.x = x_neut + point_vectors(xdp_s,xdm_s,zdp_s,zdm_s)[2]*Body.zb + Body.V0*t
    Body.z = z_neut + point_vectors(xdp_s,xdm_s,zdp_s,zdm_s)[3]*Body.zb
    
    # Body panel midpoint positions
    # Store past values for backwards differencing to get surface velocities
    Body.x_mid[1:,:] = Body.x_mid[:-1,:]
    Body.z_mid[1:,:] = Body.z_mid[:-1,:]
    # Get current midpoint positions
    Body.x_mid[0,:] = (Body.x[:-1]+Body.x[1:])/2
    Body.z_mid[0,:] = (Body.z[:-1]+Body.z[1:])/2
    
    # Collocation points are the points where impermeable boundary condition is forced
    # They should be shifted inside or outside of the boundary depending on the dirichlet or neumann condition
    # Shifting surface collocation points some percent of the height from the neutral axis
    # Normal vectors point outward but positive S is inward, so the shift must be subtracted from the panel midpoints
    Body.x_col = Body.x_mid[0,:] - S*panel_vectors(Body.x,Body.z)[2]*np.absolute(Body.zb_col)
    Body.z_col = Body.z_mid[0,:] - S*panel_vectors(Body.x,Body.z)[3]*np.absolute(Body.zb_col)
    
# This method calculates the actual surface positions of the airfoil for each time step
# Using the neutral axis and appropriate normal vectors of each point on the neutral axis
# This class also calculates the velocity of the panel midpoints for each time step
def surface_kinematics(Body,DSTEP,TSTEP,t,i,DEL_T):
# Uses neutral_axis(), point_vectors()
# Uses Body.(xb_col, zb_col, x_mid, z_mid, V0)
# Gets Body.(vx, vz)
# Others: xtpneut, ztpneut, xtpdp, ztpdp, xtpdm, ztpdm, xtmneut, ztmneut, xtmdp, ztmdp, xtmdm, ztmdm, xctp, xctm, zctp, zctm
    
    if i == 1:
        # Panel midpoint velocity calculations
        # Calculating the surface positions at tplus(tp) and tminus(tm) for every timestep
        (xtpneut,ztpneut) = neutral_axis(Body,Body.xb_col,0,TSTEP,t)
        (xtpdp,ztpdp) = neutral_axis(Body,Body.xb_col,DSTEP,TSTEP,t)
        (xtpdm,ztpdm) = neutral_axis(Body,Body.xb_col,-DSTEP,TSTEP,t)
        (xtmneut,ztmneut) = neutral_axis(Body,Body.xb_col,0,-TSTEP,t)
        (xtmdp,ztmdp) = neutral_axis(Body,Body.xb_col,DSTEP,-TSTEP,t)
        (xtmdm,ztmdm) = neutral_axis(Body,Body.xb_col,-DSTEP,-TSTEP,t)
        
        # Displaced airfoil's panel midpoints for times tplus(tp) and tminus(tm)      
        xctp = xtpneut + point_vectors(xtpdp,xtpdm,ztpdp,ztpdm)[2]*Body.zb_col + Body.V0*t
        xctm = xtmneut + point_vectors(xtmdp,xtmdm,ztmdp,ztmdm)[2]*Body.zb_col + Body.V0*t
            
        zctp = ztpneut + point_vectors(xtpdp,xtpdm,ztpdp,ztpdm)[3]*Body.zb_col
        zctm = ztmneut + point_vectors(xtmdp,xtmdm,ztmdp,ztmdm)[3]*Body.zb_col
        
        # Velocity calculations on the surface panel midpoints
        Body.vx = (xctp - xctm)/(2*TSTEP)
        Body.vz = (zctp - zctm)/(2*TSTEP)
        
    elif i == 2:
        # First-order backwards differencing of body collocation point positions
        Body.vx = (Body.x_mid[0,:]-Body.x_mid[1,:])/DEL_T - Body.V0
        Body.vz = (Body.z_mid[0,:]-Body.z_mid[1,:])/DEL_T
    
    else:
        # Second-order backwards differencing of body collocation point positions
        Body.vx = (3*Body.x_mid[0,:]-4*Body.x_mid[1,:]+Body.x_mid[2,:])/(2*DEL_T) - Body.V0
        Body.vz = (3*Body.z_mid[0,:]-4*Body.z_mid[1,:]+Body.z_mid[2,:])/(2*DEL_T)
    
def edge_shed(Body,Edge,i,DEL_T):
# Uses Body.(x, z, x_neut, z_neut, V0) and Edge.CE
# Gets Edge.(x, z)
# Others: none
    
    if i == 0:
        pass
    
    # The trailing edge panel whose length can be varied to optimize the solution
    else:
        Edge.x[0] = Body.x[0]
        Edge.z[0] = Body.z[0]
        Edge.x[1] = Body.x[0] + Edge.CE*panel_vectors(Body.x_neut,Body.z_neut)[0][0]*Body.V0*DEL_T
        Edge.z[1] = Body.z[0] + Edge.CE*panel_vectors(Body.x_neut,Body.z_neut)[1][0]*Body.V0*DEL_T
    
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
        (xp1,xp2,zp) = transformation(Wake.x[1:i],Wake.z[1:i],Body.x,Body.z)
        
        # Angle of normal vector with respect to global z-axis
        (nx,nz) = panel_vectors(Body.x,Body.z)[2:4]
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
        xp = np.repeat(Wake.x[1:i,np.newaxis].T,NI,0) - np.repeat(Body.x[:,np.newaxis],NT,1)
        zp = np.repeat(Wake.z[1:i,np.newaxis].T,NI,0) - np.repeat(Body.z[:,np.newaxis],NT,1)
        
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