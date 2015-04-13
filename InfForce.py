import numpy as np
from normalvectors import PanelVectors

# velocity and velocity potential equations are defined in panel coordinates so a transformation should be done
# each row of Xpleft/Xpright/Zp is an influence, and each column is a target
# Ni is N influences, Nt is N targets
# Xi/Zi is X/Z of influences, Xt/Zt is X/Z of target points
def Transformation(Xt,Zt,Xi,Zi):
# returns Xpleft, Xpright, Zp
# others: Nt, Ni, tx, tz, nx, nz, X2, x_tile, z_tile, tx_tile, tz_tile
    
    Nt=np.size(Xt)
    Ni=np.size(Xi)-1
    
    (tx,tz,nx,nz)=PanelVectors(Xi,Zi)[:-1]
    
    # intermediary variables to reduce number of tile/repeat operations
    # from normalvectors: tx==nz, tz==-nx
    x_tile = np.repeat(Xt[:,np.newaxis].T,Ni,0) - np.repeat(Xi[:-1,np.newaxis],Nt,1)
    z_tile = np.repeat(Zt[:,np.newaxis].T,Ni,0) - np.repeat(Zi[:-1,np.newaxis],Nt,1)
    tx_tile= np.repeat(tx[:,np.newaxis],Nt,1)
    tz_tile= np.repeat(tz[:,np.newaxis],Nt,1)
    
    # transforming left side collocation points from global to local coordinates
    Xpleft = x_tile*tx_tile + z_tile*tz_tile
    Zp = x_tile*(-tz_tile) + z_tile*tx_tile
    
    # transforming right side panel points into local coordinate system
    X2 = (Xi[1:]-Xi[:-1])*tx + (Zi[1:]-Zi[:-1])*tz
    Xpright = Xpleft - np.repeat(X2[:,np.newaxis],Nt,1)
    
    return (Xpleft,Xpright,Zp)

# this method constructs the influence coefficent matrices  
def InfluenceMatrices(Body,Edge,Wake,dstep,tstep,S,i,counter):
# uses Transformation
# uses Body.(N, X, Z, X_col, Z_col, V0, Vx, Vz)
# uses Edge.(X, Z)
# uses Wake.(X, Z)
# gets Body.(phis, sigma, phidb, phidw)
# others: Xpleft, Xpright, Zp, nx, nz, tx, tz, Xpleft_e, Xpright_e, Zp_e
    
    if i>0:
        # tangential and normal body panels and panel length calculations
        (Xpleft,Xpright,Zp)=Transformation(Body.X_col,Body.Z_col,Body.X,Body.Z)
    
        # body source singularities influencing the body
        # transpose so that row elements represent the effect on the (row number)th panel
        Body.phis=np.transpose( (Xpleft * np.log(Xpleft**2 + Zp**2) - Xpright * np.log(Xpright**2 + Zp**2) \
        + 2*Zp*(np.arctan2(Zp , Xpright) - np.arctan2(Zp , Xpleft)))/(4*np.pi) )
        
        # body source strength calculations
        (nx,nz)=PanelVectors(Body.X,Body.Z)[2:4]
        Body.sigma = nx*(Body.V0 + Body.Vx) + nz*Body.Vz   # normal vector pointing outward(overall sigma pointing outward)
        
        # body doublet singularities influencing body itself
        # and transpose similar to phis
        Body.phidb = np.transpose( -(np.arctan2(Zp , Xpright)\
        - np.arctan2(Zp , Xpleft))/(2*np.pi) )
        
        # edge doublet influencing the body
        (Xpleft,Xpright,Zp)=Transformation(Body.X_col,Body.Z_col,Edge.X,Edge.Z)
        
        if i>1: # no wake panels until i==2
            # wake doublets influencing the body
            (Xpleft_w,Xpright_w,Zp_w)=Transformation(Body.X_col,Body.Z_col,Wake.X[:i],Wake.Z[:i])
            # join edge and wake doublet influences into a single matrix
            Xpleft=np.insert(Xpleft, Edge.N, Xpleft_w, axis=0)
            Xpright=np.insert(Xpright, Edge.N, Xpright_w, axis=0)
            Zp=np.insert(Zp, Edge.N, Zp_w, axis=0)
    
        Body.phidw = np.transpose( -(np.arctan2(Zp , Xpright)\
        - np.arctan2(Zp , Xpleft))/(2*np.pi) )
        
import time
def Kutta(Body,Edge,Wake,i,delt,switch_Kutta):
# uses Body.(N, phidw, phis, sigma, phidb), Wake.mu, CpCalc()
# gets Edge.(mu, gamma), Body.(phidinv, mu, mu_past, gamma)
# others: n_iter, mu_guess, delta_Cp, Rhs, slope
    
    if i==0:
        pass
    else:
        n_iter=0
        mu_guess = np.empty(2) # [0] is current guess, [1] is previous
        delta_Cp = np.empty(2) # [0] is current delta_Cp, [1] is previous
        
        while True:
            n_iter+=1
            
            if n_iter==1:
                # begin with explicit Kutta condition as first guess
                # construct the augmented body matrix by combining body and trailing edge panel arrays
                C = np.zeros((Body.N,Body.N))
                C[:,0] = -Body.phidw[:,0] # C contains the explicit Kutta condition (first and last columns)
                C[:,-1] = Body.phidw[:,0]
                Body.phidinv=np.linalg.inv(Body.phidb + C)
                # get Rhs
                if i==1:
                    Rhs = -np.dot(Body.phis,Body.sigma)
                else:
                    Rhs = -np.dot(Body.phis,Body.sigma) - np.dot(Body.phidw[:,1:],Wake.mu[:i-1])
                # solve for body doublet strengths using explicit Kutta
                Body.mu=np.dot(Body.phidinv,Rhs)   
                # first mu_guess (from explicit Kutta)
                mu_guess[0] = Body.mu[-1]-Body.mu[0]
                
            else:
                if n_iter==2: # make a second initial guess
                    # update phidinv so it no longer includes explicit Kutta condition
                    Body.phidinv=np.linalg.inv(Body.phidb)
                    
                    mu_guess[1] = mu_guess[0]
                    delta_Cp[1] = delta_Cp[0]
                    mu_guess[0] = 0.8*mu_guess[1] # multiply first (explicit) guess by arbitrary constant to get second guess
                    
                else: # Newton method to get delta_Cp == 0
                    # get slope, which is change in delta_Cp divided by change in mu_guess
                    slope = (delta_Cp[0]-delta_Cp[1])/(mu_guess[0]-mu_guess[1])
                    mu_guess[1] = mu_guess[0]
                    delta_Cp[1] = delta_Cp[0]
                    mu_guess[0] = mu_guess[1] - delta_Cp[0]/slope
                
                # form right-hand side including mu_guess as an influence
                if i==1:
                    Rhs = -np.dot(Body.phis,Body.sigma) - np.squeeze(np.dot(Body.phidw,mu_guess[0]))
                else:              
                    Rhs = -np.dot(Body.phis,Body.sigma) - np.dot(Body.phidw,np.insert(Wake.mu[:i-1],0,mu_guess[0]))
                
                Body.mu=np.dot(Body.phidinv,Rhs)
                
            
            CpCalc(Body,i,delt)
            if switch_Kutta == 0:
                break
#            delta_Cp[0] = np.absolute(Body.Cp[-1]-Body.Cp[0])
#            print n_iter, delta_Cp[0]
            time.sleep(1)
            if delta_Cp[0]<0.0001:
                break
                    
        # mu_past used in differencing for CpCalc
        Body.mu_past[1,:] = Body.mu_past[0,:]
        Body.mu_past[0,:] = Body.mu
        
        Edge.mu = mu_guess[0]
        Edge.gamma[0]=-Edge.mu
        Edge.gamma[1]=Edge.mu
        
        # get gamma of body panels for use in wake rollup
        Body.gamma[0]=-Body.mu[0]
        Body.gamma[1:-1]=Body.mu[:-1]-Body.mu[1:]
        Body.gamma[-1]=Body.mu[-1]

#this method calculates the pressure coefficient for each body panel
#def CpCalc(Body,i,delt):
## uses Body.(N, X, Z)
## gets Body.Cp
## others: dmu_dl, PanelVectors, dmu_dt
#    
#    if i>0:
#        
#        (tx,tz,nx,nz,lpanel)=PanelVectors(Body.X,Body.Z)
#        
#        # tangential panel velocity dmu/dl, first-order differencing
#        dmu_dl=np.empty(Body.N)
#        dmu_dl[0]   = (Body.mu[1]-Body.mu[0]) / (lpanel[0]/2 + lpanel[1]/2)
#        dmu_dl[1:-1]= (Body.mu[2:]-Body.mu[:-2]) / (lpanel[:-2]/2 + lpanel[1:-1] + lpanel[2:]/2)
#        dmu_dl[-1]  = (Body.mu[-1]-Body.mu[-2]) / (lpanel[-2]/2 + lpanel[-1]/2)
#        
#        # potential change dmu/dt, second-order differencing after first time step
#        if i==1:
#            dmu_dt = (Body.mu - Body.mu_past[0,:])/delt
#        else:
#            dmu_dt = (3*Body.mu - 4*Body.mu_past[0,:] + Body.mu_past[1,:])/(2*delt)
#        
#        # Katz-Plotkin eqn 9.29 (steady Cp calculation, does not include dmu_dt)
##        Body.Cp = 1 - ( ((Body.V0+Body.Vx)*tx + Body.Vz*tz - dmu_dl)**2 ) / ((Body.V0 + Body.Vx)**2 + Body.Vz**2)
#        
#        # Katz-Plotkin eqn 13.28 (full unsteady Cp calculation)
#        Body.Cp = 1 - ( ((Body.V0+Body.Vx)*tx - dmu_dl)**2 + (Body.Vz*tz-Body.sigma)**2 + 2*dmu_dt ) / ((Body.V0 + Body.Vx)**2 + Body.Vz**2)

def CpCalc(Body,i,delt):
# uses Body.(N, X, Z)
# gets Body.Cp
# others: dmu_dl, PanelVectors, dmu_dt
    
    if i>0:
        
        (tx,tz,nx,nz,lpanel)=PanelVectors(Body.X,Body.Z)
        
        # tangential panel velocity dmu/dl, first-order differencing
        dmu_dl=np.empty(Body.N)
        dmu_dl[0]   = (Body.mu[0]-Body.mu[1]) / (lpanel[0]/2 + lpanel[1]/2)
        dmu_dl[1:-1]= (Body.mu[:-2]-Body.mu[2:]) / (lpanel[:-2]/2 + lpanel[1:-1] + lpanel[2:]/2)
        dmu_dl[-1]  = (Body.mu[-2]-Body.mu[-1]) / (lpanel[-2]/2 + lpanel[-1]/2)
        
        # potential change dmu/dt, second-order differencing after first time step
        if i==1:
            dmu_dt = (Body.mu - Body.mu_past[0,:])/delt
        else:
            dmu_dt = (3*Body.mu - 4*Body.mu_past[0,:] + Body.mu_past[1,:])/(2*delt)
        
        # unsteady Cp calculation from Matlab code
        Qinf = np.absolute(Body.V0)
        Qt = dmu_dl - (Body.V0 + Body.Vx)*tx - Body.Vz*tz
        
        Cp_s = 1. - Qt**2/Qinf**2
        Cp_us= 2.*dmu_dt/Qinf**2 + (Body.Vx**2 + Body.Vz**2)/Qinf**2 - 2.*Body.Vx/Qinf
        
        Body.Cp = Cp_s + Cp_us

def SingularityCalculations(Body,Edge,Wake,i):
# uses Body.(phis, sigma, phidb, phidw, phidinv), Wake.mu
# gets Body.mu
# others: Rhs
    
    # forming RHS and solving for doublet strengths
    if i==0:
        pass        
    elif i==1:
        Rhs = -np.dot(Body.phis,Body.sigma)
    else:
#        Rhs = -np.dot(Body.phis,Body.sigma) - np.dot(Body.phidw[:,1:],Wake.mu)
        Rhs = -np.dot(Body.phis,Body.sigma) - np.dot(Body.phidw,np.insert(Wake.mu[:i-1],0,Edge.mu))
        
    if i>0:
        Body.mu=np.dot(Body.phidinv,Rhs)
        
def KuttaExplicit(Body,Edge,i):
# uses Body.mu
# gets Edge.(mu, gamma), Body.gamma

    if i>0:
        #determining the trailing edge panel strength from the Kutta condition
        Edge.mu = Body.mu[-1] - Body.mu[0]
        Edge.gamma[0]=-Edge.mu
        Edge.gamma[1]=Edge.mu
    
        # get gamma of body panels for use in wake rollup
        Body.gamma[0]=-Body.mu[0]
        Body.gamma[1:-1]=Body.mu[:-1]-Body.mu[1:]
        Body.gamma[-1]=Body.mu[-1]

#this method calculates the pressure on top of each panel
def Pressure(Body,Edge,Wake,dstep,tstep,delt,rho,S,i,counter):
    
    if i>0:
        InfluenceMatrices(Body,Edge,Wake,dstep,tstep,S,i,counter)
        SingularityCalculations(Body,Edge,Wake,i)
        
        tx=PanelVectors(Body.X,Body.Z)[0]  ;  tz=PanelVectors(Body.X,Body.Z)[1]
        nx=PanelVectors(Body.X,Body.Z)[2]  ;  nz=PanelVectors(Body.X,Body.Z)[3]
        lpanel=PanelVectors(Body.X,Body.Z)[4]
    
        ##minus delmu over deldistance(calculating the tangential velocity on body panels)
        mu_v_one = (Body.mu[0] - Body.mu[1])/(lpanel[1]/2 + lpanel[0]/2)        
        mu_v = (Body.mu[0:-2] - Body.mu[2:])/(lpanel[0:-2]/2 +lpanel[1:-1] + lpanel[2:]/2)
        mu_v_N = (Body.mu[-2] - Body.mu[-1])/(lpanel[-2]/2 + lpanel[-1]/2)
        mu_v = np.vstack((np.reshape(mu_v,(Body.N-2,1)) , mu_v_N))
        mu_v = np.vstack((mu_v_one , mu_v))
        
##first order differencing (differencing wrt time)
#            dummy[1,:] = mub
#
#            mu_del_t = (dummy[1,:] - dummy[0,:])/delt
#            
#            Pb[i-1,:] = -rho*mu_del_t + rho*(((Vx+V0)*nx + Vz*nz)*sigma \
#            + ((Vx+V0)*tx+Vz*tz)*np.reshape(mu_v,(1,N)))\
#            - rho*(sigma**2 + np.reshape(mu_v,(1,N))**2)/2
#         
#            dummy[0,:]=dummy[1,:]
        
#second order differencing (differencing wrt time)
        if i==1:
            Body.dummy[1,:]=Body.mu
            mu_del_t = (Body.dummy[1,:] - Body.dummy[0,:])/delt
            
            Body.Pb[i-1,:] = -rho*mu_del_t + rho*(((Body.Vx+Body.V0)*nx + Body.Vz*nz)*Body.sigma \
            + ((Body.Vx+Body.V0)*tx+Body.Vz*tz)*np.reshape(mu_v,(1,Body.N)))\
            - rho*(Body.sigma**2 + np.reshape(mu_v,(1,Body.N))**2)/2
        
        else:
            Body.dummy[2,:]=Body.mu
            mu_del_t = (3*Body.dummy[2,:] - 4*Body.dummy[1,:] + Body.dummy[0,:])/(2*delt)
            
            Body.Pb[i-1,:] = -rho*mu_del_t + rho*(((Body.Vx+Body.V0)*nx + Body.Vz*nz)*Body.sigma \
            + ((Body.Vx+Body.V0)*tx+Body.Vz*tz)*np.reshape(mu_v,(1,Body.N)))\
            - rho*(Body.sigma**2 + np.reshape(mu_v,(1,Body.N))**2)/2
            
            Body.dummy[0,:]=Body.dummy[1,:]
            Body.dummy[1,:]=Body.dummy[2,:]
       
#calculation of drag and lift forces affecting overall airfoil      
def Force(Body,i):
    
    if i==0:
        pass

    else:
        tx=PanelVectors(Body.X,Body.Z)[0]  ;  tz=PanelVectors(Body.X,Body.Z)[1]
        nx=PanelVectors(Body.X,Body.Z)[2]  ;  nz=PanelVectors(Body.X,Body.Z)[3]
        lpanel=PanelVectors(Body.X,Body.Z)[4]
                              
        Body.D[i-1] = np.dot(Body.Pb[i-1,:]*lpanel , np.reshape(tx,(Body.N,1)))\
                    + np.dot(Body.Pb[i-1,:]*lpanel , np.reshape(-tz,(Body.N,1)))

        Body.L[i-1] = np.dot(Body.Pb[i-1,:]*lpanel , np.reshape(-nz,(Body.N,1)))\
                    + np.dot(Body.Pb[i-1,:]*lpanel , np.reshape(nx,(Body.N,1)))