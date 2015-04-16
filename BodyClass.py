import numpy as np

class Body(object):
    
    def __init__(self,N,V0,epsilon,k,c,theta_max,h_c,f,phi,counter):
        
        self.N=N
        self.V0=V0
        self.epsilon=epsilon
        self.k=k
        self.c=c
        self.theta_max=theta_max
        self.h_c=h_c
        self.f=f
        self.phi=phi
        
        self.x=np.zeros(N/2)
        self.z_bot=np.zeros(N/2)
        self.z_top=np.zeros(N/2)
        self.x_col=np.zeros(N/2+1)
        self.z_col_top=np.zeros(N/2+1)
        self.z_col_bot=np.zeros(N/2+1)
        
        self.Vx=np.zeros(N)
        self.Vz=np.zeros(N)
        
        self.Xts=np.zeros(N/2+1)
        self.Xbs=np.zeros(N/2+1)
        self.Zts=np.zeros(N/2+1)
        self.Zbs=np.zeros(N/2+1)
        
        self.X=np.zeros(N+1)
        self.Z=np.zeros(N+1)
        self.X_mid=np.zeros(N)
        self.Z_mid=np.zeros(N)
        self.X_col=np.zeros(N)
        self.Z_col=np.zeros(N)
        
        self.x_neut=np.zeros(N/2)
        self.z_neut=np.zeros(N/2)
        
        self.sigma=np.zeros(N)
        self.phis=np.zeros((N,N))
        self.phidb=np.zeros((N,N))
        self.phidw=np.zeros((N,1)) # dynamic
        self.mu=np.zeros(N)
        self.gamma=np.zeros(N+1)
        
        self.P=np.zeros(N)
        self.Cp=np.zeros(N)
        self.mu_past=np.zeros((2,N))
        
#        self.D=np.zeros(counter-1)
#        self.L=np.zeros(counter-1)