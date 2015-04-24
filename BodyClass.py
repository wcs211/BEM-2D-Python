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
        
        self.x=np.zeros(N+1)
        self.z=np.zeros(N+1)
        self.x_col=np.zeros(N)
        self.z_col=np.zeros(N)
        
        self.Vx=np.zeros(N)
        self.Vz=np.zeros(N)
        
        self.X=np.zeros(N+1)
        self.Z=np.zeros(N+1)
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