import numpy as np

class Body(object):
    
    def __init__(self,N,V0,EPSILON,K,C,THETA_MAX,H_C,F,PHI,COUNTER):
        
        self.N = N
        self.V0 = V0
        self.EPSILON = EPSILON
        self.K = K
        self.C = C
        self.THETA_MAX = THETA_MAX
        self.H_C = H_C
        self.F = F
        self.PHI = PHI
        # positions in body frame of reference
        self.xb = np.zeros(N+1)
        self.zb = np.zeros(N+1)
        self.xb_col = np.zeros(N)
        self.zb_col = np.zeros(N)
        
        self.vx = np.zeros(N)
        self.vz = np.zeros(N)
        
        self.x = np.zeros(N+1)
        self.z = np.zeros(N+1)
        self.x_col = np.zeros(N)
        self.z_col = np.zeros(N)
        self.x_mid = np.zeros((3,N))
        self.z_mid = np.zeros((3,N))
        
        self.x_neut = np.zeros(N/2)
        self.z_neut = np.zeros(N/2)
        
        self.sigma = np.zeros(N)
        self.phi_s = np.zeros((N,N))
        self.phi_db = np.zeros((N,N))
        self.phi_dw = np.zeros(N)
        self.mu = np.zeros(N)
        self.gamma = np.zeros(N+1)
        
        self.p = np.zeros(N)
        self.cp = np.zeros(N)
        self.mu_past = np.zeros((2,N))
        
#        self.d = np.zeros(COUNTER-1)
#        self.l = np.zeros(COUNTER-1)