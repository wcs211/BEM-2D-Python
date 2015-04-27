import numpy as np

class Body(object):
    
    def __init__(self,N,xb,zb,xb_col,zb_col,V0,THETA_MAX,H_C,F,PHI):
        
        self.N = N
        # positions in body frame of reference
        self.xb = xb
        self.zb = zb
        self.xb_col = xb_col
        self.zb_col = zb_col
        
        # prescribed motion
        self.V0 = V0
        self.THETA_MAX = THETA_MAX
        self.H_C = H_C
        self.F = F
        self.PHI = PHI
        
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
    
    @classmethod
    # VandeVooren airfoil geometry mapping
    def from_van_de_vooren(cls, N, C, K, EPSILON, V0,THETA_MAX,H_C,F,PHI):
    # Uses self.(C, EPSILON, K, N)
    # Gets self.(xb, zb, xb_col, zb_col)
    # Others: A, phi, r1, r2, phi1, phi2
        
        A = C*((1+EPSILON)**(K-1))*(2**(-K))
        
        phi = np.linspace(0,np.pi,N/2+1)
        
        r1 = np.sqrt((A*np.cos(phi)-A)**2+(A**2)*np.sin(phi)**2)
        r2 = np.sqrt((A*np.cos(phi)-EPSILON*A)**2+(A**2)*np.sin(phi)**2)
        
        phi1 = np.arctan2((A*np.sin(phi)) , (A*np.cos(phi)-A))
        phi2 = np.arctan2(A*np.sin(phi) ,(A*np.cos(phi)-EPSILON*A))
        
        xb = ((r1**K)/(r2**(K-1)))*(np.cos(K*phi1)*np.cos((K-1)*phi2) + np.sin(K*phi1)*np.sin((K-1)*phi2))
        zb_top = ((r1**K)/(r2**(K-1)))*(np.sin(K*phi1)*np.cos((K-1)*phi2) - np.cos(K*phi1)*np.sin((K-1)*phi2))
        zb_bot = -((r1**K)/(r2**(K-1)))*(np.sin(K*phi1)*np.cos((K-1)*phi2) - np.cos(K*phi1)*np.sin((K-1)*phi2))
        
        xb = xb-xb[-1] # Carrying the leading edge to the origin
        xb[0] = C
        
        zb_top[0] = 0
        zb_bot[0] = 0
        zb_bot[-1] = 0
        
        # Merge top and bottom surfaces together
        xb = np.hstack((xb , xb[-2::-1]))
        zb = np.hstack((zb_bot , zb_top[-2::-1]))
        
        xb_col = ((xb[1:] + xb[:-1])/2)
        zb_col = ((zb[1:] + zb[:-1])/2)
        
        return Body(N,xb,zb,xb_col,zb_col,V0,THETA_MAX,H_C,F,PHI)
    
#    #Flat plate geometry
#    def flat_plate(self):        
#        self.xf = np.linspace(1,0,self.N/2+1)
#        self.zf = np.zeros(self.N/2+1)