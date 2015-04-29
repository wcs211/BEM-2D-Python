import numpy as np

class Wake(object):
    
    def __init__(self,V0,N):
        
        self.N = N
        self.V0 = V0
        self.x = np.zeros(N+1)
        self.z = np.zeros(N+1)
        self.mu = np.zeros(N)
        self.gamma = np.zeros(N+1)