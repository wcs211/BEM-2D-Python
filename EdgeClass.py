import numpy as np

class Edge(object):
    
    def __init__(self,V0,Ce,counter):
        
        self.N=1
        self.Ce=Ce
        self.mu=np.zeros(self.N)
        self.X=np.zeros(self.N+1)
        self.Z=np.zeros(self.N+1)
        self.gamma=np.zeros(self.N+1)