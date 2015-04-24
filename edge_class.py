import numpy as np

class Edge(object):
    
    def __init__(self,V0,CE,COUNTER):
        
        self.N = 1
        self.CE = CE
        self.mu = np.zeros(self.N)
        self.x = np.zeros(self.N+1)
        self.z = np.zeros(self.N+1)
        self.gamma = np.zeros(self.N+1)