from collections import namedtuple
import numpy as np

GeoVDVParameters = namedtuple('GeoVDVParameters', 'N C K EPSILON')
MotionParameters = namedtuple('MotionParameters', 'V0 THETA_MAX H_C F PHI')
SwimmerParameters = namedtuple('SimParameters', 'CE S SW_GEOMETRY SW_KUTTA')
BodyFrameCoordinates = namedtuple('BodyFrameCoordinates', 'x z x_col z_col')
AbsFrameCoordinates = namedtuple('AbsFrameCoordinates', 'x z x_col z_col x_mid z_mid x_neut z_neut')

class Swimmer(object):
    
    def __init__(self, SwimmerParameters, GeoParameters, MotionParameters):
        self.CE = SwimmerParameters.CE
        self.S = SwimmerParameters.S
        self.SW_GEOMETRY = SwimmerParameters.SW_GEOMETRY
        self.SW_KUTTA = SwimmerParameters.SW_KUTTA
        
        if self.SW_GEOMETRY == 'VDV':
            self.Body = Body.from_van_de_vooren(GeoParameters, MotionParameters)

class Body(object):
    
    def __init__(self, N, BodyFrameCoordinates, MotionParameters):
        
        self.N = N
        
        # Body-frame panel coordinates
        # x, z, x_col, z_col
        self.BF = BodyFrameCoordinates
        # Initialize absolute-frame panel coordinates
        # x, z, x_col, z_col, x_mid, z_mid, x_neut, z_neut
        self.AF = AbsFrameCoordinates(np.empty(N+1), np.empty(N+1), np.empty(N), np.empty(N),\
                                       np.zeros((3,N)), np.zeros((3,N)), np.empty(N/2), np.empty(N/2))
        
        # Prescribed motion
        self.MP = MotionParameters
        self.V0 = MotionParameters.V0 # gets used frequently enough
        
        self.vx = np.zeros(N)
        self.vz = np.zeros(N)
        
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
    def from_van_de_vooren(cls, GeoVDVParameters, MotionParameters):
        """
        Creates Body based on Van de Vooren airfoil geometry.
        Motion parameters are unused here, just getting passed through.
        """
    
        N = GeoVDVParameters.N
        C = GeoVDVParameters.C
        K = GeoVDVParameters.K
        EPSILON = GeoVDVParameters.EPSILON
        
        A = C*((1+EPSILON)**(K-1))*(2**(-K))
        
        phi = np.linspace(0,np.pi,N/2+1)
        
        r1 = np.sqrt((A*np.cos(phi)-A)**2+(A**2)*np.sin(phi)**2)
        r2 = np.sqrt((A*np.cos(phi)-EPSILON*A)**2+(A**2)*np.sin(phi)**2)
        
        phi1 = np.arctan2((A*np.sin(phi)) , (A*np.cos(phi)-A))
        phi2 = np.arctan2(A*np.sin(phi) ,(A*np.cos(phi)-EPSILON*A))
        
        x = ((r1**K)/(r2**(K-1)))*(np.cos(K*phi1)*np.cos((K-1)*phi2) + np.sin(K*phi1)*np.sin((K-1)*phi2))
        z_top = ((r1**K)/(r2**(K-1)))*(np.sin(K*phi1)*np.cos((K-1)*phi2) - np.cos(K*phi1)*np.sin((K-1)*phi2))
        z_bot = -((r1**K)/(r2**(K-1)))*(np.sin(K*phi1)*np.cos((K-1)*phi2) - np.cos(K*phi1)*np.sin((K-1)*phi2))
        
        x = x-x[-1] # Carrying the leading edge to the origin
        x[0] = C
        
        z_top[0] = 0
        z_bot[0] = 0
        z_bot[-1] = 0
        
        # Merge top and bottom surfaces together
        x = np.hstack((x , x[-2::-1]))
        z = np.hstack((z_bot , z_top[-2::-1]))
        
        x_col = ((x[1:] + x[:-1])/2)
        z_col = ((z[1:] + z[:-1])/2)
        
        BFC = BodyFrameCoordinates(x, z, x_col, z_col)
        
        return Body(N, BFC, MotionParameters)