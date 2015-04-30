from collections import namedtuple
import numpy as np
from general_functions import point_vectors, panel_vectors, archive

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
        
        # Body-frame panel coordinates:
        # x, z, x_col, z_col
        self.BF = BodyFrameCoordinates
        # Initialize absolute-frame panel coordinates:
        # x, z, x_col, z_col, x_mid, z_mid, x_neut, z_neut
        self.AF = AbsFrameCoordinates(np.empty(N+1), np.empty(N+1), np.empty(N), np.empty(N),\
                                       np.zeros((3,N)), np.zeros((3,N)), np.empty(N/2), np.empty(N/2))
        # Prescribed motion
        self.MP = MotionParameters
        self.V0 = MotionParameters.V0 # gets used frequently
        
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
        MotionParameters are unused here, just getting passed through for Body creation.
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
        
    # Updates the absolute-frame coordinates of the body
    def panel_positions(self, S, DSTEP, t):
    # Uses neutral_axis()
    # Uses Body.(xb, zb, zb_col, V0)
    # Gets Body.(x, z, x_col, z_col, x_mid, z_mid)
    # Others: x_neut, z_neut, xdp_s, zdp_s, xdm_s, zdm_s
    
        (bfx, bfz, bfx_col, bfz_col) = self.BF
        V0 = self.V0
        
        # Body surface panel endpoint calculations
        (x_neut, z_neut) = self.neutral_axis(bfx, 0, 0, t)
        
        # Infinitesimal differences on the neutral axis to calculate the tangential and normal vectors
        (xdp_s, zdp_s) = self.neutral_axis(bfx, DSTEP, 0, t)
        (xdm_s, zdm_s) = self.neutral_axis(bfx, -DSTEP, 0, t)
        
        # Absolute-frame panel endpoint positions for time t
        afx = x_neut + point_vectors(xdp_s, xdm_s, zdp_s, zdm_s)[2]*bfz + V0*t
        afz = z_neut + point_vectors(xdp_s, xdm_s, zdp_s, zdm_s)[3]*bfz
        
        # Absolute-frame panel midpoint positions
        x_mid = (afx[:-1]+afx[1:])/2
        z_mid = (afz[:-1]+afz[1:])/2
        
        # Collocation points are the points where impermeable boundary condition is forced
        # They should be shifted inside or outside of the boundary depending on the dirichlet or neumann condition
        # Shifting surface collocation points some percent of the height from the neutral axis
        # Normal vectors point outward but positive S is inward, so the shift must be subtracted from the panel midpoints
        afx_col = x_mid - S*panel_vectors(afx, afz)[2]*np.absolute(bfz_col)
        afz_col = z_mid - S*panel_vectors(afx, afz)[3]*np.absolute(bfz_col)
        
        # Archive past values of x_mid and z_mid for differencing wrt time later
        archive(self.AF.x_mid)
        archive(self.AF.z_mid)
        
        self.AF.x[:] = afx
        self.AF.z[:] = afz
        self.AF.x_col[:] = afx_col
        self.AF.z_col[:] = afz_col
        self.AF.x_mid[0,:] = x_mid
        self.AF.z_mid[0,:] = z_mid
        
    # Neutral axis is the axis which coincides with the chord line and divides the symmetric airfoil into two
    def neutral_axis(self, x, DSTEP, TSTEP, t):
    # Uses Body.(THETA_MAX, F, PHI)
    # Returns x_neut, z_neut (just pitching for now)
        
        THETA_MAX = self.MP.THETA_MAX
        F = self.MP.F
        PHI = self.MP.PHI
        
        # Pitching motion of the airfoil(x and z points of the neutral axis due to pitching motion)
        x_neut = (x+DSTEP)*np.cos(THETA_MAX*np.sin(2*np.pi*F*(t+TSTEP) + PHI))
        z_neut = (x+DSTEP)*np.sin(THETA_MAX*np.sin(2*np.pi*F*(t+TSTEP) + PHI))
    
        # Overall neutral axis position
        return(x_neut, z_neut)
        
    # This method calculates the actual surface positions of the airfoil for each time step
    # Using the neutral axis and appropriate normal vectors of each point on the neutral axis
    # This class also calculates the velocity of the panel midpoints for each time step
    def surface_kinematics(self, DSTEP, TSTEP, t, i, DEL_T):
    # Uses neutral_axis(), point_vectors()
    # Uses Body.(xb_col, zb_col, x_mid, z_mid, V0)
    # Gets Body.(vx, vz)
    # Others: xtpneut, ztpneut, xtpdp, ztpdp, xtpdm, ztpdm, xtmneut, ztmneut, xtmdp, ztmdp, xtmdm, ztmdm, xctp, xctm, zctp, zctm
            
        if i == 1:
            
            x_col = self.BF.x_col
            z_col = self.BF.z_col
            
            # Panel midpoint velocity calculations
            # Calculating the surface positions at tplus(tp) and tminus(tm) for every timestep
            (xtpneut, ztpneut) = self.neutral_axis(x_col, 0, TSTEP, t)
            (xtpdp, ztpdp) = self.neutral_axis(x_col, DSTEP, TSTEP, t)
            (xtpdm, ztpdm) = self.neutral_axis(x_col, -DSTEP, TSTEP, t)
            (xtmneut, ztmneut) = self.neutral_axis(x_col, 0, -TSTEP, t)
            (xtmdp, ztmdp) = self.neutral_axis(x_col, DSTEP, -TSTEP, t)
            (xtmdm, ztmdm) = self.neutral_axis(x_col, -DSTEP, -TSTEP, t)
            
            # Displaced airfoil's panel midpoints for times tplus(tp) and tminus(tm)      
            xctp = xtpneut + point_vectors(xtpdp, xtpdm, ztpdp, ztpdm)[2]*z_col + self.V0*t
            xctm = xtmneut + point_vectors(xtmdp, xtmdm, ztmdp, ztmdm)[2]*z_col + self.V0*t
                
            zctp = ztpneut + point_vectors(xtpdp, xtpdm, ztpdp, ztpdm)[3]*z_col
            zctm = ztmneut + point_vectors(xtmdp, xtmdm, ztmdp, ztmdm)[3]*z_col
            
            # Velocity calculations on the surface panel midpoints
            self.vx = (xctp - xctm)/(2*TSTEP)
            self.vz = (zctp - zctm)/(2*TSTEP)
            
        elif i == 2:
            # First-order backwards differencing of body collocation point positions
            self.vx = (self.AF.x_mid[0,:]-self.AF.x_mid[1,:])/DEL_T - self.V0
            self.vz = (self.AF.z_mid[0,:]-self.AF.z_mid[1,:])/DEL_T
        
        else:
            # Second-order backwards differencing of body collocation point positions
            self.vx = (3*self.AF.x_mid[0,:]-4*self.AF.x_mid[1,:]+self.AF.x_mid[2,:])/(2*DEL_T) - self.V0
            self.vz = (3*self.AF.z_mid[0,:]-4*self.AF.z_mid[1,:]+self.AF.z_mid[2,:])/(2*DEL_T)
    
    # Calculate pressure distribution on the body
    def pressure(self, RHO, i, DEL_T):
    # Uses Body.(N, x, z, V0, vx, vz, sigma, mu, mu_past), panel_vectors()
    # Gets Body.cp and Body.p
    # Others: dmu_dl, dmu_dt, qpx_tot, qpz_tot
        
        if i > 0:
            
            (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x,self.AF.z)
            
            # Tangential panel velocity dmu/dl, first-order differencing
            dmu_dl = np.empty(self.N)
            dmu_dl[0] = (self.mu[0]-self.mu[1]) / (lpanel[0]/2 + lpanel[1]/2)
            dmu_dl[1:-1] = (self.mu[:-2]-self.mu[2:]) / (lpanel[:-2]/2 + lpanel[1:-1] + lpanel[2:]/2)
            dmu_dl[-1] = (self.mu[-2]-self.mu[-1]) / (lpanel[-2]/2 + lpanel[-1]/2)
            
            # Potential change dmu/dt, second-order differencing after first time step
            if i == 1:
                dmu_dt = (self.mu - self.mu_past[0,:])/DEL_T
            else:
                dmu_dt = (3*self.mu - 4*self.mu_past[0,:] + self.mu_past[1,:])/(2*DEL_T)
            
            # Unsteady pressure calculation (from Matlab code)
            qpx_tot = dmu_dl*tx + self.sigma*nx
            qpz_tot = dmu_dl*tz + self.sigma*nz
            
            self.p = -RHO*(qpx_tot**2 + qpz_tot**2)/2. + RHO*dmu_dt + RHO*(qpx_tot*(self.V0+self.vx) + qpz_tot*self.vz)
            self.cp = self.p / (0.5*RHO*self.V0**2)
           
    # Calculation of drag and lift forces affecting overall airfoil      
    def force(self, i):
    # Uses Body.(x, z, p, N), panel_vectors()
    # Gets Body.(drag, lift)
    # Others: tx, tz, nx, nz, lpanel
        
        if i == 0:
            pass
        else:
            (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x, self.AF.z)
                                  
            Body.drag[i-1] = np.dot(self.p[i-1,:]*lpanel, np.reshape(tx,(self.N,1)))\
                          + np.dot(self.p[i-1,:]*lpanel, np.reshape(-tz,(self.N,1)))
    
            self.lift[i-1] = np.dot(self.p[i-1,:]*lpanel, np.reshape(-nz,(self.N,1)))\
                          + np.dot(self.p[i-1,:]*lpanel, np.reshape(nx,(self.N,1)))