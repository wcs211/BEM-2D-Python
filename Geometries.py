import numpy as np
    
# VandeVooren airfoil type
def VanDeVooren(Body):
# uses Body.(c, epsilon, k, N)
# gets Body.(x, z, x_col, z_col)
# others: a, phi, r1, r2, phi1, phi2
    
    a=Body.c*((1+Body.epsilon)**(Body.k-1))*(2**(-Body.k))
    
    phi = np.linspace(0,np.pi,Body.N/2+1)
    
    r1= np.sqrt((a*np.cos(phi)-a)**2+(a**2)*np.sin(phi)**2)
    r2= np.sqrt((a*np.cos(phi)-Body.epsilon*a)**2+(a**2)*np.sin(phi)**2)
    
    phi1=np.arctan2((a*np.sin(phi)) , (a*np.cos(phi) -a))
    phi2=np.arctan2(a*np.sin(phi) ,(a*np.cos(phi) -Body.epsilon*a))
    
    Body.x=((r1**Body.k)/(r2**(Body.k-1)))*(np.cos(Body.k*phi1)*np.cos((Body.k-1)*phi2) + np.sin(Body.k*phi1)*np.sin((Body.k-1)*phi2))
    Body.z_top=((r1**Body.k)/(r2**(Body.k-1)))*(np.sin(Body.k*phi1)*np.cos((Body.k-1)*phi2) - np.cos(Body.k*phi1)*np.sin((Body.k-1)*phi2))
    Body.z_bot=-((r1**Body.k)/(r2**(Body.k-1)))*(np.sin(Body.k*phi1)*np.cos((Body.k-1)*phi2) - np.cos(Body.k*phi1)*np.sin((Body.k-1)*phi2))
    
    Body.x=Body.x-Body.x[-1]           #carrying the leading edge to the origin
    Body.x[0]=Body.c
    
    Body.z_top[0]=0
    Body.z_bot[0]=0
    Body.z_bot[-1]=0
    
    # merge top and bottom surfaces together
    Body.x=np.hstack((Body.x , Body.x[-2::-1]))
    Body.z=np.hstack((Body.z_bot , Body.z_top[-2::-1]))
    
    Body.x_col=((Body.x[1:] + Body.x[:-1])/2)
    Body.z_col=((Body.z[1:] + Body.z[:-1])/2)

#flat plate geometry
def FlatPlate(Body):
    
    Body.x_f=np.linspace(1,0,Body.N/2+1)
    Body.z_f=np.zeros(Body.N/2+1)