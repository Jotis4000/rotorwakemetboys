import numpy as np
import matplotlib.pyplot as plt
import defgeom
import calcloads
import getforces
import momcorrections

### Geometry Configuration

start = 0.25

x = np.linspace(start,1,100)

R = 0.7
B = 6
tdist = -np.radians(50)*x+np.radians(35) # for x>0.25
pitch = np.radians(46) # x=0.7
cdist = 0.18-0.06*x # for x>0.25

airfoil = "data/ARAD8pct_polar.txt"

### Simulation Configuration

U0 = 60
RPM = 1200
h = 2000
rho = 1.00649
fi = 0
ry = 0

maxiter = 101
tol = 1e-6

### Calculate geometry specification

c, theta, sigma = defgeom.defGeom(R,B,x,start,tdist,pitch,cdist)
Vax,Vtgt,Veff,phi = defgeom.calcPhi(x,R,RPM,U0)

### Obtain forces



### Main Iteration Loop

results =np.zeros([len(x)-1,7])

for i in range(len(x)-1):

    results[i,:] = calcloads.calculate_element_loads2(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,phi[i],airfoil) # dT, dQ, a, a_prime, phi, alpha


plt.plot(x[:99],results[:,2],label="a")
plt.plot(x[:99],results[:,3],label="a'")
# plt.plot(x[:99],results[:,6],label="F")
plt.legend()
plt.show()