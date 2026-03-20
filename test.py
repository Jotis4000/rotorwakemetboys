import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import calcloads
import defgeom

x = np.linspace(0.25,1,50)
xi = [0.25,0.5,0.75,1.0]
yi = [10,25,18,8]
airfoil = "data/ARAD8pct_polar.txt"

# test = sp.interpolate.pchip_interpolate(xi,yi,x)
# plt.plot(x,test)
# plt.show()

r = x*0.7
B = 6
Cl = 0.4
U0 = 60
RPM = 1200
R=0.7
start=0.25

tdist = np.radians(-32.43*x+20.01)
pitch = np.radians(23.5)
cdist = -8/75*x+.246666666666667

tdist = np.radians(sp.interpolate.pchip_interpolate(xi,[33.0872443,26.61575089,20.80975062,8.55404546],x))
cdist = sp.interpolate.pchip_interpolate(xi,[0.2,0.22,0.15,0.14],x)
pitch = np.radians(0.03436332)

plt.plot(x,tdist+pitch) #+pitch
plt.show()
plt.plot(x,cdist)
plt.show()

c, theta, sigma = defgeom.defGeom(R, B, x, start, tdist, pitch, cdist)

results = np.zeros((len(x)-1,7))
for i in range(len(x)-1):
    results[i,:] = calcloads.calculate_turbine2(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,airfoil) # dT, dQ, a, a_prime, phi, alpha]

print(results[:,1])

plt.plot(x[:49],results[:,1])
plt.show()