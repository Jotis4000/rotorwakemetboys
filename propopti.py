import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import defgeom
import calcloads

def getData(x,R,B,airfoil,tdist,pitch,cdist,U0,rho,RPM):
    
    ### Calculate geometry specification

    c, theta, sigma = defgeom.defGeom(R,B,x,start,tdist,pitch,cdist)

    results = np.zeros((len(x)-1,7))
    for i in range(len(x)-1):
        results[i,:] = calcloads.calculate_turbine(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,airfoil) # dT, dQ, a, a_prime, phi, alpha

    return results

### Geometry Configuration

start = 0.25
x = np.linspace(start,1,100)
R = 0.7
B = 6
airfoil = "data/ARAD8pct_polar.txt"

tdist = -np.radians(50)*x+np.radians(35) # for x>0.25
pitch = np.radians(46) # x=0.7
cdist = 0.18-0.06*x # for x>0.25

### Simulation Configuration

U0 = 60
h = 2000
rho = 1.00649
RPM = 1200

maxiter = 101
tol = 1e-6

results = getData(x,R,B,airfoil,tdist,pitch,cdist,U0,rho,RPM)

tempT=0
tempQ=0
for j in range(len(results[:,1])):
    tempT+=results[j,0]*(x[j+1]-x[j])
    tempQ+=results[j,1]*(x[j+1]-x[j])*R

T=tempT
Q=tempQ

print(T)
print(Q)

fig2 = plt.figure(figsize=(8,5))
plt.title(r"Spanwise Distribution of $a$ and $a'$")
plt.plot(x[:99],results[:,2],label=r"$a$",color="black")
plt.plot(x[:99],results[:,3],label=r"$a'$",color="red")
plt.xlabel("$r/R$ [-]")
plt.ylabel("Induction [-]")
plt.legend()
plt.show()

fig3 = plt.figure(figsize=(8,5))
plt.title(r"Spanwise Distribution of $dT$ and $dQ$")
plt.plot(x[:99],results[:,0],label=r"$dT$",color="orange")
plt.plot(x[:99],results[:,1],label=r"$dQ$",color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel("Force or Moment [N or Nm]") ## THIS IS RETARDED
plt.legend()
plt.show()

### Main Iteration Loop

# results = np.zeros(len(x)-1,7)
# results2 = np.zeros([len(x)-1,7])
# results3 = np.zeros([len(x)-1,7])

