import numpy as np
import scipy as sp
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
h = 2000
rho = 1.00649
fi = 0
ry = 0

J = np.array([1.6,2.0,2.4])
RPM = U0/(J*2*R)*60
print("Advance Ratio J: "+str(J))
print("RPM: "+str(RPM))

maxiter = 101
tol = 1e-6

### Calculate geometry specification

c, theta, sigma = defgeom.defGeom(R,B,x,start,tdist,pitch,cdist)
# Vax,Vtgt,Veff,phi = defgeom.calcPhi(x,R,RPM,U0)

### Obtain forces



### Main Iteration Loop

results = np.zeros([3,len(x)-1,7])
# results2 = np.zeros([len(x)-1,7])
# results3 = np.zeros([len(x)-1,7])

for j in range(len(J)):
    for i in range(len(x)-1):

        results[j,i,:] = calcloads.calculate_element_loads3(x[i]*R,R,start*R,c[i],theta[i],U0,RPM[j]/60*2*np.pi,sigma[i],B,J[j],airfoil) # dT, dQ, a, a_prime, phi, alpha

plt.rcParams.update({
    "text.usetex": False, # Set to True if you have LaTeX installed on your system
    "font.family": "serif",
    "font.size": 12,
    "axes.grid": True,
    "grid.linestyle": '--',
    "grid.alpha": 0.7
})

# print(results[1,1,6])
# plt.plot(x[:99],np.degrees(results[0,:,5]),label=r"$a$",color="black")

fig = plt.figure(figsize=(8,5))
plt.title(r"Spanwise Distribution of $\alpha$ and $\phi$")
plt.plot(x[:99],np.degrees(results[0,:,5]),label=r"$\alpha$",color="black")
plt.plot(x[:99],np.degrees(results[0,:,4]),label=r"$\phi$",color="red")
plt.xlabel("$r/R$ [-]")
plt.ylabel("Angle ($\degree$)")
plt.legend()
plt.show()

fig2 = plt.figure(figsize=(8,5))
plt.title(r"Spanwise Distribution of $a$ and $a'$")
plt.plot(x[:99],results[0,:,2],label=r"$a$",color="black")
plt.plot(x[:99],results[0,:,3],label=r"$a'$",color="red")
plt.xlabel("$r/R$ [-]")
plt.ylabel("Induction [-]")
plt.legend()
plt.show()

fig3 = plt.figure(figsize=(8,5))
plt.title(r"Spanwise Distribution of $dT$ and $dQ$")
plt.plot(x[:99],results[0,:,0],label=r"$dT$",color="orange")
plt.plot(x[:99],results[0,:,1],label=r"$dQ$",color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel("Force or Moment [N or Nm]") ## THIS IS RETARDED
plt.legend()
plt.show()

T = np.zeros(len(J))
Q = np.zeros(len(J))
# print(len(J))
for i in range(len(J)):
    tempT=0
    tempQ=0
    for j in range(len(results[i,:,1])):
        tempT+=results[i,j,0]*(x[j+1]-x[j])*R
        T[i]=tempT
        tempQ+=results[i,j,1]*(x[j+1]-x[j])*R
        Q[i]=tempQ


    # T[i] = sp.integrate.quad(x*R,results[i,:,1])
    
    # sum(results[i,:,0])
    # Q[i] = sum(results[i,:,1])

    # print(results[i,:,0])
    # print(i)

plt.plot(J,T,label="Thrust")
plt.plot(J,Q,label="Torque")
plt.xlabel("Advance Ratio J [-]")
plt.ylabel("Force or Moment [N or N]") ## THIS IS RETARDED FIND BETTER WAY TO PLOT
plt.legend()
plt.show()