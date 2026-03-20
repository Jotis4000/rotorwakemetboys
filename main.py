import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import defgeom
import calcloads
import getforces
import momcorrections

### Geometry Configuration

start = 0.25

el = 50
x = np.linspace(start,1,el)

R = 0.7
B = 6
tdist = -np.radians(50)*x+np.radians(35) # for x>0.25
pitch = np.radians(46) # x=0.7
cdist = 0.18-0.06*x # for x>0.25

# xi = [0.25,0.5,0.75,1.0]
# tdist = np.radians(sp.interpolate.pchip_interpolate(xi,[33.0872443,26.61575089,20.80975062,8.55404546],x))
# cdist = sp.interpolate.pchip_interpolate(xi,[0.2,0.22,0.15,0.14],x)
# pitch = np.radians(30)

# xi = [0.25,0.5,0.7,1.0]
# tdist = np.radians(sp.interpolate.pchip_interpolate(xi,[10.90035694,3.7214586,0,-12.62668013],x))
# cdist = sp.interpolate.pchip_interpolate(xi,[0.2,0.22,0.18,0.14],x)
# pitch = np.radians(57.10674679) # 21.10674679

airfoil = "data/ARAD8pct_polar.txt"

### Simulation Configuration

U0 = 60
h = 2000
rho = 1.00649
fi = 0
ry = 0

J = np.array([1.6,2.0,2.4]) #1.95
RPM = U0/(J*2*R)*U0
print("Advance Ratio J: "+str(J))
print("RPM: "+str(RPM))

maxiter = 101
tol = 1e-6

### Calculate geometry specification

c, theta, sigma = defgeom.defGeom(R,B,x,start,tdist,pitch,cdist)

### Main Iteration Loop

results = np.zeros([3,len(x)-1,7])

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

T = np.zeros(len(J))
Q = np.zeros(len(J))
for i in range(len(J)):
    tempT=0
    tempQ=0
    for j in range(len(results[i,:,1])):
        tempT+=results[i,j,0]*(x[j+1]-x[j])*R
        T[i]=tempT
        tempQ+=results[i,j,1]*(x[j+1]-x[j])*R
        Q[i]=tempQ

P = Q*RPM/60*2*np.pi
CP = P*2/(rho*np.pi*R**2*U0**3)
print("Thrust: "+str(T))
print("Torque: "+str(Q))
print("Power: "+str(P))
print("CP: "+str(CP))
print("Eta: "+str(T*U0/P))

# print(results[1,1,6])
# plt.plot(x[:99],np.degrees(results[0,:,5]),label=r"$a$",color="black")

### ALPHA

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $\alpha$")
plt.plot(x[:(el-1)],np.degrees(results[0,:,5]),label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[:(el-1)],np.degrees(results[1,:,5]),label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[:(el-1)],np.degrees(results[2,:,5]),label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Angle of Attack $\alpha$ [$\degree$]")
plt.legend()
plt.show()

### PHI

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $\phi$")
plt.plot(x[:(el-1)],np.degrees(results[0,:,4]),label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[:(el-1)],np.degrees(results[1,:,4]),label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[:(el-1)],np.degrees(results[2,:,4]),label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Inflow Angle $\phi$ [$\degree$]")
plt.legend()
plt.show()

### a and a'

fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(18,6))
plt.suptitle(r"Spanwise Distribution of $a$ and $a'$")
ax1.plot(x[:(el-1)],results[0,:,2],label=r"Axial $a$",color="black")
ax1.plot(x[:(el-1)],results[0,:,3],label=r"Azim. $a'$",color="red")
ax1.set_title("J="+str(J[0]))
ax2.plot(x[:(el-1)],results[1,:,2],label=r"Axial $a$",color="black")
ax2.plot(x[:(el-1)],results[1,:,3],label=r"Azim. $a'$",color="red")
ax2.set_title("J="+str(J[1]))
ax3.plot(x[:(el-1)],results[2,:,2],label=r"Axial $a$",color="black")
ax3.plot(x[:(el-1)],results[2,:,3],label=r"Azim. $a'$",color="red")
ax3.set_title("J="+str(J[2]))
ax1.set_xlabel("$r/R$ [-]")
ax1.set_ylabel("Induction [-]")
ax2.set_xlabel("$r/R$ [-]")
ax2.set_ylabel("Induction [-]")
ax3.set_xlabel("$r/R$ [-]")
ax3.set_ylabel("Induction [-]")
ax1.legend()
ax2.legend()
ax3.legend()
plt.tight_layout()
plt.show()

# fig2 = plt.figure(figsize=(8,5))
# plt.title(r"Spanwise Distribution of $a$ and $a'$")
# # plt.plot(x[:(el-1)],results[0,:,2],label=r"$a$",color="black")
# # plt.plot(x[:(el-1)],results[0,:,3],label=r"$a'$",color="red")
# # plt.plot(x[:(el-1)],results[1,:,2],label=r"$a$",color="green")
# # plt.plot(x[:(el-1)],results[1,:,3],label=r"$a'$",color="green")
# plt.plot(x[:(el-1)],results[2,:,2],label=r"$a$",color="blue")
# plt.plot(x[:(el-1)],results[2,:,3],label=r"$a'$",color="orange")
# plt.xlabel("$r/R$ [-]")
# plt.ylabel("Induction [-]")
# plt.legend()
# plt.show()

# plt.plot(x[:(el-1)],results[2,:,6],label=r"$a$",color="blue")
# plt.show()

### THRUST DIST.

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $dT$")
plt.plot(x[:(el-1)],results[0,:,0],label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[:(el-1)],results[1,:,0],label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[:(el-1)],results[2,:,0],label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Thrust Distribution [$N/m$]")
plt.legend()
plt.show()

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $dQ$")
plt.plot(x[:(el-1)],results[0,:,1],label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[:(el-1)],results[1,:,1],label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[:(el-1)],results[2,:,1],label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Torque Distribution [$Nm/m$]")
plt.legend()
plt.show()

# fig3 = plt.figure(figsize=(8,5))
# plt.title(r"Spanwise Distribution of $dT$ and $dQ$")
# plt.plot(x[:(el-1)],results[0,:,0],label=r"$dT$",color="orange")
# plt.plot(x[:(el-1)],results[0,:,1],label=r"$dQ$",color="green")
# plt.xlabel("$r/R$ [-]")
# plt.ylabel("Force or Moment [N or Nm]") ## THIS IS RETARDED
# plt.legend()
# plt.show()

plt.plot(J,T,label="Thrust")
plt.plot(J,Q,label="Torque")
plt.plot(J,P/1000,label="Power")
plt.plot()
plt.xlabel("Advance Ratio J [-]")
plt.ylabel("Force or Moment [N or N]") ## THIS IS RETARDED FIND BETTER WAY TO PLOT
plt.legend()
plt.show()