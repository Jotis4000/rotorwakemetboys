import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import defgeom
import calcloads
import momcorrections

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.size": 12,
    "axes.grid": True,
    "grid.linestyle": '--',
    "grid.alpha": 0.7
})

### Geometry Configuration

start = 0.25

el = 200
x = np.linspace(start,1,el)

R = 0.7
B = 6
tdist = -np.radians(50)*x+np.radians(35) # for x>0.25
pitch = np.radians(46) # 29
cdist = 0.18-0.06*x # for x>0.25

xi = [0.25,0.5,0.7,1.0]
tdist2 = np.radians(sp.interpolate.pchip_interpolate(xi,[20.8362307,   10.03308008,   0.,         -10.06593984],x))
cdist2 = sp.interpolate.pchip_interpolate(xi,[0.16263473,   0.18026756,   0.17807533,   0.13926443],x)
pitch2 = np.radians(30.85819216) # +13.4

airfoil = "data/ARAD8pct_polar.txt"

### Simulation Configuration

U0 = 60
h = 2000
rho = 1.00649
pinf = 79495

# J = np.array([1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4]) 
J = np.array([1.6,2.0,2.4]) # 2.13776722
RPM = U0/(J*2*R)*U0
print("Advance Ratio J: "+str(J))
print("RPM: "+str(RPM))

### Calculate geometry specification

c, theta, sigma = defgeom.defGeom(R,B,x,start,tdist,pitch,cdist)
c2, theta2, sigma2 = defgeom.defGeom(R,B,x,start,tdist2,pitch2,cdist2)

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))

ax1.plot(x,np.degrees(theta),color="black",label="Original")
ax1.plot(x,np.degrees(tdist2)+np.degrees(pitch2),color="red",label="Optimized")
ax1.set_title("Blade Twist Distribution")
ax1.set_xlabel("Radial Position [$r/R$]")
ax1.set_ylabel("Twist Angle [$\degree$]")
ax1.grid(True)
ax1.legend()

ax2.plot(x,cdist,color="black",label="Original")
ax2.plot(x,cdist2,color="red",label="Optimized")
ax2.set_title("Blade Chord Distribution")
ax2.set_xlabel("Radial Position [$r/R$]")
ax2.set_ylabel("Chord Normalized Length [$-$]")
ax2.set_ylim([0,0.2])
ax2.grid(True)
ax2.legend()
plt.tight_layout()

### Main Iteration Loop

results = np.zeros([len(J),len(x)-1,8])
results2 = np.zeros([len(J),len(x)-1,8])

for j in range(len(J)):
    for i in range(len(x)-1):

        results[j,i,:] = calcloads.calculate_element_loads3(x[i]*R,R,start*R,c[i],theta[i],U0,RPM[j]/60*2*np.pi,sigma[i],B,J[j],airfoil) # dT, dQ, a, a_prime, phi, alpha
        results2[j,i,:] = calcloads.calculate_element_loads3(x[i]*R,R,start*R,c2[i],theta2[i],U0,RPM[j]/60*2*np.pi,sigma2[i],B,J[j],airfoil) # dT, dQ, a, a_prime, phi, alpha

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

fig, ax1 = plt.subplots(figsize=(8,5))

ax1.set_xlabel('$r/R$ [-]')
ax1.set_ylabel(r'Angle of Attack $\alpha$ [-]', color="red")
ax1.plot(x[1:(el-1)],np.degrees(results[1,:,5][     1:]), color="red")
ax1.tick_params(axis='y', labelcolor="red")

ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

ax2.set_ylabel(r'L/D [-]', color="blue")  # we already handled the x-label with ax1
ax2.plot(x[1:(el-1)], results[1,:,7][1:], color="blue")
ax2.tick_params(axis='y', labelcolor="blue")

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

### ALPHA

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $\alpha$")
plt.plot(x[1:(el-1)],np.degrees(results[0,:,5][1:]),label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[1:(el-1)],np.degrees(results[1,:,5][1:]),label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[1:(el-1)],np.degrees(results[2,:,5][1:]),label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Angle of Attack $\alpha$ [$\degree$]")
plt.legend()

### PHI

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $\phi$")
plt.plot(x[1:(el-1)],np.degrees(results[0,:,4][1:]),label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[1:(el-1)],np.degrees(results[1,:,4][1:]),label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[1:(el-1)],np.degrees(results[2,:,4][1:]),label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Inflow Angle $\phi$ [$\degree$]")
plt.legend()

### a and a'

fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(18,6))
plt.suptitle(r"Spanwise Distribution of $a$ and $a'$")
ax1.plot(x[1:(el-1)],results[0,:,2][1:],label=r"Axial $a$",color="black")
ax1.plot(x[1:(el-1)],results[0,:,3][1:],label=r"Azim. $a'$",color="red")
ax1.set_title("J="+str(J[0]))
ax2.plot(x[1:(el-1)],results[1,:,2][1:],label=r"Axial $a$",color="black")
ax2.plot(x[1:(el-1)],results[1,:,3][1:],label=r"Azim. $a'$",color="red")
ax2.set_title("J="+str(J[1]))
ax3.plot(x[1:(el-1)],results[2,:,2][1:],label=r"Axial $a$",color="black")
ax3.plot(x[1:(el-1)],results[2,:,3][1:],label=r"Azim. $a'$",color="red")
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

### THRUST DIST.

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $dT$")
plt.plot(x[1:(el-1)],results[0,:,0][1:],label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[1:(el-1)],results[1,:,0][1:],label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[1:(el-1)],results[2,:,0][1:],label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Thrust Distribution [$N/m$]")
plt.legend()
plt.tight_layout()

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $dQ$")
plt.plot(x[1:(el-1)],results[0,:,1][1:],label=r"$J=$"+str(J[0]),color="black")
plt.plot(x[1:(el-1)],results[1,:,1][1:],label=r"$J=$"+str(J[1]),color="red")
plt.plot(x[1:(el-1)],results[2,:,1][1:],label=r"$J=$"+str(J[2]),color="green")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Torque Distribution [$Nm/m$]")
plt.legend()
plt.tight_layout()
plt.show()

# Final Plots

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $dT$")
plt.plot(x[2:(el-1)],results[1,:,0][2:],label=r"Original",color="red")
plt.plot(x[2:(el-1)],results2[1,:,0][2:],label=r"Optimized",color="darkred")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Thrust Distribution [$N/m$]")
plt.legend()
plt.tight_layout()

fig1 = plt.figure(figsize=(10,6))
plt.title(r"Spanwise Distribution of $dQ$")
plt.plot(x[2:(el-1)],results[1,:,1][2:],label=r"Original",color="red")
plt.plot(x[2:(el-1)],results2[1,:,1][2:],label=r"Optimized",color="darkred")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Torque Distribution [$Nm/m$]")
plt.legend()
plt.tight_layout()
plt.show()

fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(14,5))

ax1.set_xlabel('Advance Ratio J [-]')
ax1.set_ylabel(r'Thrust [$N$]')
ax1.plot(J,T, color="darkblue", label="Thrust")
ax1.tick_params(axis='y')
ax1.legend()

ax2.set_xlabel('Advance Ratio J [-]')
ax2.set_ylabel(r'Torque [$Nm$]')
ax2.plot(J,Q, color="red",label="Torque")
ax2.tick_params(axis='y')
ax2.legend()

ax3.set_xlabel('Advance Ratio J [-]')
ax3.set_ylabel(r'Power [$kW$]')
ax3.plot(J,P/1000, color="orange",label="Power")
ax3.tick_params(axis='y')
ax3.legend()
fig.tight_layout()
plt.show()