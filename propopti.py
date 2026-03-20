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
        results[i,:] = calcloads.calculate_turbine2(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,airfoil) # dT, dQ, a, a_prime, phi, alpha

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

# xi = [0.25,0.5,0.75,1.0]
# tdist = np.radians(sp.interpolate.pchip_interpolate(xi,[33.0872443,26.61575089,20.80975062,8.55404546],x))
# cdist = sp.interpolate.pchip_interpolate(xi,[0.2,0.22,0.15,0.14],x)
# pitch = np.radians(0.03436332) # np.radians(40) # np.radians(0.03436332)

xi = [0.25,0.5,0.7,1.0]
tdist = np.radians(sp.interpolate.pchip_interpolate(xi,[45.53810998, 30.60890109,  0.,         13.92684386],x))
cdist = sp.interpolate.pchip_interpolate(xi,[0.2,
  0.22,        0.18,        0.14],x)
pitch = np.radians(-1.18840555) # np.radians(40) # np.radians(0.03436332)

# xi = [0.25,0.5,0.7,1.0]
# # tdist = np.radians(sp.interpolate.pchip_interpolate(xi,[51.69690749, 49.64579799,  0,          9.61612702],x))
# tdist = np.array([
#     1.33292302, 1.32472074, 1.24249181, 1.19926608, 1.16942446, 1.14603042,
#     1.12632687, 1.10897514, 1.09323882, 1.07867566, 1.06500255, 1.0520291,
#     1.03962218, 1.02768569, 1.01614798, 1.00495593, 0.99406742, 0.98344972,
#     0.97307671, 0.96292733, 0.95298437, 0.94323362, 0.93366323, 0.92426325,
#     0.91502524, 0.90594197, 0.89700725, 0.88821571, 0.87956265, 0.87104399,
#     0.86265609, 0.85439576, 0.84626017, 0.83824675, 0.83028055, 0.82229656,
#     0.81430683, 0.80632256, 0.79835429, 0.79041195, 0.78250493, 0.77464216,
#     0.76683217, 0.75908267, 0.75140255, 0.74379899, 0.73627945, 0.72885124,
#     0.72152155, 0.71429751, 0.70718618, 0.70019461, 0.69332988, 0.68659912,
#     0.68000954, 0.67356846, 0.66728336, 0.66116192, 0.65521201, 0.64944179,
#     0.64383913, 0.63832416, 0.63288755, 0.62752894, 0.62224806, 0.61704475,
#     0.61191895, 0.60687071, 0.60190021, 0.59700779, 0.59219395, 0.58745939,
#     0.582805, 0.57823197, 0.57374175, 0.56933616, 0.56501742, 0.56078827,
#     0.55665203, 0.55261275, 0.54867539, 0.54484604, 0.54113219, 0.53754314,
#     0.53409054, 0.53078912, 0.52765776, 0.52472102, 0.52201143, 0.51957307,
#     0.51746721, 0.51578185, 0.51464858, 0.51427495, 0.51500862, 0.51749352,
#     0.52308851, 0.53536101, 0.56875695
# ])
# cdist = sp.interpolate.pchip_interpolate(xi,[0.14570071, 0.16377962,  0.10663266,  0.09252311],x)
# pitch = np.radians(16.31213639) # np.radians(40) # np.radians(0.03436332)

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
P = Q*RPM/60*2*np.pi

CP = P*2/(rho*np.pi*R**2*U0**3)

print(T)
print(Q)
print(CP)

fig1 = plt.figure(figsize=(8,5))
plt.plot(x[:99],results[:,6])
plt.show()

fig2 = plt.figure(figsize=(8,5))
plt.title(r"Spanwise Distribution of $a$ and $a'$")
plt.plot(x[:99],results[:,2],label=r"$a$",color="black")
plt.plot(x[:99],results[:,3],label=r"$a'$",color="red")
plt.xlabel("$r/R$ [-]")
plt.ylabel("Induction [-]")
plt.legend()
plt.show()

fig = plt.figure(figsize=(8,5))
plt.title(r"Spanwise Distribution of $\alpha$ and $\phi$")
plt.plot(x[:(99)],np.degrees(results[:,5]),label=r"$\alpha$",color="black")
plt.plot(x[:(99)],np.degrees(results[:,4]),label=r"$\phi$",color="red")
plt.xlabel("$r/R$ [-]")
plt.ylabel("Angle ($\degree$)")
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

