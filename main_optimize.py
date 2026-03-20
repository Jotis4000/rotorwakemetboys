from scipy.optimize import minimize,dual_annealing,differential_evolution
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import defgeom
import calcloads
import getforces
import momcorrections
from objectivefunction import objectivefunction,objectivefunctionQuad,objectivefunctionQuadNew,objectivefunctionBez

### Geometry Configuration

start = 0.25

el = 50
x = np.linspace(start,1,el)

R = 0.7
B = 6
airfoil = "data/ARAD8pct_polar.txt"

### Simulation Configuration

U0 = 60
h = 2000
rho = 1.00649
fi = 0
ry = 0

J = 2.4
RPM = U0/(J*2*R)*60

tol = 1e-4

### Optimization

initial_guess = [50, 35, 46, 0.18, 0.06]

def monitor_progress(xk):
    print(f"Iteration complete. Current variables: {xk}")

bounds = [
    (10, 80),    
    (5, 45),    
    (0, 90),    
    (0.1, 0.5),  
    (0.01, 0.1) 
]

boundsQuad = [
    (-2,2),
    (-2,2),    
    (0, 1),    
    (0, 45),
    (-2,-1),
    (0,2),    
    (0.1, 0.5)  
]

boundsQuadNew = [
    (20,80),
    (-10,40),    
    (-20,20),    
    (-10,20),
    (0.1,0.18),
    (0.05,0.2),    
    (0.01, 0.08)  
]

boundsBez = [
    (20,60),
    (10,40),
    (0,0),    
    (-20,20),    
    (-20,40),
    (0.1,0.2),
    (0.1,0.22),
    (0.05,0.18),    
    (0.01, 0.14)  
]

# boundsBez = [
#     (0,80),
#     (0,60),
#     (0,0),    
#     (-20,20),    
#     (-20,40),
#     (0.1,0.2),
#     (0.1,0.22),
#     (0.05,0.18),    
#     (0.01, 0.14)  
# ]

# res = minimize(
#     objectivefunction, 
#     initial_guess, 
#     args= (R, B, start, U0, RPM, J, x, airfoil), 
#     method='Nelder-Mead',
#     bounds=bounds,
#     callback=monitor_progress,
#     options={'disp': True}
# )

xi = [0.25,0.5,0.7,1.0]

res = differential_evolution(objectivefunctionBez,boundsBez,args=(R, B, start, U0, RPM, J, x, airfoil, xi),strategy='best1bin', 
                                maxiter=5, popsize=15, polish=True, disp=True)

print(f"Optimal variables: {res.x}")
print(f"Maximum PowerCoefficient: {-res.fun*2/(rho*np.pi*R**2*U0**3)}")

twist = sp.interpolate.pchip_interpolate(xi,[res.x[0]+res.x[4],res.x[1]+res.x[4],res.x[2]+res.x[4],res.x[3]+res.x[4]],x)
chord = sp.interpolate.pchip_interpolate(xi,[res.x[5],res.x[6],res.x[7],res.x[8]],x)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(x,twist)
ax1.set_title("Blade Twist Distribution")
ax1.set_xlabel("Radial Position (r/R)")
ax1.set_ylabel("Twist Angle (deg)")
ax1.grid(True)
ax1.legend()

ax2.plot(x,chord)
ax2.set_title("Blade Chord Distribution")
ax2.set_xlabel("Radial Position (r/R)")
ax2.set_ylabel("Chord Length (m)")
ax2.grid(True)
ax2.legend()

# Plot Twist on the first axis (left)
# ax1.plot(x, twist, label="Twist", color='blue')
# ax1.set_title("Blade Twist Distribution")
# ax1.set_xlabel("Radial Position (r/R)")
# ax1.set_ylabel("Twist Angle (deg)")
# ax1.grid(True)
# ax1.legend()

# # Plot Chord on the second axis (right)
# ax2.plot(x, chord, label="Chord", color='green')
# ax2.set_title("Blade Chord Distribution")
# ax2.set_xlabel("Radial Position (r/R)")
# ax2.set_ylabel("Chord Length (m)")
# ax2.grid(True)
# ax2.legend()

# Adjust layout so labels don't overlap
plt.tight_layout()
plt.show()