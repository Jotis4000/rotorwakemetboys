from scipy.optimize import minimize
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import defgeom
import calcloads
import getforces
import momcorrections
from objectivefunction import objectivefunction

### Geometry Configuration

start = 0.25

x = np.linspace(start,1,100)

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

res = minimize(
    objectivefunction, 
    initial_guess, 
    args= (R, B, start, U0, RPM, J, x, airfoil), 
    method='powell',
    bounds=bounds,
    callback=monitor_progress,
    options={'disp': True}
)

print(f"Optimal variables: {res.x}")
print(f"Maximum PowerCoefficient: {-res.fun*2/rho*np.pi*R**2*U0**3}")
