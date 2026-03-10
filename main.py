import numpy as np
import matplotlib.pyplot as plt

### Geometry Configuration

r = np.linspace(0,1,100)

R = 0.7
B = 6
start = 0.25
twist = -np.radians(50)*r+np.radians(35) # for r>0.25
pitch = np.radians(46) # r=0.7
cdist = 0.18-0.06*r # for r>0.25

airfoil = "data/ARAD8pct_polar.txt"

### Simulation Configuration

U0 = 60
RPM = 1200
h = 2000
fi = 0
ry = 0