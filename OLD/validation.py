import numpy as np

start = 0.25

x = np.linspace(start,1,20)

R = 0.7
B = 6
tdist = -np.radians(50)*x+np.radians(35) # for x>0.25
pitch = np.radians(46) # x=0.7
cdist = 0.18-0.06*x # for x>0.25

print(np.degrees(tdist+pitch))