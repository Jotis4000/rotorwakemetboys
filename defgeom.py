import numpy as np

# x = np.linspace(0.25,1,100)

# R = 0.7
# B = 6
# start = 0.25
# twist = -np.radians(50)*x+np.radians(35) # for x>0.25
# pitch = np.radians(46) # x=0.7
# cdist = 0.18-0.06*x # for x>0.25

def defGeom(R, B, x, start, tdist, pitch, cdist):

    r = x*R

    c = R*cdist                 # Local Chord
    theta = tdist+pitch         # Local Pitch Angle
    sigma = B*c/(2*np.pi*r)     # Local Solidity

    return c, theta, sigma

# c, theta, sigma = defGeom(R,B,x,start,twist,pitch,cdist)

# import matplotlib.pyplot as plt
# plt.plot(x,c)
# plt.show()
# plt.plot(x,np.degrees(theta))
# plt.show()
# plt.plot(x,sigma)
# plt.show()