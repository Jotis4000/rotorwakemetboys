import numpy as np

def defGeom(R, B, x, start, tdist, pitch, cdist):

    r = x*R

    c = R*cdist                 # Local Chord
    theta = tdist+pitch         # Local Pitch Angle
    sigma = B*c/(2*np.pi*r)     # Local Solidity

    return c, theta, sigma

def calcPhi(x,R,RPM,U0):

    Vax = U0
    Vtgt = x*(2*np.pi*R)*RPM/60
    Veff = np.sqrt(Vax**2+Vtgt**2)
    phi = np.arctan(Vax/Vtgt)

    return Vax,Vtgt,Veff,phi