import numpy as np

def calcPres(r,pinf,rho,U0,dT,a,dQ):

    p1 = pinf+0.5*rho*U0**2
    p2 = p1

    dr  = np.zeros(len(dT))
    # dTr = np.zeros(len(dT))
    dpr  = np.zeros(len(dT))
    for i in range(len(r)-1):

        dr[i] = r[i+1]-r[i]
        dpr[i] = 6*dT[i]/(np.pi*(r[i+1]-r[i])**2)

    # dpr = 6*dT/(2*np.pi*r[:199])
    # dpr = 2*rho*U0**2*a*(1+a)

    p3 = p2+dpr
    print(dpr)

    p4 = 0

    p3int=0
    for i in range(len(dpr)):
        p3int+=dpr[i]*dr[i]

    # dTcheck = 
    print(p3int/(np.pi*0.7**2))

    return p1,p2,p3,p4