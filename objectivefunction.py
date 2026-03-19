import numpy as np
import defgeom
import calcloads
import getforces
import momcorrections

def calculate_CP(results, R, RPM, x):
    
    omega = RPM / 60 * 2 * np.pi
    total_Q = 0

    dr = np.diff(x) * R
    local_Q = results[:-1, 1] * dr
    total_Q = np.sum(local_Q)
    
    return total_Q * omega

def objectivefunction(params, R, B, start, U0, RPM, J, x, airfoil):

    v1, v2, v3, v4, v5 = params

    # Design Variables
    tdist = -np.radians(v1) * x + np.radians(v2) 
    pitch = np.radians(v3)  
    cdist = v4 - v5 * x  

    # Define Geometry
    c, theta, sigma = defgeom.defGeom(R, B, x, start, tdist, pitch, cdist)

    n = len(x) - 1
    results = np.zeros((n, 7))

    # Run BEM
    for i in range(n):
        r_local = x[i] * R
        results[i, :] = calcloads.calculate_element_loads3(r_local,R,start * R,c[i],theta[i],U0,RPM / 60 * 2 * np.pi,sigma[i],B,J,airfoil)

    # Calc power
    dr = np.diff(x * R)
    omega = RPM / 60 * 2 * np.pi
    Q = np.sum(results[:, 1] * dr)
    P = omega * Q

    # # Penalties for invalid geometry or NaNs
    # if np.any(c <= 0) or np.any(sigma <= 0) or np.isnan(P) or np.isinf(P):
    #     return 1e9

    return -P


