import numpy as np
import scipy as sp
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
        results[i, :] = calcloads.calculate_turbine2(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,airfoil)

    # Calc power
    dr = np.diff(x * R)
    omega = RPM / 60 * 2 * np.pi
    Q = np.sum(results[:, 1] * dr)
    P = omega * Q

    # # Penalties for invalid geometry or NaNs
    # if np.any(c <= 0) or np.any(sigma <= 0) or np.isnan(P) or np.isinf(P):
    #     return 1e9

    return -P

def objectivefunctionQuad(params, R, B, start, U0, RPM, J, x, airfoil):

    v1, v2, v3, v4, v5, v6, v7  = params

    # Design Variables
    tdist = np.radians(v1) * x**2 + np.radians(v2) *x +np.radians(v3)
    pitch = np.radians(v4)  
    cdist = v7 + v6 * x + v5 * x**2
    

    # Define Geometry
    c, theta, sigma = defgeom.defGeom(R, B, x, start, tdist, pitch, cdist)

    n = len(x) - 1
    results = np.zeros((n, 7))

    # Run BEM
    for i in range(n):
        r_local = x[i] * R
        results[i, :] = calcloads.calculate_turbine2(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,airfoil)

    # Calc power
    dr = np.diff(x * R)
    omega = RPM / 60 * 2 * np.pi
    Q = np.sum(results[:, 1] * dr)
    P = omega * Q

    # # Penalties for invalid geometry or NaNs
    # if np.any(c <= 0) or np.any(sigma <= 0) or np.isnan(P) or np.isinf(P):
    #     return 1e9

    return -P

def objectivefunctionQuadNew(params, R, B, start, U0, RPM, J, x, airfoil):

    v1, v2, v3, v4, v5, v6, v7  = params

    xi = [0.25,0.4,1.0]
    # Design Variables
    tcoeff = np.polyfit(xi,[v1,v2,v3],2)
    ccoeff = np.polyfit(xi,[v5,v6,v7],2)

    tdist = np.radians(tcoeff[0]*x**2+tcoeff[1]*x+tcoeff[2])
    cdist = ccoeff[0]*x**2+ccoeff[1]*x+ccoeff[2]
    # tdist = np.radians(v1) * x**2 + np.radians(v2) *x +np.radians(v3)
    pitch = np.radians(v4)  
    # cdist = v7 + v6 * x + v5 * x**2
    
    # Define Geometry
    c, theta, sigma = defgeom.defGeom(R, B, x, start, tdist, pitch, cdist)

    n = len(x) - 1
    results = np.zeros((n, 7))

    # Run BEM
    for i in range(n):
        r_local = x[i] * R
        results[i, :] = calcloads.calculate_turbine2(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,airfoil)

    # Calc power
    dr = np.diff(x * R)
    omega = RPM / 60 * 2 * np.pi
    Q = np.sum(results[:, 1] * dr)
    P = omega * Q

    # # Penalties for invalid geometry or NaNs
    # if np.any(c <= 0) or np.any(sigma <= 0) or np.isnan(P) or np.isinf(P):
    #     return 1e9

    return -P

def objectivefunctionBez(params, R, B, start, U0, RPM, J, x, airfoil, xi):

    v1, v2, v3, v4, v5, v6, v7, v8, v9  = params

    # Design Variables
    tdist = np.radians(sp.interpolate.pchip_interpolate(xi,[v1,v2,v3,v4],x))
    cdist = sp.interpolate.pchip_interpolate(xi,[v6,v7,v8,v9],x)

    pitch = np.radians(v5)
    
    # Define Geometry
    c, theta, sigma = defgeom.defGeom(R, B, x, start, tdist, pitch, cdist)

    n = len(x) - 1
    results = np.zeros((n, 7))

    # Run BEM
    for i in range(n):
        r_local = x[i] * R
        # results[i, :] = calcloads.calculate_turbine2(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,airfoil)
        results[i, :] = calcloads.calculate_element_loads3(x[i]*R,R,start*R,c[i],theta[i],U0,RPM/60*2*np.pi,sigma[i],B,J,airfoil)

    # Calc power
    dr = np.diff(x * R)
    omega = RPM / 60 * 2 * np.pi
    Q = np.sum(results[:, 1] * dr)
    P = omega * Q

    # # Penalties for invalid geometry or NaNs
    if (np.degrees(results[:,5])<-10).any():
        # print("AH")
        return 1e9

    # print(str(B*c/(2*np.pi*R)))

    return P