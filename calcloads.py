import numpy as np

def load_airfoil_data(path):

    data = np.loadtxt(path, skiprows=2)
    
    alpha_polar_deg = data[:, 0]
    cl_polar = data[:, 1]
    cd_polar = data[:, 2]
    
    return alpha_polar_deg, cl_polar, cd_polar

def getforces(alpha_rad, alpha_polar_deg, cl_polar, cd_polar):
    alpha_deg = np.degrees(alpha_rad)
    
    # Interpolate cl and cd
    cl_val = np.interp(alpha_deg, alpha_polar_deg, cl_polar)
    cd_val = np.interp(alpha_deg, alpha_polar_deg, cd_polar)
    
    return cl_val, cd_val


def calcloads(phi, theta, getforces_func):

    alpha = theta - phi
    
    Cl, Cd = getforces_func(alpha)
    
    Cn = Cl * np.cos(phi) - Cd * np.sin(phi)
    Ct = Cl * np.sin(phi) + Cd * np.cos(phi)
    
    return Cl, Cd, Cn, Ct, alpha

def calculate_element_loads(r_local, R, R_root, chord, theta, U0, Omega, B, getforces_func):

    a = 0.0
    a_prime = 0.0
    tolerance = 1e-5
    max_iter = 100
    iteration = 0
    error = 1.0
    
    #solidity
    sigma = (B * chord) / (2 * np.pi * r_local)
    
    while error > tolerance and iteration < max_iter:
        a_old = a
        aprime_old = a_prime
        
        # Flow angle (phi)
        phi = np.arctan((U0 * (1 + a)) / (Omega * r_local * (1 - a_prime)))
        
        #Angle of attack (alpha)
        alpha = theta - phi
        
        # get cl and cd
        Cl, Cd = getforces_func(alpha)
        
        # axial and azimuthal force coefficients
        Cn = Cl * np.cos(phi) -Cd *np.sin(phi)
        Ct = Cl *np.sin(phi)+ Cd * np.cos(phi)
        
        #Prandtl Tip and Root Corrections
        sin_phi = abs(np.sin(phi)) + 1e-8  # Avoid division by zero

        #tip correction
        exp_arg_tip = -(B * (R - r_local)) / (2 * r_local * sin_phi)
        F_tip = (2 / np.pi) * np.arccos(np.exp(exp_arg_tip))

        #root correction
        exp_arg_root = -(B * (r_local - R_root)) / (2 * r_local * sin_phi)
        F_root = (2 / np.pi) * np.arccos(np.exp(exp_arg_root))    

        #combined correction
        F = F_tip * F_root

        #Calculate new induction factors
        # Adding a small number (1e-8) to denominators to prevent division by zero
        a = (sigma * Cn) / (4 * F * (np.sin(phi)**2) - sigma * Cn + 1e-8)
        a_prime = (sigma * Ct) / (4 * F * np.sin(phi) * np.cos(phi) + sigma * Ct + 1e-8)
        
        # 7. Convergence check
        error = abs(a - a_old) + abs(a_prime - aprime_old)
        iteration += 1
        
    # Calculate loads
    V_rel =(U0*(1+a))/np.sin(phi)
    rho = 1.0065 # Isa adjusted for 2000m altitude
    dT = 0.5 * rho * (V_rel**2) * chord * Cn
    dQ = 0.5 * rho * (V_rel**2) * chord * Ct * r_local
    
    return dT, dQ, a, a_prime, phi, alpha
