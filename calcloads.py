import numpy as np
import matplotlib.pyplot as plt
import momcorrections

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


# def calcloads(phi, theta, getforces_func):

#     alpha = theta - phi
    
#     Cl, Cd = getforces_func(alpha)
    
#     Cn = Cl * np.cos(phi) - Cd * np.sin(phi)
#     Ct = Cl * np.sin(phi) + Cd * np.cos(phi)
    
#     return Cl, Cd, Cn, Ct, alpha

# def calcas(Cl,Cd):

    

#     return a, aprime

def calculate_element_loads2(r_local, R, R_root, chord, theta, U0, Omega, sigma, B, phi, foilpath):

    a = 0.3
    a_prime = 0.0
    tolerance = 1e-6
    max_iter = 1000
    iteration = 0
    error = 1.0
    rho = 1.0065
    lada = Omega*R/U0

    print(lada)

    alpha_polar_deg, cl_polar, cd_polar = load_airfoil_data(foilpath)

    while error > tolerance and iteration < max_iter:
        a_old = a
        aprime_old = a_prime

        # Flow angle (phi)
        phi = np.arctan((U0 * (1 + a)) / (Omega * r_local * (1 - a_prime)))

        #Angle of attack (alpha)
        alpha = theta - phi
        
        # get cl and cd
        Cl, Cd = getforces(alpha, alpha_polar_deg, cl_polar, cd_polar)
        Vax = U0*(1-a_old)
        Vtan = Omega*r_local*(1+aprime_old)
        Vp = np.sqrt(Vax**2+Vtan**2)

        lift = 0.5*chord*rho*Vp**2*Cl
        drag = 0.5*chord*rho*Vp**2*Cd

        Faz = lift*np.sin(phi)-drag*np.cos(phi)
        Fax = lift*np.cos(phi)+drag*np.sin(phi)

        CT = (Fax*B)/(0.5*rho*U0**2*2*np.pi*r_local)
        a = 0.5*(1-np.sqrt(1-CT))
        aprime = (Faz*B)/(2*rho*(2*np.pi*r_local)*U0**2*(1-a)*lada*r_local/R)

        #Prandtl Tip and Root Corrections
        sin_phi = abs(np.sin(phi)) + 1e-8  # Avoid division by zero

        #tip correction
        exp_arg_tip = -(B * (R - r_local)) / (2 * r_local * sin_phi)
        F_tip = (2 / np.pi) * np.arccos(np.exp(exp_arg_tip))

        #root correction
        exp_arg_root = -(B * (r_local - R_root)) / (2 * r_local * sin_phi)
        # print(np.exp(exp_arg_root))
        F_root = (2 / np.pi) * np.arccos(np.exp(exp_arg_root))    

        #combined correction
        F = F_tip * F_root

        a = a/F
        aprime = aprime/F

        a = 0.25*a+0.75*a_old
        a_prime = 0.25*a_prime+0.75*aprime_old

        # convergence check
        error = abs(a - a_old) + abs(a_prime - aprime_old)
        # print(error)
        iteration += 1

    print("Finished in "+str(iteration)+" iterations.")

    # print(a)  
    # print(a_prime)
    # print(alpha)

    # Calculate loads
    V_rel =(U0*(1+a))/np.sin(phi)
    rho = 1.0065 # Isa adjusted for 2000m altitude
    dT = 1
    dQ = 1

    return dT, dQ, a, aprime, phi, alpha, F

def calculate_element_loads(r_local, R, R_root, chord, theta, U0, Omega, sigma, B, phi, foilpath):

    a = 0.3
    a_prime = 0.0
    tolerance = 1e-6
    max_iter = 1000
    iteration = 0
    error = 1.0
    
    
    #solidity
    # sigma = (B * chord) / (2 * np.pi * r_local)
    
    alpha_polar_deg, cl_polar, cd_polar = load_airfoil_data(foilpath)

    while error > tolerance and iteration < max_iter:
        a_old = a
        aprime_old = a_prime
        
        # Flow angle (phi)
        phi = np.arctan((U0 * (1 + a)) / (Omega * r_local * (1 - a_prime)))
        
        #Angle of attack (alpha)
        alpha = theta - phi
        
        # get cl and cd
        Cl, Cd = getforces(alpha, alpha_polar_deg, cl_polar, cd_polar)
        
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
        # print(np.exp(exp_arg_root))
        F_root = (2 / np.pi) * np.arccos(np.exp(exp_arg_root))    

        #combined correction
        F = F_tip * F_root

        # F, ftip, froot = momcorrections.PrandtlTipLossCorrection(r_local/R, R_root/R, Omega*R/U0, B, a, a_prime)

        #Calculate new induction factors
        # Adding a small number (1e-8) to denominators to prevent division by zero
        a = (sigma * Cn) / (4 * F * (np.sin(phi)**2) - sigma * Cn + 1e-8)
        a_prime = (sigma * Ct) / (4 * F * np.sin(phi) * np.cos(phi) + sigma * Ct + 1e-8)
        
        a = 0.25*a+0.75*a_old
        a_prime = 0.25*a_prime+0.75*aprime_old

        # 7. Convergence check
        error = abs(a - a_old) + abs(a_prime - aprime_old)
        # print(error)
        iteration += 1
    
    print("Finished in "+str(iteration)+" iterations.")

    # print(a)  
    # print(a_prime)
    # print(alpha)

    # Calculate loads
    V_rel =(U0*(1+a))/np.sin(phi)
    rho = 1.0065 # Isa adjusted for 2000m altitude
    dT = 0.5 * rho * (V_rel**2) * chord * Cn
    dQ = 0.5 * rho * (V_rel**2) * chord * Ct * r_local
    
    return dT, dQ, a, a_prime, phi, alpha, F
