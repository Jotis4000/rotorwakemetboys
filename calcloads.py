import numpy as np
import matplotlib.pyplot as plt
import momcorrections

def load_airfoil_data(path):

    data = np.loadtxt(path, skiprows=2)
    
    alpha_polar_deg = data[:, 0]
    cl_polar = data[:, 1]
    cd_polar = data[:, 2]
    
    return alpha_polar_deg, cl_polar, cd_polar

# alpha_polar_deg, cl_polar, cd_polar = load_airfoil_data('data/ARAD8pct_polar.txt')
# plt.plot(alpha_polar_deg,cl_polar)
# # plt.plot(alpha_polar_deg,cl_polar/cd_polar)
# plt.show()

def getforces(alpha_rad, alpha_polar_deg, cl_polar, cd_polar):
    alpha_deg = np.degrees(alpha_rad)
    
    # Interpolate cl and cd
    cl_val = np.interp(alpha_deg, alpha_polar_deg, cl_polar)
    cd_val = np.interp(alpha_deg, alpha_polar_deg, cd_polar)
    
    return cl_val, cd_val


def corrections(phi,B,R,r_local,R_root):

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

    return F

def calculate_turbine3tokyodrift(r_local, R, R_root, chord, U0, Omega, B, foilpath):

    a = 0.3
    a_prime = 0.0
    tolerance = 1e-6
    max_iter = 1000
    iteration = 0
    error = 1.0
    phi = 0
    rho = 1.0065 # Isa adjusted for 2000m altitude
    lada = Omega*R/U0
    alpha = np.radians(-9)
    Cl = -0.5365
    Cd =  0.01863
    
    # alpha_polar_deg, cl_polar, cd_polar = load_airfoil_data(foilpath)

    while error > tolerance and iteration < max_iter:
        a_old = a
        aprime_old = a_prime
        
        Vax = U0*(1-a)
        Vtan = Omega*r_local*(1+a_prime)
        Vp = np.sqrt(Vax**2+Vtan**2)

        phi = np.arctan2((U0 * (1 - a)),(Omega * r_local * (1 + a_prime)))

        theta = alpha+phi
        
        lift = 0.5*chord*rho*Vp**2*Cl
        drag = 0.5*chord*rho*Vp**2*Cd

        Faz = lift*np.sin(phi)-drag*np.cos(phi)
        Fax = lift*np.cos(phi)+drag*np.sin(phi)

        # Calculate inductions w corrections

        F = corrections(phi,B,R,r_local,R_root)
        if F==0:
            F=0.27

        CT = (Fax*B)/(0.5*rho*U0**2*2*np.pi*r_local)
        CT1 = 1.816
        CT2 = 2*np.sqrt(CT1)-CT1
        if CT<CT2:
            a_new = 0.5-np.sqrt(1-CT)/2
        else:
            a_new = 1+(CT-CT1)/(4*np.sqrt(CT1)-4)

        a_prime_new = (Faz*B)/(2*rho*2*np.pi*r_local*U0**2*(1-a_new)*lada*r_local/R)

        a_new = a_new/F
        a_prime_new = a_prime_new/F

        # Apply relaxation to aid convergence
        a = 0.25 * a_new + 0.75 * a_old
        a_prime = 0.25 * a_prime_new + 0.75 * aprime_old
        if a>=0.95: a=0.95
        if a_prime>0.95: a_prime=0

        error = abs(a - a_old) + abs(a_prime - aprime_old)
        # print(error)
        iteration += 1

    # Calculate loads
    dT = Fax*B
    dQ = Faz*B*r_local

    # print(dT)
    # print(dQ)

    return dT, dQ, a, a_prime, phi, theta, F

def calculate_turbine2(r_local, R, R_root, chord, theta, U0, Omega, sigma, B, foilpath):

    a = 0.3
    a_prime = 0.0
    tolerance = 1e-6
    max_iter = 1000
    iteration = 0
    error = 1.0
    phi = 0
    rho = 1.0065 # Isa adjusted for 2000m altitude
    lada = Omega*R/U0
    
    alpha_polar_deg, cl_polar, cd_polar = load_airfoil_data(foilpath)

    while error > tolerance and iteration < max_iter:
        a_old = a
        aprime_old = a_prime
        
        Vax = U0*(1-a)
        Vtan = Omega*r_local*(1+a_prime)
        Vp = np.sqrt(Vax**2+Vtan**2)

        phi = np.arctan2((U0 * (1 - a)),(Omega * r_local * (1 + a_prime)))
        # alpha = phi - theta
        alpha = theta-phi ###CHANGED COMPARE TO ELEMENT LOAD 3

        Cl, Cd = getforces(alpha, alpha_polar_deg, cl_polar, cd_polar)

        lift = 0.5*chord*rho*Vp**2*Cl
        drag = 0.5*chord*rho*Vp**2*Cd

        Faz = lift*np.sin(phi)-drag*np.cos(phi)
        Fax = lift*np.cos(phi)+drag*np.sin(phi)

        # Calculate inductions w corrections
        
        F = corrections(phi,B,R,r_local,R_root)
        if F==0:
            F=0.27

        CT = (Fax*B)/(0.5*rho*U0**2*2*np.pi*r_local)
        CT1 = 1.816
        CT2 = 2*np.sqrt(CT1)-CT1
        if CT<CT2:
            a_new = 0.5-np.sqrt(1-CT)/2
        else:
            a_new = 1+(CT-CT1)/(4*np.sqrt(CT1)-4)

        a_prime_new = (Faz*B)/(2*rho*2*np.pi*r_local*U0**2*(1-a_new)*lada*r_local/R)

        a_new = a_new/F
        a_prime_new = a_prime_new/F

        # Apply relaxation to aid convergence
        a = 0.25 * a_new + 0.75 * a_old
        a_prime = 0.25 * a_prime_new + 0.75 * aprime_old
        if a>=0.95: a=0.95
        if a_prime>0.95: a_prime=0

        error = abs(a - a_old) + abs(a_prime - aprime_old)
        # print(error)
        iteration += 1

    # Calculate loads
    dT = Fax*B
    dQ = Faz*B*r_local

    # print(dT)
    # print(dQ)

    return dT, dQ, a, a_prime, phi, alpha, F

def calculate_turbine(r_local, R, R_root, chord, theta, U0, Omega, sigma, B, foilpath):

    a = 0.3
    a_prime = 0.0
    tolerance = 1e-6
    max_iter = 1000
    iteration = 0
    error = 1.0
    phi = 0
    
    alpha_polar_deg, cl_polar, cd_polar = load_airfoil_data(foilpath)

    while error > tolerance and iteration < max_iter:
        a_old = a
        aprime_old = a_prime
        
        phi = np.arctan((U0 * (1 - a)) / (Omega * r_local * (1 + a_prime)))
        alpha = phi - theta

        # get cl and cd
        Cl, Cd = getforces(alpha, alpha_polar_deg, cl_polar, cd_polar)
        
        F = corrections(phi,B,R,r_local,R_root)
        if F==0:
            F=0.27

        # Normal force (thrust direction) - Drag adds to thrust
        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)

        # Tangential force (torque direction) - Lift drives rotation, Drag opposes it!
        Ct = Cl * np.sin(phi) - Cd * np.cos(phi)

        # 1. Update a_prime (Tangential induction)
        # Prevent division by zero
        if Ct != 0:
            a_prime_new = 1.0 / (((4 * F * np.sin(phi) * np.cos(phi)) / (sigma * Ct)) - 1.0)
        else:
            a_prime_new = 0.0

        # 2. Update a (Axial induction) with Glauert Correction
        # Calculate local thrust coefficient to check if correction is needed
        CT_local = (sigma * Cn * (1 - a)**2) / (np.sin(phi)**2)

        if CT_local > 0.96: 
            # Buhl's empirical Glauert correction for high induction
            a_new = (1 / F) * (0.143 + np.sqrt(0.0203 - 0.6427 * (0.889 - CT_local)))
            print('weird stuff')
        else:
            # Standard momentum theory update
            if Cn != 0:
                a_new = 1.0 / (((4 * F * np.sin(phi)**2) / (sigma * Cn)) + 1.0)
            else:
                a_new = 0.0

        # Apply relaxation to aid convergence
        a = 0.25 * a_new + 0.75 * a_old
        a_prime = 0.25 * a_prime_new + 0.75 * aprime_old

        error = abs(a - a_old) + abs(a_prime - aprime_old)
        # print(error)
        iteration += 1

    # Calculate loads
    V_rel =(U0*(1-a))/np.sin(phi)
    rho = 1.0065 # Isa adjusted for 2000m altitude
    dT = 0.5 * rho * (V_rel**2) * chord * Cn * B
    dQ = 0.5 * rho * (V_rel**2) * chord * Ct * r_local * B

    return dT, dQ, a, a_prime, phi, alpha, F

def calculate_element_loads3(r_local, R, R_root, chord, theta, U0, Omega, sigma, B, J, foilpath):

    a = 0.3
    a_prime = 0.0
    tolerance = 1e-6
    max_iter = 1000
    iteration = 0
    error = 1.0
    phi = 0
    
    alpha_polar_deg, cl_polar, cd_polar = load_airfoil_data(foilpath)

    while error > tolerance and iteration < max_iter:
        a_old = a
        aprime_old = a_prime
        
        # Flow angle (phi)
        phi = np.arctan2((U0 * (1 + a)),(Omega * r_local * (1 - a_prime)))
        
        #Angle of attack (alpha)
        alpha = theta - phi
        
        # get cl and cd
        Cl, Cd = getforces(alpha, alpha_polar_deg, cl_polar, cd_polar)
        
        F = corrections(phi,B,R,r_local,R_root)
        if F==0:
            F=0.27
            # F = 

        f1 = sigma*(Cl*np.cos(phi))/(4*F*np.sin(phi)**2)
        f2 = sigma*(Cl)/(4*F*np.cos(phi))

        Cn = Cl * np.cos(phi) -Cd *np.sin(phi)
        Ct = Cl *np.sin(phi)+ Cd * np.cos(phi)

        a = f1/(1-f1)
        a_prime = f2/(1+f2)

        a = 0.25*a+0.75*a_old
        a_prime = 0.25*a_prime+0.75*aprime_old

        if a>=0.95: a=0.95
        if a_prime>0.95: a_prime=0

        # 7. Convergence check
        error = abs(a - a_old) + abs(a_prime - aprime_old)
        # print(error)
        iteration += 1
    
    # print("Finished in "+str(iteration)+" iterations.")

    # print(a)  
    # print(a_prime)
    # print(alpha)

    # Calculate loads
    V_rel =(U0*(1+a))/np.sin(phi)
    rho = 1.0065 # Isa adjusted for 2000m altitude
    dT = 0.5 * rho * (V_rel**2) * chord * Cn * B
    dQ = 0.5 * rho * (V_rel**2) * chord * Ct * r_local * B
    # dT=dQ=1
    
    return dT, dQ, a, a_prime, phi, alpha, F

def calculate_element_loads2(r_local, R, R_root, chord, theta, U0, Omega, sigma, B, J, foilpath):

    a = 0.3
    a_prime = 0.0
    tolerance = 1e-6
    max_iter = 1000
    iteration = 0
    error = 1.0
    rho = 1.0065
    lada = Omega*R/U0

    # print(lada)

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
        Vax = U0*(1+a_old)
        Vtan = Omega*r_local*(1-aprime_old)
        Vp = np.sqrt(Vax**2+Vtan**2)

        lift = 0.5*chord*rho*Vp**2*Cl
        drag = 0.5*chord*rho*Vp**2*Cd

        Faz = lift*np.sin(phi)-drag*np.cos(phi)
        Fax = lift*np.cos(phi)+drag*np.sin(phi)

        F = corrections(phi,B,R,r_local,R_root)
        if F==0:
            F=0.27

        CT = (Fax*B)/(0.5*rho*U0**2*2*np.pi*r_local)
        a = 0.5*(np.sqrt(1+CT/F)-1)
        aprime = (Faz*B)/(2*rho*(2*np.pi*r_local)*U0**2*(1-a)*lada*r_local/R)

        # #Prandtl Tip and Root Corrections
        # sin_phi = abs(np.sin(phi)) + 1e-8  # Avoid division by zero

        # #tip correction
        # exp_arg_tip = -(B * (R - r_local)) / (2 * r_local * sin_phi)
        # F_tip = (2 / np.pi) * np.arccos(np.exp(exp_arg_tip))

        # #root correction
        # exp_arg_root = -(B * (r_local - R_root)) / (2 * r_local * sin_phi)
        # # print(np.exp(exp_arg_root))
        # F_root = (2 / np.pi) * np.arccos(np.exp(exp_arg_root))    

        # #combined correction
        # F = F_tip * F_root

        # a = a/F
        # aprime = aprime/F

        # a = 0.5*(np.sqrt(1+CT/F)-1)
        # aprime = (Faz*B)/(2*rho*(2*np.pi*r_local)*U0**2*(1-a)*lada*r_local/R)

        a = 0.25*a+0.75*a_old
        a_prime = 0.25*a_prime+0.75*aprime_old
        if a>=0.95: a=0.95

        Cn = Cl * np.cos(phi) - Cd *np.sin(phi)
        Ct = Cl *np.sin(phi) + Cd * np.cos(phi)

        # convergence check
        error = abs(a - a_old) + abs(a_prime - aprime_old)
        # print(error)
        iteration += 1

    # print("Finished in "+str(iteration)+" iterations.")

    # print(a)  
    # print(a_prime)
    # print(alpha)

    # Calculate loads
    V_rel =(U0*(1+a))/np.sin(phi)
    rho = 1.0065 # Isa adjusted for 2000m altitude
    # dT = 0.5 * rho * (V_rel**2) * chord * Cn * B
    # dQ = 0.5 * rho * (V_rel**2) * chord * Ct * r_local * B
    dT = Fax * B
    dQ = Faz * B * r_local

    return dT, dQ, a, aprime, phi, alpha, F

def calculate_element_loads(r_local, R, R_root, chord, theta, U0, Omega, sigma, B, J, foilpath):

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
        a = (sigma * Cn) / (4 * F * (np.sin(phi)**2) - sigma * Cn)
        a_prime = (sigma * Ct) / (4 * F * np.sin(phi) * np.cos(phi) + sigma * Ct)
        
        a = 0.25*a+0.75*a_old
        a_prime = 0.25*a_prime+0.75*aprime_old

        # if a>1: a=1
        # if a_prime>1: a_prime=1

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
    dT = 0.5 * rho * (V_rel**2) * chord * Cn * B
    dQ = 0.5 * rho * (V_rel**2) * chord * Ct * r_local * B
    
    return dT, dQ, a, a_prime, phi, alpha, F
