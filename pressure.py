import numpy as np
import matplotlib.pyplot as plt
import defgeom
import calcloads

start = 0.25
el = 200
x = np.linspace(start, 1, el)

R = 0.7
B = 6
tdist = -np.radians(50)*x + np.radians(35) # for x>0.25
pitch = np.radians(46) 
cdist = 0.18 - 0.06*x # for x>0.25
airfoil = "data/ARAD8pct_polar.txt"

U0 = 60
rho = 1.00649
pinf = 79495

J = np.array([1.6, 2.0, 2.4]) 
RPM = U0/(J*2*R)*U0 

c, theta, sigma = defgeom.defGeom(R, B, x, start, tdist, pitch, cdist)


results = np.zeros([len(J), len(x)-1, 8])

for j in range(len(J)):
    for i in range(len(x)-1):
        results[j,i,:] = calcloads.calculate_element_loads3(
            x[i]*R, R, start*R, c[i], theta[i], U0, RPM[j]/60*2*np.pi, sigma[i], B, J[j], airfoil
        )

###Pressure Calculations and plotting
plt.rcParams.update({
    "text.usetex": False, 
    "font.family": "serif",
    "font.size": 12,
    "axes.grid": True,
    "grid.linestyle": '--',
    "grid.alpha": 0.7
})

r_local = x[:(el-1)] * R
D = 2 * R

# upwind
p1 = pinf + 0.5 * rho * U0**2
# 2. rotor upwind (same as freestream total pressure)
p2 = p1

figpres, axs = plt.subplots(1, 3, figsize=(18, 5))
plt.suptitle("Stagnation Pressure Distribution across Advance Ratios", fontsize=14)

for idx, j_val in enumerate(J):
    
    # Calculate Omega dynamically for the current J
    Omega = 2 * np.pi * (U0 / (j_val * D))
    
    # Extract induction factors and tip correction
    a = results[idx, :, 2]
    a_prime = results[idx, :, 3]  
    F = results[idx, :, 6]

    # Calculate Pressure Jumps
    dp_static = 2 * rho * U0**2 * a * (1 + a) * F
    dp_swirl = 0.5 * rho * (2 * a_prime * Omega * r_local)**2
    dp_total = dp_static + dp_swirl

    p3 = p2 + dp_total
    p4 = p3 

    # Streamtube radii via Continuity
    r_upstream = r_local * np.sqrt(1 + a)               # Upwind
    r_contracted = r_local * np.sqrt((1 + a) / (1 + 2 * a)) # Downwind

    # --- PLOTTING ---
    axs[idx].plot(r_upstream, np.full_like(r_upstream, p1)/1000, label="1. Infinity Upwind", color="black", linestyle="-", linewidth=4)
    axs[idx].plot(r_local, np.full_like(r_local, p2)/1000, label="2. Rotor Upwind", color="orange", linestyle="--", linewidth=2)
    axs[idx].plot(r_local, p3/1000, label="3. Rotor Downwind", color="red", linestyle="-", linewidth=4)
    axs[idx].plot(r_contracted, p4/1000, label="4. Infinity Downwind", color="blue", linestyle="--", linewidth=2)


    r_edge_down = r_contracted[-1]
    if r_edge_down < R:
        axs[idx].plot([r_edge_down, R], [p1/1000, p1/1000], color="blue", linestyle="--", linewidth=2)
    else: 
        axs[idx].plot([R, r_edge_down], [p1/1000, p1/1000], color="blue", linestyle="--", linewidth=2)

    axs[idx].set_title(f"$J = {j_val:.2f}$")
    axs[idx].set_xlabel("Radial Position $r$ [m]")
    axs[idx].set_xlim([0.1, max(r_upstream[-1], R) * 1.05])
    
    if idx == 0:
        axs[idx].set_ylabel("Total (Stagnation) Pressure [kPa]")
        axs[idx].legend(loc="lower left", fontsize=9)

    # Sanity Checks
    thrust_integrated = np.trapz(dp_static * 2 * np.pi * r_local, x=r_local)
    print(f"J = {j_val:.2f} | Integrated Thrust from Pressure: {thrust_integrated:.2f} N")

plt.tight_layout()
plt.show()