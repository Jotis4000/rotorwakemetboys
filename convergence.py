import numpy as np
import matplotlib.pyplot as plt
import defgeom
import calcloads
from scipy.interpolate import make_interp_spline

start = 0.25
R = 0.7
B = 6
U0 = 60
airfoil = "data/ARAD8pct_polar.txt"
pitch = np.radians(46)

J_values = np.array([1.6, 2.0, 2.4])
RPM_values = U0 / (J_values * 2 * R) * 60

#number of annuli
N_values = np.array([10, 20, 30, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800])

thrust_linear = np.zeros((len(J_values), len(N_values)))
thrust_cosine = np.zeros((len(J_values), len(N_values)))

for j_idx, J in enumerate(J_values):
    omega = RPM_values[j_idx] / 60 * 2 * np.pi
    print(f"J = {J}...")
    
    for n_idx, N in enumerate(N_values):
        #linear spacing
        x_lin = np.linspace(start, 1, N)
        dr_lin = np.diff(x_lin) * R
        
        tdist_lin = -np.radians(50) * x_lin + np.radians(35)
        cdist_lin = 0.18 - 0.06 * x_lin
        
        c_lin, theta_lin, sigma_lin = defgeom.defGeom(R, B, x_lin, start, tdist_lin, pitch, cdist_lin)
        
        T_lin_total = 0
        for i in range(N-1):
            res = calcloads.calculate_element_loads3(x_lin[i]*R, R, start*R, c_lin[i], theta_lin[i], 
                                                     U0, omega, sigma_lin[i], B, J, airfoil)
            T_lin_total += res[0] * dr_lin[i] # dT * dr
            
        thrust_linear[j_idx, n_idx] = T_lin_total
        
        #cosine spacing
        beta = np.linspace(0, np.pi, N)
        x_cos = start + (1 - start) * 0.5 * (1 - np.cos(beta))
        dr_cos = np.diff(x_cos) * R
        
        tdist_cos = -np.radians(50) * x_cos + np.radians(35)
        cdist_cos = 0.18 - 0.06 * x_cos
        
        c_cos, theta_cos, sigma_cos = defgeom.defGeom(R, B, x_cos, start, tdist_cos, pitch, cdist_cos)
        
        T_cos_total = 0
        for i in range(N-1):
            res = calcloads.calculate_element_loads3(x_cos[i]*R, R, start*R, c_cos[i], theta_cos[i], 
                                                     U0, omega, sigma_cos[i], B, J, airfoil)
            T_cos_total += res[0] * dr_cos[i] # dT * dr
            
        thrust_cosine[j_idx, n_idx] = T_cos_total


plt.rcParams.update({
    "text.usetex": False, 
    "font.family": "serif",
    "font.size": 12,
    "axes.grid": True,
    "grid.linestyle": '--',
    "grid.alpha": 0.7
})

fig, axs = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
fig.suptitle("Convergence History of Total Thrust for Different Advance Ratios", fontsize=16, y=1.05)

for j_idx in range(len(J_values)):
    axs[j_idx].plot(N_values, thrust_linear[j_idx, :], marker='o', color='black', label="Constant Spacing")
    axs[j_idx].plot(N_values, thrust_cosine[j_idx, :], marker='s', color='red', linestyle='--', label="Cosine Spacing")
    
    axs[j_idx].set_title(f"$J = {J_values[j_idx]}$")
    axs[j_idx].set_xlabel("Number of Annuli ($N$)")
    if j_idx == 0:
        axs[j_idx].set_ylabel("Total Thrust [N]")
    axs[j_idx].legend()

plt.tight_layout()
plt.show()

