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
J = np.array([1.6, 2.0, 2.4]) 
RPM = (U0 / (J * 2 * R)) * 60 

print("Advance Ratios J:", J)
print("RPMs:", RPM)

c, theta, sigma = defgeom.defGeom(R, B, x, start, tdist, pitch, cdist)

results_corr = np.zeros([len(J), len(x)-1, 8])
results_uncorr = np.zeros([len(J), len(x)-1, 8])

for j in range(len(J)):
    omega = (RPM[j] / 60) * 2 * np.pi
    for i in range(len(x)-1):
        # 1. Normal model (WITH Tip Corrections)
        results_corr[j, i, :] = calcloads.calculate_element_loads3(
            x[i]*R, R, start*R, c[i], theta[i], U0, omega, sigma[i], B, J[j], airfoil)
        
        # 2. Uncorrected model (WITHOUT Tip Corrections)
        results_uncorr[j, i, :] = calcloads.calculate_element_loads_notipcorr(
            x[i]*R, R, start*R, c[i], theta[i], U0, omega, sigma[i], B, J[j], airfoil)

T_corr = np.zeros(len(J))
Q_corr = np.zeros(len(J))
T_uncorr = np.zeros(len(J))
Q_uncorr = np.zeros(len(J))

for i in range(len(J)):
    # Integrate Corrected
    T_corr[i] = np.sum(results_corr[i, :, 0] * np.diff(x) * R)
    Q_corr[i] = np.sum(results_corr[i, :, 1] * np.diff(x) * R)
    
    # Integrate Uncorrected
    T_uncorr[i] = np.sum(results_uncorr[i, :, 0] * np.diff(x) * R)
    Q_uncorr[i] = np.sum(results_uncorr[i, :, 1] * np.diff(x) * R)

print("\n--- Corrected Global Results ---")
print("Thrust:", T_corr)
print("Torque:", Q_corr)
print("\n--- Uncorrected Global Results ---")
print("Thrust:", T_uncorr)
print("Torque:", Q_uncorr)


plt.rcParams.update({
    "text.usetex": False, 
    "font.family": "serif",
    "font.size": 12,
    "axes.grid": True,
    "grid.linestyle": '--',
    "grid.alpha": 0.7
})


x_plot = x[1:(el-1)]
colors = ["black", "red", "green"]

# alpha
fig1 = plt.figure(figsize=(10,6))
plt.title(r"Influence of Tip Correction on $\alpha$")
for j in range(len(J)):
    plt.plot(x_plot, np.degrees(results_corr[j, :, 5][1:]), color=colors[j], linestyle='-', label=f"$J={J[j]:.2f}$ (Corr)")
    plt.plot(x_plot, np.degrees(results_uncorr[j, :, 5][1:]), color=colors[j], linestyle='--', label=f"$J={J[j]:.2f}$ (Uncorr)")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Angle of Attack $\alpha$ [$\degree$]")
plt.legend(ncol=3)
plt.tight_layout()

# phi
fig2 = plt.figure(figsize=(10,6))
plt.title(r"Influence of Tip Correction on $\phi$")
for j in range(len(J)):
    plt.plot(x_plot, np.degrees(results_corr[j, :, 4][1:]), color=colors[j], linestyle='-', label=f"$J={J[j]:.2f}$ (Corr)")
    plt.plot(x_plot, np.degrees(results_uncorr[j, :, 4][1:]), color=colors[j], linestyle='--', label=f"$J={J[j]:.2f}$ (Uncorr)")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Inflow Angle $\phi$ [$\degree$]")
plt.legend(ncol=3)
plt.tight_layout()

# inductions a and a'
fig3, axs = plt.subplots(1, 3, figsize=(18,6))
plt.suptitle(r"Influence of Tip Correction on Inductions $a$ and $a'$")
for j in range(len(J)):
    # Corrected
    axs[j].plot(x_plot, results_corr[j, :, 2][1:], label=r"Axial $a$ (Corr)", color="black", linestyle='-')
    axs[j].plot(x_plot, results_corr[j, :, 3][1:], label=r"Azim. $a'$ (Corr)", color="red", linestyle='-')
    # Uncorrected
    axs[j].plot(x_plot, results_uncorr[j, :, 2][1:], label=r"Axial $a$ (Uncorr)", color="black", linestyle='--')
    axs[j].plot(x_plot, results_uncorr[j, :, 3][1:], label=r"Azim. $a'$ (Uncorr)", color="red", linestyle='--')
    
    axs[j].set_title(f"J={J[j]:.2f}")
    axs[j].set_xlabel("$r/R$ [-]")
    axs[j].set_ylabel("Induction [-]")
    axs[j].legend()
plt.tight_layout()

# thrust
fig4 = plt.figure(figsize=(10,6))
plt.title(r"Influence of Tip Correction on $dT$")
for j in range(len(J)):
    plt.plot(x_plot, results_corr[j, :, 0][1:], color=colors[j], linestyle='-', label=f"$J={J[j]:.2f}$ (Corr)")
    plt.plot(x_plot, results_uncorr[j, :, 0][1:], color=colors[j], linestyle='--', label=f"$J={J[j]:.2f}$ (Uncorr)")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Thrust Distribution [$N/m$]")
plt.legend(ncol=3)
plt.tight_layout()

# torque
fig5 = plt.figure(figsize=(10,6))
plt.title(r"Influence of Tip Correction on $dQ$")
for j in range(len(J)):
    plt.plot(x_plot, results_corr[j, :, 1][1:], color=colors[j], linestyle='-', label=f"$J={J[j]:.2f}$ (Corr)")
    plt.plot(x_plot, results_uncorr[j, :, 1][1:], color=colors[j], linestyle='--', label=f"$J={J[j]:.2f}$ (Uncorr)")
plt.xlabel("$r/R$ [-]")
plt.ylabel(r"Torque Distribution [$Nm/m$]")
plt.legend(ncol=3)
plt.tight_layout()

plt.show()