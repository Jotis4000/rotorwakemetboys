import numpy as np
import matplotlib.pyplot as plt

def PrandtlTipLossCorrection(x, x_root, lada, B, a, ap):
# 
#   x: r_R 
#   x_root: location of the root vortex 
#   lada: tip speed ratio
#   B: number of blades
#   a: axial induction
#   ap: azimuthal induction factor

    exptip = -(B/2) * ((1-x)/x) * np.sqrt(1 + (lada**2 * x**2)/((1-a)**2))
    ftip = 2 / np.pi * np.arccos(np.exp(exptip))

    if np.isnan(ftip):
        ftip = 0.0

    exproot = -(B/2) * ((x-x_root)/x) * np.sqrt(1 + (lada**2 * x**2)/((1-a)**2))
    froot = 2 / np.pi * np.arccos(np.exp(exproot))

    if np.isnan(froot):
        froot = 0.0

    # if froot < 0.000001:
    #     froot = 0.000001

    ftot = ftip * froot
    a_corr = a / ftot
    ap_corr = ap / ftot

    return ftot, ftip, froot


    # Rotor parameters
B = 3          
lada = 8.0     
x_root = 0   
a_guess = 0.3  # Assumed uniform axial induction for the sake of the plot
ap_guess = 0.01

# Create radial array, slightly offsetting from exact 0.2 and 1.0 to avoid 0.0 values inside the arccos
radii = np.linspace(x_root + 0.001, 0.999, 100)

F_tot_vals = []
F_tip_vals = []
F_root_vals = []

for x in radii:
    ftot, ftip, froot = PrandtlTipLossCorrection(x, x_root, lada, B, a_guess, ap_guess)
    F_tot_vals.append(ftot)
    F_tip_vals.append(ftip)
    print(ftip)
    F_root_vals.append(froot)

# Generate Plot

plt.figure(figsize=(8, 5))
plt.plot(radii, F_tip_vals, linestyle='--', color='blue')
plt.plot(radii, F_root_vals, linestyle='--', color='red')
# plt.plot(radii, F_tot_vals, label='Total Loss Factor ($F_{tot}$)', linewidth=2, color='black')

plt.title('Prandtl Tip and Root Loss Corrections')
plt.xlabel('Non-dimensional radius ($r/R$)')
plt.ylabel('Correction Factor ($F$)')
plt.xlim([0, 1])
plt.ylim([0, 1.1])
plt.grid(True)
plt.legend()
plt.show()
