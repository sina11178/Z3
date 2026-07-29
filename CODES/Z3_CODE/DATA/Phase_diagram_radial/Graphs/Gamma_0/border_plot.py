from matplotlib.colors import TwoSlopeNorm
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from scipy.ndimage import gaussian_filter, map_coordinates
 
# ── 1. Grid — must match the heatmap script exactly ──────────────────────────
L6, L7, L8 = 6, 7, 8
Nr = 40
N_phi = 40
r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2 * np.pi / 3, N_phi)
R, Phi = np.meshgrid(r_vals, phi_vals, indexing='ij')
 
# ── 2. Load data — both r-major, shape (Nr, N_phi) ────────────────────────────
tau_6 = 1 / np.load('lifetimes_L6_gamma_NPY.npy')
tau_7 = 1 / np.load('lifetimes_L7_gamma_NPY.npy')   # shape (Nr, N_phi)
tau_8 = 1 / np.load('lifetimes_L8_gamma_NPY.npy')   # shape (Nr, N_phi)

tau_6[0][0] = np.exp(6)
tau_7[0][0] = np.exp(7)
tau_8[0][0] = np.exp(8)

'''
tau_6[0][0] = np.exp(L6)
tau_7[0][0] = np.exp(L7)
tau_8[0][0] = np.exp(L8)
'''

assert tau_6.shape == (Nr, N_phi), f"tau_7 shape {tau_6.shape} != ({Nr},{N_phi})"
assert tau_7.shape == (Nr, N_phi), f"tau_7 shape {tau_7.shape} != ({Nr},{N_phi})"
assert tau_8.shape == (Nr, N_phi), f"tau_8 shape {tau_8.shape} != ({Nr},{N_phi})"
 
# ── 3. Scale: f_L = ln(tau) / sqrt(L) ─────────────────────────────────────────
floor = 1e-14
pow = 0.5
f6 = np.log(np.clip(tau_6, floor, None)) / np.power(L6, pow)
f7 = np.log(np.clip(tau_7, floor, None)) / np.power(L7, pow)
f8 = np.log(np.clip(tau_8, floor, None)) / np.power(L8, pow)

 
# ── 4. Signed difference field: > 0 where L8 scaling exceeds L7, < 0 where it
#      falls below. The phase border is the zero-crossing of this field.
#Delta_1 = f8 - f7   # shape (Nr, N_phi)
#Delta_2 = f7 - f6
#Delta = Delta_2 - Delta_1
Delta = f8 - f7
Delta = gaussian_filter(Delta, sigma=0.5)
plt.imshow(f8[::-1, :], cmap='RdBu')
plt.colorbar()
plt.show()
minimum_delta = np.min(np.abs(Delta), axis=0)

plt.plot(minimum_delta)
plt.show()

# ---------------------------------------------
# Plot
# ---------------------------------------------
fig = plt.figure(figsize=(8,7))
ax = plt.subplot(projection='polar')


vmax = np.nanpercentile(np.abs(Delta), 98)
norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)

r_cut = 1
R_plot = R[r_cut:, :]
Phi_plot = Phi[r_cut:, :]
Delta_plot = Delta[r_cut:, :]


# Filled contour of the difference
c = ax.contourf(
    Phi_plot,
    R_plot,
    Delta_plot,
    levels=50,
    cmap='RdBu_r',
    norm=norm, extend='both'
)

# Zero contour = estimated phase boundary

boundary = ax.contour(
    Phi_plot,
    R_plot,
    Delta_plot,
    levels=[0],
    colors='black',
    linewidths=2
)

phi_line_1 = np.pi/4  # radians
#phi_line_2 = 0.9*np.pi/4  # radians
phi_line_3 = np.pi/6  # radians
phi_line_4 = np.pi/3  # radians

print(phi_line_4)

#print(phi_line)
ax.plot([phi_line_1, phi_line_1], [0, 1], linestyle='--', color='green', linewidth=2, label=r'$\phi=pi/4$ rad')
#ax.plot([phi_line_2, phi_line_2], [0, 1], linestyle='--', color='red', linewidth=2, label=r'$\phi=0.9*pi/4$ rad')
ax.plot([phi_line_3, phi_line_3], [0, 1], linestyle='--', color='blue', linewidth=2, label=r'$\phi=pi/6$ rad')
ax.plot([phi_line_4, phi_line_4], [0, 1], linestyle='--', color='black', linewidth=2, label=r'$\phi=pi/3$ rad')
ax.legend(loc='upper right')
ax.set_rlim(0, 0.8)
plt.colorbar(c, pad=0.1, label=r'$\Delta = \ln(\tau)/\sqrt{L}|_{L=8}-\ln(\tau)/\sqrt{L}|_{L=7}$')

ax.set_title("Estimated Phase Boundary")

plt.show()