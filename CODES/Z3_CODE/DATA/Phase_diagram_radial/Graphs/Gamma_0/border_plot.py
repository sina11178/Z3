


import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from scipy.ndimage import gaussian_filter, map_coordinates
 
# ── 1. Grid — must match the heatmap script exactly ──────────────────────────
L7, L8 = 7, 8
Nr = 40
N_phi = 40
r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2 * np.pi / 3, N_phi)
R, Phi = np.meshgrid(r_vals, phi_vals, indexing='ij')
 
# ── 2. Load data — both r-major, shape (Nr, N_phi) ────────────────────────────
tau_7 = 1 / np.load('lifetimes_L7_gamma_NPY.npy')   # shape (Nr, N_phi)
tau_8 = 1 / np.load('lifetimes_L8_gamma_NPY.npy')   # shape (Nr, N_phi)
 
assert tau_7.shape == (Nr, N_phi), f"tau_7 shape {tau_7.shape} != ({Nr},{N_phi})"
assert tau_8.shape == (Nr, N_phi), f"tau_8 shape {tau_8.shape} != ({Nr},{N_phi})"
 
# ── 3. Scale: f_L = ln(tau) / sqrt(L) ─────────────────────────────────────────
floor = 1e-14
f7 = np.log(np.clip(tau_7, floor, None)) / np.sqrt(L7)
f8 = np.log(np.clip(tau_8, floor, None)) / np.sqrt(L8)
 
# ── 4. Signed difference field: > 0 where L8 scaling exceeds L7, < 0 where it
#      falls below. The phase border is the zero-crossing of this field.
Delta = f8 - f7  # shape (Nr, N_phi)
 
Delta = gaussian_filter(Delta, sigma=1)

# ---------------------------------------------
# Plot
# ---------------------------------------------
fig = plt.figure(figsize=(8,7))
ax = plt.subplot(projection='polar')

# Filled contour of the difference
c = ax.contourf(
    Phi,
    R,
    Delta,
    levels=50,
    cmap='seismic'
)

# Zero contour = estimated phase boundary
boundary = ax.contour(
    Phi,
    R,
    Delta,
    levels=[0],
    colors='black',
    linewidths=3
)

plt.colorbar(c, pad=0.1, label=r'$\Delta = \ln(\tau)/\sqrt{L}|_{L=8}-\ln(\tau)/\sqrt{L}|_{L=7}$')

ax.set_title("Estimated Phase Boundary")

plt.show()
