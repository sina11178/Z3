import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import binary_fill_holes
from skimage import measure

# ── 1. Load data ──────────────────────────────────────────────────────────────
L1, L2, L3 = 4, 6, 8





tau1 = np.load('lifetimes_L4_gamma_NPY.npy')   # shape (N_phi, N_r) — origin top-left
tau2 = np.load('lifetimes_L6_gamma_NPY.npy')
tau3 = np.load('lifetimes_L8_gamma_NPY.npy')

if tau1[0,0] == np.inf:
    tau1[0,0] = np.exp(25*np.sqrt(L1))
if tau2[0,0] == np.inf:
    tau2[0,0] = np.exp(25*np.sqrt(L2)) 
if tau3[0,0] == np.inf:
    tau3[0,0] = np.exp(25*np.sqrt(L3))

# Flip so origin is bottom-left (standard plot orientation)
#tau1, tau2, tau3 = np.flipud(tau1), np.flipud(tau2), np.flipud(tau3)

# ── 2. Compute log(τ)/L for each size ─────────────────────────────────────────
# Clip to avoid log(0) issues; tweak the floor (1e-10) if needed
min_scaling = 0.001  # Set a minimum scaling value to identify in-phase region
floor = 1e-14
f1 = np.log(np.clip(tau1, floor, None)) / np.sqrt(L1)
f2 = np.log(np.clip(tau2, floor, None)) / np.sqrt(L2)
f3 = np.log(np.clip(tau3, floor, None)) / np.sqrt(L3)

print("Lifetime of (0,0) for L4: " + str(f1[0,0]))
print("Lifetime of (0,0) for L6: " + str(f2[0,0]))
print("Lifetime of (0,0) for L8: " + str(f3[0,0]))

# ── 3. Classify: in-phase where f is flat (not decaying) across L ─────────────
# "Flat" = f3 ≥ f2 and f2 ≥ f1*(1 - tol)
# Increase tol (0.0 → 0.2) to be more permissive near the boundary
tol = 0
#ratio_32 = f3 / np.clip(f2, floor, None)  # f3/f2 ≈ 1 in phase
#ratio_21 = f2 / np.clip(f1, floor, None)  # f2/f1 ≈ 1 in phase
# in_phase_raw = (np.abs(f3/f2)<10) & (np.abs(f3/f2)>0.7) & (np.abs(f2/f1)<10) & (np.abs(f2/f1)>0.7) & (f3 >= min_scaling)  # Also require f3 above a minimum to exclude very short lifetimes
in_phase_raw = (np.abs(f3/f2)>=1) & (np.abs(f2/f1)>=1) #& (f3 >= min_scaling)  # Also require f3 above a minimum to exclude very short lifetimes)
#in_phase_raw =(f3 >= f2*(1-tol)) & (f2 >= f1 * (1 - tol))  & (f3 >= min_scaling)  # Also require f3 above a minimum to exclude very short lifetimes
#in_phase_raw = (np.abs(ratio_32 - 1) < tol) & (np.abs(ratio_21 - 1) < tol) & (f3 >= min_scaling)  # Also require f3 above a minimum to exclude very short lifetimes
# ── 4. Clean up: fill holes, remove isolated pixels ──────────────────────────
from scipy.ndimage import binary_opening, binary_closing

struct = np.ones((1, 1))                        # 2×2 structuring element
in_phase_clean = binary_closing(in_phase_raw, structure=struct)   # close gaps
in_phase_clean = binary_opening(in_phase_clean, structure=struct)  # remove specks
in_phase_clean = binary_fill_holes(in_phase_clean)                 # fill interior

binary_matrix = in_phase_clean.astype(int)
#binary_matrix = in_phase_raw.astype(int)  # Use raw classification without morphological operations
# ── 5. Extract boundary contour ───────────────────────────────────────────────
# Build your axis grids to match the plot (adjust these to your actual grid)
N_x, N_y = tau1.shape
x_vals = np.linspace(0, 1, N_x)   # f/J * cos(φ) axis
y_vals = np.linspace(0, 1, N_y)   # f/J * sin(φ) axis

contours = measure.find_contours(binary_matrix, level=0.5)

# ── 6. Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# — 6a. Raw log(τ)/L for largest system —
ax = axes[0]
im = ax.imshow(f3, origin='lower', aspect='auto',
               extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]],
               cmap='Blues_r')
plt.colorbar(im, ax=ax, label=r'$\log(\tau)/L$')
ax.set_title(f'$\\log(\\tau)/L$,  $L={L3}$')
ax.set_xlabel(r'$f/J \cdot \cos(\varphi)$')
ax.set_ylabel(r'$f/J \cdot \sin(\varphi)$')

# — 6b. Binary phase matrix —
ax = axes[1]
ax.imshow(binary_matrix, origin='lower', aspect='auto',
          extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]],
          cmap='Blues_r', vmin=0, vmax=1)
ax.set_title('Binary phase matrix')
ax.set_xlabel(r'$f/J \cdot \cos(\varphi)$')
ax.set_ylabel(r'$f/J \cdot \sin(\varphi)$')

# — 6c. Boundary overlaid on L3 decay rate —
ax = axes[2]
# Load or recompute decay rate (= 1/tau); reuse tau3 here
decay3 = 1.0 / np.clip(tau3, floor, None)
im = ax.imshow(decay3, origin='lower', aspect='auto',
               extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]],
               cmap='Blues_r')
plt.colorbar(im, ax=ax, label=r'Decay Rate $\lambda$')

# Convert contour pixel coords → axis coords
for contour in contours:
    # contour[:,0] = row index (phi), contour[:,1] = col index (r)
    r_contour   = x_vals[0]   + contour[:, 1] / (N_x   - 1) * (x_vals[-1]   - x_vals[0])
    y_contour = y_vals[0] + contour[:, 0] / (N_y - 1) * (y_vals[-1] - y_vals[0])
    ax.plot(r_contour, y_contour, 'y-', linewidth=2, label='Phase boundary')

ax.set_title(f'Boundary overlaid on decay rate  $L={L3}$')
ax.set_xlabel(r'$f/J \cdot \cos(\varphi)$')
ax.set_ylabel(r'$f/J \cdot \sin(\varphi)$')
ax.legend()

plt.tight_layout()
plt.savefig('phase_boundary.png', dpi=150)
plt.show()

# ── 7. Save outputs ───────────────────────────────────────────────────────────
np.save('binary_phase_matrix.npy', binary_matrix)

# Save boundary coordinates as (r, phi) pairs
for i, contour in enumerate(contours):
    r_c   = x_vals[0]   + contour[:, 1] / (N_x   - 1) * (x_vals[-1]   - x_vals[0])
    y_c = y_vals[0] + contour[:, 0] / (N_y - 1) * (y_vals[-1] - y_vals[0])
    np.savetxt(f'boundary_contour_{i}.csv',
               np.column_stack([r_c, y_c]),
               header='r, y', delimiter=',')

print("Binary matrix shape:", binary_matrix.shape)
print("Number of boundary contours found:", len(contours))