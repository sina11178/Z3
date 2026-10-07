
import matplotlib.ticker as ticker
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
phi_line_2 = 1.25  # radians
phi_line_3 = np.pi/6  # radians
phi_line_4 = np.pi/3  # radians

print(phi_line_4)

#print(phi_line)
#ax.plot([phi_line_1, phi_line_1], [0, 1], linestyle='--', color='green', linewidth=2, label=r'$\phi=pi/4$ rad')
#ax.plot([phi_line_2, phi_line_2], [0, 1], linestyle='--', color='red', linewidth=2, label=r'$\phi=0.9*pi/4$ rad')
#ax.plot([phi_line_3, phi_line_3], [0, 1], linestyle='--', color='blue', linewidth=2, label=r'$\phi=pi/6$ rad')
#ax.plot([phi_line_4, phi_line_4], [0, 1], linestyle='--', color='black', linewidth=2, label=r'$\phi=pi/3$ rad')
#ax.legend(loc='upper right')
f_point_1 = 4
phi_point_1 = 10

f_point_2 = 3
phi_point_2 = 16

f_point_3 = 8
phi_point_3 = 18

r_1 = r_vals[f_point_1]
phi_1 = phi_vals[phi_point_1]

r_2 = r_vals[f_point_2]
phi_2 = phi_vals[phi_point_2]

r_3 = r_vals[f_point_3]
phi_3 = phi_vals[phi_point_3]

ax.set_rlim(0, 0.4)
ax.scatter(phi_1, r_1, color="#14D714", s=30,zorder=10)
ax.scatter(phi_2, r_2, color='blue', s=30,zorder=10)
ax.scatter(phi_3, r_3, color='purple', s=30,zorder=10)

ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
cbar = plt.colorbar(c, pad=0.1)
cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
angles_deg = np.arange(0, 360, 30)  # your current tick spacing
ax.set_thetagrids(angles_deg, labels=[
    '0', r'$\pi/6$', r'$\pi/3$', r'$\pi/2$', r'$2\pi/3$', r'$5\pi/6$',
    r'$\pi$', r'$7\pi/6$', r'$4\pi/3$', r'$3\pi/2$', r'$5\pi/3$', r'$11\pi/6$'
])

#ax.set_title("Estimated Phase Boundary")

plt.savefig('phase_boundary_estimation.png', dpi=600, transparent=True)
plt.savefig('phase_boundary_estimation.pdf', transparent=True)
plt.savefig('phase_boundary_estimation.svg', transparent=True)

plt.show()
# import numpy as np
# import matplotlib.pyplot as plt
# from skimage import measure
# from scipy.ndimage import gaussian_filter, map_coordinates
 
# # ── 1. Grid — must match the heatmap script exactly ──────────────────────────
# L7, L8 = 7, 8
# Nr = 40
# N_phi = 40
# r_vals = np.linspace(0, 1, Nr)
# phi_vals = np.linspace(0, 2 * np.pi / 3, N_phi)
 
# # ── 2. Load data — both r-major, shape (Nr, N_phi) ────────────────────────────
# tau_7 = 1 / np.load('lifetimes_L7_gamma_NPY.npy')   # shape (Nr, N_phi)
# tau_8 = 1 / np.load('lifetimes_L8_gamma_NPY.npy')   # shape (Nr, N_phi)
 
# assert tau_7.shape == (Nr, N_phi), f"tau_7 shape {tau_7.shape} != ({Nr},{N_phi})"
# assert tau_8.shape == (Nr, N_phi), f"tau_8 shape {tau_8.shape} != ({Nr},{N_phi})"
 
# # ── 3. Scale: f_L = ln(tau) / sqrt(L) ─────────────────────────────────────────
# floor = 1e-14
# f7 = np.log(np.clip(tau_7, floor, None)) / np.sqrt(L7)
# f8 = np.log(np.clip(tau_8, floor, None)) / np.sqrt(L8)
 
# # ── 4. Signed difference field: > 0 where L8 scaling exceeds L7, < 0 where it
# #      falls below. The phase border is the zero-crossing of this field.
# diff_raw = f8 - f7   # shape (Nr, N_phi)
 
# # Where |diff| is essentially 0 (numerical noise, not a real crossing), the raw
# # field flips sign pixel-to-pixel and find_contours traces every one of those
# # flips. A light Gaussian smooth suppresses single-pixel noise while leaving a
# # real, spatially-extended boundary intact.
# smooth_sigma = 1.0   # increase if you still see noisy contours; decrease if the
#                       # real boundary starts getting blurred/shifted
# diff = gaussian_filter(diff_raw, sigma=smooth_sigma)
 
# # ── 5. Contour the zero-crossing directly (sub-pixel accurate, no binarizing)─
# raw_contours = measure.find_contours(diff, level=1e-14)
 
# # ── 5b. Drop small contours: a real phase boundary spans many grid cells, a
# #      noise artifact from near-zero data is usually a tiny loop or wisp.
# min_span_pixels = 5   # a contour must span at least this many pixels in
#                        # either row or column direction to be kept; tune to taste
# def contour_span(contour):
#     row_span = contour[:, 0].max() - contour[:, 0].min()
#     col_span = contour[:, 1].max() - contour[:, 1].min()
#     return max(row_span, col_span)
 
# contours = [c for c in raw_contours if contour_span(c) >= min_span_pixels]
# print(f"Kept {len(contours)} of {len(raw_contours)} contours after filtering")
 
# def contour_to_polar(contour):
#     """(row, col) index-space contour points -> (phi, r) data coords."""
#     r_idx = contour[:, 0]
#     phi_idx = contour[:, 1]
#     r_c = np.interp(r_idx, np.arange(Nr), r_vals)
#     phi_c = np.interp(phi_idx, np.arange(N_phi), phi_vals)
#     return phi_c, r_c
 
# # ── 5c. Directional filtering: keep only crossings where, moving in the
# #      direction of INCREASING r, diff goes from >=0 down to <0 (i.e. drops
# #      from positive into zero-or-negative). This discards the reverse
# #      crossing (negative -> positive), which is what produces re-entrant
# #      loops like the "petal" seen around a locally positive island, OR
# #      conversely discards the main arc if that one happens to run the other
# #      way — check which one you get and flip `keep_direction` if needed.
# keep_direction = 'pos_to_neg'   # 'pos_to_neg' or 'neg_to_pos'
# step = 0.5                      # sample offset (in grid-index units) used to
#                                  # probe diff just before/after the crossing
 
# def classify_contour_direction(contour, field, step=0.5, n_samples=10):
#     """Vote across several points along the contour on whether, moving toward
#     increasing row index (r), the field goes positive->negative or negative->positive."""
#     n = len(contour)
#     idxs = np.linspace(0, n - 1, min(n_samples, n)).astype(int)
#     votes = []
#     for k in idxs:
#         row, col = contour[k]
#         v_before = map_coordinates(field, [[row - step], [col]], order=1, mode='nearest')[0]
#         v_after = map_coordinates(field, [[row + step], [col]], order=1, mode='nearest')[0]
#         if v_before >= 0 and v_after < 0:
#             votes.append('pos_to_neg')
#         elif v_before < 0 and v_after >= 0:
#             votes.append('neg_to_pos')
#     if not votes:
#         return 'ambiguous'
#     return max(set(votes), key=votes.count)
 
# print("\nContour directions (as r increases):")
# directional_contours = []
# for i, c in enumerate(contours):
#     direction = classify_contour_direction(c, diff, step=step)
#     print(f"  [{i}] direction={direction}")
#     if direction == keep_direction:
#         directional_contours.append(c)
 
# print(f"\nKeeping {len(directional_contours)} contour(s) with direction '{keep_direction}'")
 
# # ── 5c. Characterize each surviving contour so we can pick out specific ones
# #      (e.g. the small "petal" island near r~0) rather than treating them all
# #      as one boundary.
# def contour_properties(contour):
#     phi_c, r_c = contour_to_polar(contour)
#     is_closed = np.hypot(r_c[0] - r_c[-1], phi_c[0] - phi_c[-1]) < 1e-6
#     return {
#         'r_min': r_c.min(), 'r_max': r_c.max(),
#         'phi_min': phi_c.min(), 'phi_max': phi_c.max(),
#         'n_points': len(r_c),
#         'closed': is_closed,
#     }
 
# print("\nContour summary:")
# for i, c in enumerate(contours):
#     p = contour_properties(c)
#     print(f"  [{i}] r=({p['r_min']:.3f},{p['r_max']:.3f}) "
#           f"phi=({p['phi_min']:.3f},{p['phi_max']:.3f}) "
#           f"n={p['n_points']} closed={p['closed']}")
 
# # A closed loop confined to small r (and small phi) is the "petal" — pick it
# # out by its bounding box rather than by hardcoding an index, since the index
# # order isn't guaranteed to stay stable if the data changes.
# petal_r_max = 0.3   # tune: upper bound on r for what counts as "near origin"
# petal_indices = [i for i, c in enumerate(contours)
#                   if contour_properties(c)['closed']
#                   and contour_properties(c)['r_max'] < petal_r_max]
# print(f"\nCandidate petal contour(s): {petal_indices}")
 
# # ── 6. Plot ────────────────────────────────────────────────────────────────────
# fig, axes = plt.subplots(1, 2, figsize=(13, 6), subplot_kw={'projection': 'polar'})
 
# R, PHI = np.meshgrid(r_vals, phi_vals, indexing='ij')  # shape (Nr, N_phi)
 
# def style_polar_ax(ax, title):
#     ax.set_ylim(0, 1)
#     ax.set_thetamin(0)
#     ax.set_thetamax(120)
#     ax.set_title(title)
 
# # 6a. Signed difference field with boundary overlaid
# ax = axes[0]
# vmax = np.abs(diff).max()
# im = ax.pcolormesh(PHI, R, diff, shading='auto', cmap='RdBu', vmin=-vmax, vmax=vmax)
# fig.colorbar(im, ax=ax, label=r'$f_8 - f_7$')
# # all zero-crossings, faint, for reference
# for contour in contours:
#     phi_c, r_c = contour_to_polar(contour)
#     ax.plot(phi_c, r_c, color='gray', linewidth=1, alpha=0.5)
# # only the crossings matching keep_direction, highlighted
# for contour in directional_contours:
#     phi_c, r_c = contour_to_polar(contour)
#     ax.plot(phi_c, r_c, 'k-', linewidth=2.5)
# style_polar_ax(ax, f"$f_8 - f_7$, direction='{keep_direction}' highlighted")
 
# # 6a-ii. Highlight the petal contour(s) in a distinct color on the same axis
# for i in petal_indices:
#     phi_c, r_c = contour_to_polar(contours[i])
#     ax.plot(phi_c, r_c, color='magenta', linewidth=2.5, zorder=5, linestyle='--')
 
# # 6b. L=8 decay rate with directionally-filtered boundary overlaid
# ax = axes[1]
# decay_8 = 1.0 / np.clip(tau_8, floor, None)
# im = ax.pcolormesh(PHI, R, decay_8, shading='auto', cmap='Blues_r')
# fig.colorbar(im, ax=ax, label=r'Decay Rate $\lambda$, $L=8$')
# for contour in directional_contours:
#     phi_c, r_c = contour_to_polar(contour)
#     ax.plot(phi_c, r_c, 'y-', linewidth=2.5)
# style_polar_ax(ax, f'Boundary overlaid on decay rate, $L={L8}$')
 
# plt.tight_layout()
# plt.savefig('phase_border_crossing.png', dpi=150)
# plt.show()

# '''
# import numpy as np
# import matplotlib.pyplot as plt
# from skimage import measure
 
# # ── 1. Grid — must match the heatmap script exactly ──────────────────────────
# L7, L8 = 7, 8
# Nr = 40
# N_phi = 40
# r_vals = np.linspace(0, 1, Nr)
# phi_vals = np.linspace(0, 2 * np.pi / 3, N_phi)
 
# # ── 2. Load data — both r-major, shape (Nr, N_phi) ────────────────────────────
# tau_7 = 1 / np.load('lifetimes_L7_gamma_NPY.npy')   # shape (Nr, N_phi)
# tau_8 = 1 / np.load('lifetimes_L8_gamma_NPY.npy')   # shape (Nr, N_phi)
 
# assert tau_7.shape == (Nr, N_phi), f"tau_7 shape {tau_7.shape} != ({Nr},{N_phi})"
# assert tau_8.shape == (Nr, N_phi), f"tau_8 shape {tau_8.shape} != ({Nr},{N_phi})"
 
# # ── 3. Scale: f_L = ln(tau) / sqrt(L) ─────────────────────────────────────────
# floor = 1e-14
# f7 = np.log(np.clip(tau_7, floor, None)) / np.sqrt(L7)
# f8 = np.log(np.clip(tau_8, floor, None)) / np.sqrt(L8)
 
# # ── 4. Signed difference field: > 0 where L8 scaling exceeds L7, < 0 where it
# #      falls below. The phase border is the zero-crossing of this field.
# diff = f8 - f7   # shape (Nr, N_phi)
 
# # ── 5. Contour the zero-crossing directly (sub-pixel accurate, no binarizing)─
# contours = measure.find_contours(diff, level=0.0)
 
# def contour_to_polar(contour):
#     """(row, col) index-space contour points -> (phi, r) data coords."""
#     r_idx = contour[:, 0]
#     phi_idx = contour[:, 1]
#     r_c = np.interp(r_idx, np.arange(Nr), r_vals)
#     phi_c = np.interp(phi_idx, np.arange(N_phi), phi_vals)
#     return phi_c, r_c
 
# # ── 6. Plot ────────────────────────────────────────────────────────────────────
# fig, axes = plt.subplots(1, 2, figsize=(13, 6), subplot_kw={'projection': 'polar'})
 
# R, PHI = np.meshgrid(r_vals, phi_vals, indexing='ij')  # shape (Nr, N_phi)
 
# def style_polar_ax(ax, title):
#     ax.set_ylim(0, 1)
#     ax.set_thetamin(0)
#     ax.set_thetamax(120)
#     ax.set_title(title)
 
# # 6a. Signed difference field with boundary overlaid
# ax = axes[0]
# vmax = np.abs(diff).max()
# im = ax.pcolormesh(PHI, R, diff, shading='auto', cmap='RdBu', vmin=-vmax, vmax=vmax)
# fig.colorbar(im, ax=ax, label=r'$f_8 - f_7$')
# for contour in contours:
#     phi_c, r_c = contour_to_polar(contour)
#     ax.plot(phi_c, r_c, 'k-', linewidth=2)
# style_polar_ax(ax, r'$f_8 - f_7$ with phase border')
 
# # 6b. L=8 decay rate with boundary overlaid
# ax = axes[1]
# decay_8 = 1.0 / np.clip(tau_8, floor, None)
# im = ax.pcolormesh(PHI, R, decay_8, shading='auto', cmap='Blues_r')
# fig.colorbar(im, ax=ax, label=r'Decay Rate $\lambda$, $L=8$')
# for contour in contours:
#     phi_c, r_c = contour_to_polar(contour)
#     ax.plot(phi_c, r_c, 'y-', linewidth=2)
# style_polar_ax(ax, f'Boundary overlaid on decay rate, $L={L8}$')
 
# plt.tight_layout()
# plt.savefig('phase_border_crossing.png', dpi=150)
# plt.show()
# '''