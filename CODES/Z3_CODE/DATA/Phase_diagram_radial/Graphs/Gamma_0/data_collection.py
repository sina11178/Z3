import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit



def make_p0(t, y):
    B0 = np.median(y[-max(5, len(y)//10):])       # tail plateau estimate
    A0 = y[0] - B0                                  # initial drop amplitude
    A0 = max(A0, 1e-3)                              # avoid zero amplitude

    # estimate the timescale from where y crosses the midpoint
    mid = B0 + A0/2
    idx = np.argmin(np.abs(y - mid))
    t_half = t[idx]
    gamma0 = 1/t_half if t_half > 0 else 1e-6

    return [np.clip(A0, 0, 1), np.clip(B0, 0, 0.5), np.clip(gamma0, -1, 1)]




def exp_model(t, A, B, gamma):
    return B + A*np.exp(-t*gamma)


def func(r, phi, time_list, l, cutoff):
    if os.path.isfile(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy') == True:
        data = np.load(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy')[1]
        try:
    # Fit on log-scale time for numerical stability
            time_list_new = time_list[cutoff:]


            data_2 = np.load(f'../../extra_time_points/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy')[1]
            data_3 = np.load(f'../../extra_time_points_V2/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy')[1]
           
            

            data = np.concatenate((data, data_2))
            data = np.concatenate((data, data_3))

            added_time = np.exp(np.arange(25, 40, 0.5))
            #added_time = np.exp(np.arange(25, 30, 0.5))
            time_list_new = np.concatenate((time_list_new, added_time))
            
            y = data[cutoff_index:]

            mask = np.isfinite(time_list_new) & np.isfinite(y)
            time_list_new = time_list_new[mask]
            y = y[mask]

            p0 = make_p0(time_list_new, y)
            popt, _ = curve_fit(
                exp_model,
                time_list_new,
                y,
                p0=p0, #[data[-1], 0, 1, (1-data[-1]), 1 / time_list[len(time_list)//2]],#[data_current[-1], data_current[0] - data_current[-1], 1 / time_list_new[len(time_list_new)//2]],
                #bounds = ([0, 0, -2], [1, 0.5, 10]),
                maxfev=10000
            )
            fit_curve = exp_model(time_list_new, *popt)
            A, B, gamma = popt
            return gamma
        except RuntimeError:
            print(f"Fit did not converge for L={l}")



gamma_array = [0]
L_array = [6, 7, 8]

Nr = 40
N_phi = 40

r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2*np.pi/3, N_phi)

R, PHI = np.meshgrid(r_vals, phi_vals, indexing='ij')

time_list = np.exp(np.arange(-10, 25, 0.5))
#time_list = np.exp(np.arange(-10, 30, 0.5))

cutoff_index = 0

for l in L_array:
    Z = np.zeros((Nr, N_phi))
    for i, r in enumerate(r_vals):
        for j, phi in enumerate(phi_vals):
            Z[i, j] = func(r, phi, time_list, l, cutoff_index)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    c = ax.pcolormesh(PHI, R, Z, shading='auto', cmap='viridis')

    fig.colorbar(c, ax=ax, label='Decay Rate γ')
    np.save(f'lifetimes_L{l}_gamma_NPY.npy', Z)

    ax.set_ylim(0, 1)
    ax.set_thetamin(0)
    ax.set_thetamax(120)   # pi/3 in degrees
    ax.set_title(f'Polar plot of K for L={l}')
    plt.savefig(f'phase_diagram_L{l}.png', dpi=150)

    plt.show()