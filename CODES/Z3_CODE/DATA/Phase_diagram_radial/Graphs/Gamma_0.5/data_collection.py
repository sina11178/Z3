


import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit



def exp_model(t, A, B, gamma):
    return B + A*np.exp(-t*gamma)


def func(r, phi, time_list, l):
    if os.path.isfile(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy') == True:
        data = np.load(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy')[1]
        try:
    # Fit on log-scale time for numerical stability
            popt, _ = curve_fit(
                exp_model,
                time_list_new,
                data[cutoff_index:],
                #p0= [data[-1], 0, 1, (1-data[-1]), 1 / time_list[len(time_list)//2]],#[data_current[-1], data_current[0] - data_current[-1], 1 / time_list_new[len(time_list_new)//2]],
                bounds = ([0, 0, 0], [1, 0.5, 10]),
                maxfev=10000
            )
            fit_curve = exp_model(time_list, *popt)
            A, B, gamma = popt
            return gamma
        except RuntimeError:
            print(f"Fit did not converge for L={l}")



gamma_array = [0.5]
L_array = [7, 8]

Nr = 40
N_phi = 40

r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2*np.pi/3, N_phi)

R, PHI = np.meshgrid(r_vals, phi_vals, indexing='ij')

time_list = np.exp(np.arange(-10, 25, 0.5))

cutoff_index = 0
time_list_new = time_list[cutoff_index:]

for l in L_array:
    Z = np.zeros((Nr, N_phi))
    for i, r in enumerate(r_vals):
        for j, phi in enumerate(phi_vals):
            Z[i, j] = func(r, phi, time_list_new, l)
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
