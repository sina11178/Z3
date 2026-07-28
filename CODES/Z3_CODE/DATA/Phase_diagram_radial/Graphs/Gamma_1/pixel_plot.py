import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit

f_point = 2
phi_point = 14

def exp_model(t, A, B, gamma):
    return B + A*np.exp(-t*gamma)

def func(r, phi, time_list, cutoff_index, l):
    if os.path.isfile(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy') == True:
        data = np.load(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(r,3)}_phi_{round(phi,3)}_NPY.npy')[1]
        try:
    # Fit on log-scale time for numerical stability
            time_list_new = time_list[cutoff_index:]
            popt, _ = curve_fit(
                exp_model,
                time_list_new,
                data[cutoff_index:],
                #p0= [data[-1], 0, 1, (1-data[-1]), 1 / time_list[len(time_list)//2]],#[data_current[-1], data_current[0] - data_current[-1], 1 / time_list_new[len(time_list_new)//2]],
                bounds = ([0, 0, -1], [1, 0.5, 1]),
                maxfev=10000
            )
            fit_curve = exp_model(time_list, *popt)
            A, B, gamma = popt
            return fit_curve
        except RuntimeError:
            print(f"Fit did not converge for L={l}")

Nr = 40
N_phi = 40
r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2 * np.pi / 3, N_phi)
gamma_array = [1]
L_array = [6, 7, 8]

time_list = np.exp(np.arange(-10, 25, 0.5))

plt.figure()
data_sets = []
data_set_curve = []
cutoff = 26
for l in L_array:
    index = 1
    params = [r_vals[f_point], phi_vals[phi_point]]
    if os.path.isfile(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy') == True:
        data = np.load(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]
        data_sets.append(data)
        data_set_curve.append(func(params[0], params[1], time_list, cutoff, l))


    else:
        print(f"File not found for pixel")

log_time_list = np.arange(-10, 25, 0.5)

for i, l in enumerate(L_array):
    plt.plot(time_list, data_sets[i], marker= 'o', label = f"L = {l}")
    plt.plot(time_list, data_set_curve[i], linestyle='--', label = f"L = {l}, approximate curve")

plt.xscale('log')
plt.xlabel('ln(t)')
plt.ylabel("Autocorr")
plt.legend()

plt.title(f"Gamma = 1, f = {r_vals[f_point]}, φ = {phi_vals[phi_point]}")
plt.tight_layout()
plt.show()
    



#print(phase_diagram)