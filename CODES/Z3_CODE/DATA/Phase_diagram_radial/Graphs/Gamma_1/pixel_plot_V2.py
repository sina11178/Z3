import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit

f_point = 1
phi_point = 22

Nr = 40
N_phi = 40
r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2 * np.pi / 3, N_phi)
gamma_array = [1]
L_array = [4, 5, 6, 7, 8]

time_list = np.exp(np.arange(-10, 25, 0.5))

plt.figure()
data_sets = []
data_set_curve = []
cutoff = 0
for l in L_array:
    index = 1
    params = [r_vals[f_point], phi_vals[phi_point]]
    if os.path.isfile(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy') == True:
        data = np.load(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]

        data_2 = np.load(f'../../extra_time_points/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]
        data_3 = np.load(f'../../extra_time_points_V2/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]
        
        data = np.concatenate((data, data_2))
        data = np.concatenate((data, data_3))
        data_sets.append(data)
        #data_set_curve.append(func(params[0], params[1], time_list, cutoff, l))


    else:
        print(f"File not found for pixel")

log_time_list = np.arange(-10, 25, 0.5)

for i, l in enumerate(L_array):
    added_time_exp = np.exp(np.arange(25, 40, 0.5))
    time_list_78 = np.concatenate((time_list, added_time_exp))
    plt.plot(time_list_78, data_sets[i], label = f"L = {l}")
    #plt.plot(time_list_78, data_set_curve[i], linestyle='--', label = f"L = {l}, approximate curve")

        

plt.xscale('log')
plt.xlabel('ln(t)')
plt.ylabel("Autocorr")
plt.legend()

plt.title(f"Gamma = 1, f = {r_vals[f_point]}, φ = {phi_vals[phi_point]}")

plt.savefig('pixel_plot_V2.png', dpi=600, transparent=True)
plt.savefig('pixel_plot_V2.pdf', transparent=True)
plt.show()
    



#print(phase_diagram)