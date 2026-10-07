import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
from matplotlib import colormaps
from matplotlib.lines import Line2D   # put this with your imports


f_point = 4
phi_point = 10
phi_point_2 = 29

Nr = 40
N_phi = 40
r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2 * np.pi / 3, N_phi)
gamma_array = [0.5]
L_array = [4, 5, 6, 7, 8]

colors = ["#12CBC4", "#009432", "#0652DD", "#FFC312", "#EA2027"]  # colorblind-friendly
colors = colormaps['plasma'].resampled(len(L_array)+1).colors


phi_val_check = phi_vals[phi_point_2]

time_list = np.exp(np.arange(-10, 25, 0.5))

time_list_45 = np.exp(np.arange(-10, 45, 0.5))

#plt.figure()
plt.figure(figsize=(4, 3))
data_sets = []
data_set_curve = []
cutoff = 0
for l in L_array:
    index = 1
    params = [r_vals[f_point], phi_vals[phi_point]]
    if os.path.isfile(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy') == True:
        data = np.load(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]

        if l != 4 and l != 5:
            data_2 = np.load(f'../../extra_time_points/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]
            data_3 = np.load(f'../../extra_time_points_V2/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]

            data = np.concatenate((data, data_2))
            data = np.concatenate((data, data_3))
            data_sets.append(data)
        else:
            data_sets.append(data)
        #data_set_curve.append(func(params[0], params[1], time_list, cutoff, l))


    else:
        print(f"File not found for pixel")

log_time_list = np.arange(-10, 25, 0.5)

for i, l in enumerate(L_array):
    if l != 4 or l != 5:
        added_time_exp = np.exp(np.arange(25, 40, 0.5))
        time_list_78 = np.concatenate((time_list, added_time_exp))
        plt.plot(time_list_78, data_sets[i], color = colors[i], label = f"L = {l}", linewidth=2.5)
        #plt.plot(time_list_78, data_set_curve[i], linestyle='--', label = f"L = {l}, approximate curve")
    else:
        plt.plot(time_list_45, data_sets[i], color = colors[i], label = f"L = {l}", linewidth=2.5)
        #plt.plot(time_list, data_set_curve[i], linestyle='--', label = f"L = {l}, approximate curve")


data_sets = []
data_set_curve = []
cutoff = 0
for l in L_array:
    index = 1
    params = [r_vals[f_point], phi_vals[phi_point_2]]
    if os.path.isfile(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy') == True:
        data = np.load(f'../../L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]

        if l != 4 and l != 5:
            data_2 = np.load(f'../../extra_time_points/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]
            data_3 = np.load(f'../../extra_time_points_V2/L{l}/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]

            data = np.concatenate((data, data_2))
            data = np.concatenate((data, data_3))
            data_sets.append(data)
        else:
            data_sets.append(data)
        #data_set_curve.append(func(params[0], params[1], time_list, cutoff, l))


    else:
        print(f"File not found for pixel")

log_time_list = np.arange(-10, 25, 0.5)

for i, l in enumerate(L_array):
    if l != 4 or l != 5:
        added_time_exp = np.exp(np.arange(25, 40, 0.5))
        time_list_78 = np.concatenate((time_list, added_time_exp))
        plt.plot(time_list_78, data_sets[i], color = colors[i], ls = '--', linewidth=2.5)
        #plt.plot(time_list_78, data_set_curve[i], linestyle='--', label = f"L = {l}, approximate curve")
    else:
        plt.plot(time_list_45, data_sets[i], color = colors[i], ls = '--', linewidth=2.5)
        #plt.plot(time_list, data_set_curve[i], linestyle='--', label = f"L = {l}, approximate curve")





  # instead of ax = axes[i, j]
#plt.figure(figsize=(14, 5))
plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
plt.setp(plt.gca().spines.values(), linewidth=2.5)
plt.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator())
plt.gca().yaxis.set_minor_locator(ticker.AutoMinorLocator())

plt.tick_params(axis="x", direction="in", which='both')
plt.tick_params(axis="y", direction="in", which='both')
plt.tick_params(axis="x", which='major', length=8.0, width=2)
plt.tick_params(axis="x", which='minor', length=5.0, width=2)
plt.tick_params(axis="y", which='major', length=8.0, width=2)
plt.tick_params(axis="y", which='minor', length=5.0, width=2)

plt.tick_params(axis="x", labelsize=18)
plt.tick_params(axis="y", labelsize=18)

plt.ylim(-0, 1.05)
plt.yticks([ 0, 0.5, 1], ['', '', ''])


plt.xscale('log')
plt.xticks([1e-3, 1e3, 1e9, 1e15])#, labels=['', '', '', ''])
#plt.xlabel('ln(t)')
#plt.ylabel("Autocorr")
#plt.legend()
#plt.legend()

#plt.title(f"Gamma = 1, f = {r_vals[f_point]}, φ = {phi_vals[phi_point]}")

plt.savefig('pixel_plot_symm_check.png', dpi=600, transparent=True)
plt.savefig('pixel_plot_symm_check.pdf', transparent=True)
plt.savefig('pixel_plot_symm_check.svg', transparent=True)
plt.show()
    



#print(phase_diagram)