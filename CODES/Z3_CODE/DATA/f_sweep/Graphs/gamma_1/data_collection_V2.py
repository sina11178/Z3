import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap

f_point = 9


f_array = np.linspace(0, 0.5, 11)
phi_array = [0, np.pi/6]
gamma_array = [1]
L_array = [4 ,5, 6, 7, 8]

time_list = np.exp(np.arange(-10, 25, 0.5))

plt.figure()
data_sets = []
for l in L_array:
    index = 1
    params = [f_array[f_point], phi_array[1]]
    if os.path.isfile(f'../../data/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy') == True:
        data = np.load(f'../../data/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]
        data_sets.append(data)


    else:
        print(f"File not found for pixel")

log_time_list = np.arange(-10, 25, 0.5)

for i, l in enumerate(L_array):
    plt.plot(time_list, data_sets[i], label = f"L = {l}")

plt.xscale('log')
plt.xlabel('ln(t)')
plt.ylabel("Autocorr")
plt.legend()

plt.title(f"Gamma = 1, f = {f_array[f_point]}, φ = {phi_array[1]}")
plt.tight_layout()
plt.show()
    



#print(phase_diagram)