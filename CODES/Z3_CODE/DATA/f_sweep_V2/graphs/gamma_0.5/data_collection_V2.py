import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit

f_point =3


f_array = np.linspace(0, 1, 21)
phi_array = [0, np.pi/6]
gamma_array = [0]
L_array = [4, 5, 6, 7, 8]
#L_array = [7]
time_list = np.exp(np.arange(-9, 27, 0.5))

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

#log_time_list = np.arange(-9, 27, 0.5)

def exp_model(t, A, B, beta):
    return A + B * np.exp(-t * beta)

#def exp_model_new(t, c, a, g_1, g_0):
#    out_coeff = 1/(1+a)
#    initial_decay = a + np.exp(-t * g_1)

#    return out_coeff * initial_decay * np.exp(-t * g_0) + c

def exp_model_new(t, c, a_1, g_1, a_2, g_2):
    term_1 = a_1 * np.exp(-t*g_1)
    term_2 = a_2 * np.exp(-t*g_2)

    return c + term_1 + term_2

def detect_plateau_end(time_list_new, data_current, flatness_threshold=.1):
    """Returns the index where the data starts dropping significantly."""
    # Normalize data to [0,1] range for threshold to be meaningful
    #data_norm = (data_current - data_current.min()) / (data_current.max() - data_current.min())

    slopes = []
    
    for j in range(1, len(data_current)):
        diff_data = data_current[j] - data_current[j-1]
        diff_time_log = np.log(time_list_new[j]) - np.log(time_list_new[j-1])
        slope = abs(diff_data/diff_time_log)
        slopes.append(slope)
    print(slopes)
    binary_slopes = np.array(slopes) > flatness_threshold
    print(binary_slopes)

    transitions_off = np.where(np.diff(binary_slopes.astype(int)) == 1)[0]
    transitions = np.where(np.diff(binary_slopes.astype(int)) == -1)[0]
    # np.diff == 1 means went from 0 (False) to 1 (True)
    # the index returned is the last False before the transition

    num_transitions = len(transitions)
    print(f"Number of True->False transitions: {num_transitions}")
    print(f"Indices of last True before each transition: {transitions}")

    true_regions = []
    for start, end in zip(transitions, transitions_off):
        length = -end + start
        if length >= 3:
            true_regions.append((start, end, length))

    print(f"Number of significant transitions: {len(true_regions)}")
    for start, end, length in true_regions:
        print(f"True region from index {start+1} to {end}, length {length}")


    if len(true_regions) == 1:
        return 0  # 2 plateaus, no truncation
    elif len(true_regions) >= 2:
        return true_regions[1][0]  # cut before second decay
    else:
        raise ValueError('Houston, we have a problem...')
    #    if slope > flatness_threshold:
    #        print(j-1)
    #        return max(0, j-1)  # return index just before the drop
    #return 0  # no plateau found

cutoff_index = 16
time_list_new = time_list[cutoff_index:]

for i, l in enumerate(L_array):
    data_current = data_sets[i]
    data_current = data_current[cutoff_index:]
    color = plt.gca()._get_lines.get_next_color()
    plt.plot(time_list, data_sets[i], color = color, label = f"L = {l}")

    '''
    try:
        plateau_end = detect_plateau_end(time_list_new, data_current)

        if plateau_end > 2:  # meaningful plateau exists
            fit_time = time_list_new[plateau_end:]
            fit_data = data_current[plateau_end:]
        else:  # no plateau, fit everything
            fit_time = time_list_new
            fit_data = data_current

        popt, _ = curve_fit(
            exp_model,
            fit_time,
            fit_data,
            p0=[fit_data[-1], fit_data[0] - fit_data[-1], 1 / fit_time[len(fit_time)//2]],
            maxfev=100000
        )
        fit_curve = exp_model(time_list, *popt)
        A, B, beta = popt
        plt.plot(time_list, fit_curve, color=color, linestyle='--',
                 label=f"L={l} fit (A={A:.2f}, B={B:.2f}, 1/β={1/beta:.2e})")
    except RuntimeError:
        print(f"Fit did not converge for L={l}")

        
    '''
    
    try:
        # Fit on log-scale time for numerical stability
        popt, _ = curve_fit(
            exp_model_new,
            time_list_new,
            data_current,
            #p0= [data_current[-1], 0, 1, (1-data_current[-1]), 1 / time_list_new[len(time_list_new)//2]],#[data_current[-1], data_current[0] - data_current[-1], 1 / time_list_new[len(time_list_new)//2]],
            bounds = ([0, 0, 0, 0, 0], [.5, 1, 10, 1, 10]),
            maxfev=10000
        )
        fit_curve = exp_model_new(time_list, *popt)
        c, a_1, g_1, a_2, g_2,  = popt
        plt.plot(time_list, fit_curve, color=color, linestyle='--',
                    label=f"L={l} fit (a_1={a_1:.2f}, γ_1 ={g_1:.2f}, a_2 = {a_2:.2f}, γ_2={g_2:.2e}, c = {c: .2f})")
    except RuntimeError:
        print(f"Fit did not converge for L={l}")
    

plt.xscale('log')
plt.xlabel('ln(t)')
plt.ylabel("Autocorr")
plt.legend()

plt.title(f"Gamma = 1, f = {f_array[f_point]}, φ = {phi_array[1]}")
plt.tight_layout()
plt.show()
    



#print(phase_diagram)