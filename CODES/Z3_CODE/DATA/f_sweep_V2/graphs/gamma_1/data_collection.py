import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit


f_array = np.linspace(0, 1, 21)
phi_array = [0, np.pi/6]

gamma_array = [1]
L_array = [4, 5, 6, 7, 8, 9]

time_list = np.exp(np.arange(-9, 27, 0.5))


def exp_model_new(t, c, a_1, g_1, a_2, g_2):
    term_1 = a_1 * np.exp(-t*g_1)
    term_2 = a_2 * np.exp(-t*g_2)

    return c + term_1 + term_2

def exp_model(t, A, B, gamma):
    return B + A*np.exp(-t*gamma)




fig, axes = plt.subplots(1, 2, figsize=(8, 6))

power = 1/10
cutoff_index = 16
time_list_new = time_list[cutoff_index:]

for l in L_array:
    scaled_lifetime_phi_0 = []
    scaled_lifetime_phi_pi_6 = []
    for j in range(42):
        index = j
        params = [f_array[index%21], phi_array[int(index/21)%2]]
        if os.path.isfile(f'../../data/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy') == True:
            data = np.load(f'../../data/gamma_{gamma_array[0]}/autocorr_gamma_{gamma_array[0]}_L{l}_f_{round(params[0],3)}_phi_{round(params[1],3)}_NPY.npy')[1]
            #if data[0] <0.01:
            #    phase_diagram[int(j / 31), int(j % 31)] = None
            #    continue
            #target = data[0] / np.e

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
            except RuntimeError:
                print(f"Fit did not converge for L={l}")

            
            tau = 1/gamma


            if int(index/21)%2 == 0:
                scaled_lifetime_phi_0.append(np.log(tau)/np.power(l, power))
            elif int(index/21)%2 == 1:
                scaled_lifetime_phi_pi_6.append(np.log(tau)/np.power(l, power))
            else:
                raise ValueError("Houston, we have a problem")


        else:
            print(f"File not found for pixel {j}")

    
    if l == 9:
        scaled_lifetime_phi_pi_6[0] = 27
        #f_array = np.delete(f_array, 6)
        #axes[0].plot(f_array[1:].pop(6), scaled_lifetime_phi_0[1:], label=f"L={l}")
        axes[1].plot(f_array[1:], scaled_lifetime_phi_pi_6[1:], label=f"L={l}")
    else:
        scaled_lifetime_phi_0[0] = 27
        scaled_lifetime_phi_pi_6[0] = 27
        axes[0].plot(f_array[1:], scaled_lifetime_phi_0[1:], label=f"L={l}")
        axes[1].plot(f_array[1:], scaled_lifetime_phi_pi_6[1:], label=f"L={l}")


axes[0].set_xlabel("f")
axes[0].set_ylabel(f"ln(tau)/L^{power}")
axes[0].set_title(f"φ = {phi_array[0]}")
axes[0].legend()


axes[1].set_xlabel("f")
axes[1].set_ylabel(f"ln(tau)/L^{power}")
axes[1].set_title(f"φ = {phi_array[1]}")
axes[1].legend()

fig.suptitle(f"P_Hermicity = {gamma_array[0]}")

plt.tight_layout()
plt.show()
    



#print(phase_diagram)