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
L_array = [4, 5, 6, 7, 8]

time_list = np.exp(np.arange(-9, 27, 0.5))


def exp_model_new(t, c, a_1, g_1, a_2, g_2):
    term_1 = a_1 * np.exp(-t*g_1)
    term_2 = a_2 * np.exp(-t*g_2)

    return c + term_1 + term_2



fig, axes = plt.subplots(1, 2, figsize=(8, 6))


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
                    exp_model_new,
                    time_list,
                    data,
                    #p0= [data[-1], 0, 1, (1-data[-1]), 1 / time_list[len(time_list)//2]],#[data_current[-1], data_current[0] - data_current[-1], 1 / time_list_new[len(time_list_new)//2]],
                    bounds = ([0, 0, 0, 0, 0], [.5, 1, 10, 1, 10]),
                    maxfev=10000
                )
                fit_curve = exp_model_new(time_list, *popt)
                c, a_1, g_1, a_2, g_2,  = popt
            except RuntimeError:
                print(f"Fit did not converge for L={l}")

            print(g_1)
            print(g_2)
            print(a_1)
            print(a_2)
            threshold = 0.5
            if np.abs(a_2) > threshold and np.abs(a_1) > threshold:
                if g_1 >= g_2:
                    tau = 1/g_2
                else:
                    tau = 1/g_1
            elif np.abs(a_2) > threshold and np.abs(a_1) < threshold:
                tau = 1/g_2
            elif np.abs(a_2) < threshold and np.abs(a_1) > threshold:
                tau = 1/g_1
            else:
                if g_1 >= g_2:
                    tau = 1/g_2
                else:
                    tau = 1/g_1
            
            #tau = (a_1/g_1) + (a_2/g_2)
            #tau = tau / (a_1+a_2)
            f_1 = a_1/(a_1+a_2)
            f_2 = a_2/(a_1+a_2)
            thresh = 0.1

            if f_1< thresh:
                tau = 1/g_2
            elif f_2 < thresh:
                tau =1/g_1
            else:
                if g_1 <= g_2:
                    tau = 1/g_1
                else:
                    tau = 1/g_2


            '''
            if int(index/21)%2 == 1 and index%21 >= 3:
                data_new = data[19:]
                time_list_new = time_list[19:]


                target = data_new[0]/np.e

                # Find index just after crossing
                idx = np.argmax(data_new < target)

                # Two points straddling the threshold
                t0, t1 = time_list_new[idx-1], time_list_new[idx]
                d0, d1 = data_new[idx-1], data_new[idx]

                # Linear interpolation to find exact crossing time
                if d1 == d0:
                    tau = t0 
                else:
                    tau = (t0) + (target - d0) * (t1 - t0) / (d1 - d0)

            else:
                target = data[0]/np.e

                # Find index just after crossing
                idx = np.argmax(data < target)

                # Two points straddling the threshold
                t0, t1 = time_list[idx-1], time_list[idx]
                d0, d1 = data[idx-1], data[idx]

                # Linear interpolation to find exact crossing time
                if d1 == d0:
                    tau = t0
                else:
                    tau = t0 + (target - d0) * (t1 - t0) / (d1 - d0)
            
            '''


            if int(index/21)%2 == 0:
                scaled_lifetime_phi_0.append(np.log(tau)/np.sqrt(l))
            elif int(index/21)%2 == 1:
                scaled_lifetime_phi_pi_6.append(np.log(tau)/np.sqrt(l))
            else:
                raise ValueError("Houston, we have a problem")


        else:
            print(f"File not found for pixel {j}")


    scaled_lifetime_phi_0[0] = 27
    scaled_lifetime_phi_pi_6[0] = 27
    axes[0].plot(f_array[1:], scaled_lifetime_phi_0[1:], label=f"L={l}")
    axes[1].plot(f_array[1:], scaled_lifetime_phi_pi_6[1:], label=f"L={l}")


axes[0].set_xlabel("f")
axes[0].set_ylabel("ln(tau)/√L")
axes[0].set_title(f"φ = {phi_array[0]}")
axes[0].legend()


axes[1].set_xlabel("f")
axes[1].set_ylabel("ln(tau)/√L")
axes[1].set_title(f"φ = {phi_array[1]}")
axes[1].legend()

fig.suptitle(f"P_Hermicity = {gamma_array[0]}, threshold = {thresh}")

plt.tight_layout()
plt.show()
    



#print(phase_diagram)