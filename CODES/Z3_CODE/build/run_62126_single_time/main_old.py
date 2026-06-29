'''

WE WILL CALL TO RUN RESPECTIVE SIMULATION HERE

'''

from f_phi_autocorr_SIM import f_phi_autocorrelator_SIMULATION_EXCEL, f_phi_autocorrelator_SIMULATION_NPY
from size_comparison_autocorr_SIM import size_autocorrelator_SIMULATION_EXCEL, size_autocorrelator_SIMULATION_NPY
from laczos_coeff_sim import LACZOS_COEF_vary_gam_SIM
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
import pandas as pd
import sys




#f_array = np.linspace(0, 1, 21)
f_array = [0.2]

#phi_array = [0, np.pi/6]
phi_array = [np.pi/6]
#param_array=[[0.1, 0], [0.1, np.pi/6]]
#gamma_array = [0, 0.5, 1]
gamma_array = [1]
#L_array = [4, 5, 6, 7, 8]
L_array = [10]

interval = np.log(7.6977)
L9_log_time = np.log(1644027) + interval
L10_log_time = L9_log_time + interval
L11_log_time = L10_log_time + interval

times = [L10_log_time, L11_log_time]
print("Running")

index = int(sys.argv[1])
#for i in range(len(L_array)):
    #size_autocorrelator_SIMULATION_NPY(-9, 27, .5, [L_array[i]], [f_array[int(index)%21], phi_array[int(int(index)/21)%2]], gamma_array[int(int(index)/42)])
size_autocorrelator_SIMULATION_NPY(times[0], times[0]+1, 1, [L_array[0]], [f_array[index], phi_array[0]], gamma_array[0])


