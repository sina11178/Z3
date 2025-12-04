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





log_times = np.arange(0, 30, .5)

param_array=[[0.1, 0], [0.1, np.pi/6]]
gamma_array = [0, .2, .4, .6, .8, 1]
L_array = [4, 5, 6, 7, 8]

#size_autocorrelator_SIMULATION_NPY(log_times[(int(sys.argv[1]) % 10)*6], log_times[(int(sys.argv[1]) % 10)*6 + 5], .5, [L_array[0]], param_array[int(int(sys.argv[1])/60)], gamma_array[int(int(sys.argv[1])/10)%6])
LACZOS_COEF_vary_gam_SIM(100, L_array[int(sys.argv[1]) % 5], gamma_array[int(int(sys.argv[1])/5)%6], param_array[int(int(sys.argv[1])/30)])


