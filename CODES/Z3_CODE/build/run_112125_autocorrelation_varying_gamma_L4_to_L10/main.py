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







param_array=[[0.1, 0], [0.1, np.pi/6]]
gamma_array = [0, .2, .4, .6, .8, 1]
L_array = [4, 6, 8, 9, 10]

size_autocorrelator_SIMULATION_NPY(0, 30, .5, [L_array[int(sys.argv[1]) % 5]], param_array[int(sys.argv[1]%30)], gamma_array[int(int(sys.argv[1])/5)%6])



