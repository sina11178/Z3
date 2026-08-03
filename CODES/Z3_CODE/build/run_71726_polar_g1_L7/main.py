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





Nr = 40
N_phi = 40

r_vals = np.linspace(0, 1, Nr)
phi_vals = np.linspace(0, 2*np.pi/3, N_phi)

#R, Phi = np.meshgrid(r_vals, phi_vals, indexing='ij')



#param_array=[[0.1, 0], [0.1, np.pi/6]]
gamma_array = [0]
L_array = [5]

index = int(sys.argv[1])
ncores = 128

N_pixels = Nr * N_phi

start = index * N_pixels // ncores
end   = (index + 1) * N_pixels // ncores


for k in range(start, end):

    i = k // N_phi
    j = k % N_phi

    r_phi_values = [r_vals[i], phi_vals[j]]
    
    size_autocorrelator_SIMULATION_NPY(-10, 40, .5, [L_array[0]], r_phi_values, gamma_array[0])



