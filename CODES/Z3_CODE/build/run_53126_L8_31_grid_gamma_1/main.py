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





N = 31  # GRID DIMENSION
x_vals = np.linspace(0, 1, N)  # x from 0 to 1.1
y_vals = np.linspace(0, 1, N)  # y from 0 to 1.1

# Initialize matrix: shape (11, 11, 2) → each cell holds [r, phi]
grid = np.zeros((N, N, 2))

for i in range(N):       # row index → y-axis (bottom to top)
    for j in range(N):   # col index → x-axis (left to right)
        x = x_vals[j]
        y = y_vals[i]
        r   = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)   # in [0, pi/2] since x,y ≥ 0
        grid[i, j] = [r, phi]



#param_array=[[0.1, 0], [0.1, np.pi/6]]
gamma_array = [1]
L_array = [8]

index = int(sys.argv[1])

for i in range(index*4, index*4+4):

    r_phi_values = grid[int(i / 31), int(i) % 31]
    size_autocorrelator_SIMULATION_NPY(-10, 25, .5, [8], r_phi_values, gamma_array[0])



