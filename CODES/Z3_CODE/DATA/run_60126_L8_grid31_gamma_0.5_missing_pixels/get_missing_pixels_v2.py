import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap

missing_pixels = np.load("missing_pixels_L8_gamma_0.5.npy")



N = 31  # GRID DIMENSION
x_vals = np.linspace(0, 1, N)  # x from 0 to 1.1
y_vals = np.linspace(0, 1, N)  # y from 0 to 1.1

# Initialize matrix: shape (31, 31, 2)~F~R each cell holds [r, phi]
grid = np.zeros((N, N, 2))

for i in range(N):       # row index �~F~R y-axis (bottom to top)
    for j in range(N):   # col index �~F~R x-axis (left to right)
        x = x_vals[j]
        y = y_vals[i]
        r   = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)   # in [0, pi/2] since x,y �~I� 0
        grid[i, j] = [r, phi]


gamma_array = [0.5]
L_array = [8]
missing_pixels_v2 = []

phase_diagram = np.zeros((N, N))  # To store decay rates for each (r, phi) pair
#total_points = N * N  # Total number of (r, phi) pairs

for i in range(222):
     index = i
     start = index * 2
     end = min(start + 2, len(missing_pixels))
     for j in range(start, end):
        params = missing_pixels[j]
        if os.path.isfile('./b'+str(index)+'/autocorr_gamma_' + str( gamma_array[0]) + '_L' + str(L_array[0]) + '_f_' + str(round(params[0],3)) + '_phi_' + str(round(params[1],3)) + '_NPY.npy') == True:
            continue
        else:
            missing_pixels_v2.append(params)

print(len(missing_pixels_v2))
np.save("missing_pixels_L8_gamma_0.5_V2.npy", np.array(missing_pixels_v2))