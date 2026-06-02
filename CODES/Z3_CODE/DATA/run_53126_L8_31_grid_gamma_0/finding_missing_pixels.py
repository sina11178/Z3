import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from matplotlib.colors import LinearSegmentedColormap


cmap = LinearSegmentedColormap.from_list(
"cmap",
[
"#071F59", 
"#C4CFDE", 
"#FFFFFF", 
"#FFFFFF", 
"#FFFFFF",
"#FFFFFF" 
]
)




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


gamma_array = [0]
L_array = [8]
missing_pixels = []
time_list = np.exp(np.arange(-10, 25, 0.5))

phase_diagram = np.zeros((N, N))  # To store decay rates for each (r, phi) pair
total_points = N * N  # Total number of (r, phi) pairs

for i in range(241):
     index = i
     start = index * 4
     end = min(start + 4, total_points)
     for j in range(start, end):
        params = grid[int(int(j / 31)), int(j) % 31]
        if os.path.isfile('./b'+str(index)+'/autocorr_gamma_' + str( gamma_array[0]) + '_L' + str(L_array[0]) + '_f_' + str(round(params[0],3)) + '_phi_' + str(round(params[1],3)) + '_NPY.npy') == True:
            continue
        else:
            missing_pixels.append(params)


np.save("missing_pixels_L8_gamma_0.npy", np.array(missing_pixels))