import numpy as np
import os
import pandas as pd
import shutil


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
gamma_array = [0]
L_array = [8]
total_points = N * N  # Total number of (r, phi) pairs

#param_array=[[0.1, 0], [0.1, np.pi/6]]
#gamma_array = [0, .1, .5, .6, .9, 1]
#L_array = [4, 5, 6, 7, 8]


for i in range(241):
     index = i
     start = index * 4
     end = min(start + 4, total_points)
     for j in range(start, end):
        params = grid[int(int(j / 31)), int(j) % 31]
        if os.path.isfile('./b'+str(index)+'/autocorr_gamma_' + str( gamma_array[0]) + '_L' + str(L_array[0]) + '_f_' + str(round(params[0],3)) + '_phi_' + str(round(params[1],3)) + '_NPY.npy') == True:

            #p = np.load('./b'+str(i)+'/size_autocorr_comparison_gamma_' + str(gamma_list[int(i/4)]) + '_NPY.npy')
            #data_list.append(p[1, :])
            f = params[0]
            phi = round(params[1], 3)

            destination = '../PHASE_DIAGRAM_DATA/L'+str(L_array[0])+'/gamma_'+str(gamma_array[0])+'/'
            f_n = './b'+str(index)+'/autocorr_gamma_' + str(gamma_array[0]) + '_L' + str(L_array[0]) + '_f_' + str(round(params[0],3)) + '_phi_' + str(round(params[1],3)) + '_NPY.npy'
            os.makedirs(destination, exist_ok=True)  # <-- add this
            shutil.copy(f_n, destination)
            print("copied")
        else:
            print(f"File not found: {index}")


print("successful")