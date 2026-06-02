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


gamma_array = [1]
L_array = [4]

time_list = np.exp(np.arange(-10, 25, 0.5))

phase_diagram = np.zeros((N, N))  # To store decay rates for each (r, phi) pair
total_points = N * N  # Total number of (r, phi) pairs

for i in range(121):
     index = i
     start = index * 8
     end = min(start + 8, total_points)
     for j in range(start, end):
        params = grid[int(int(j / 31)), int(j) % 31]
        if os.path.isfile('./b'+str(index)+'/autocorr_gamma_' + str( gamma_array[0]) + '_L' + str(L_array[0]) + '_f_' + str(round(params[0],3)) + '_phi_' + str(round(params[1],3)) + '_NPY.npy') == True:
            data = np.load('./b'+str(index)+'/autocorr_gamma_' + str(gamma_array[0]) + '_L' + str(L_array[0]) + '_f_' + str(round(params[0],3)) + '_phi_' + str(round(params[1],3)) + '_NPY.npy')[1]
            if data[0] <0.01:
                phase_diagram[int(j / 31), int(j % 31)] = None
                continue
            #target = data[0] / np.e
            target = 1/np.e

            # Find index just after crossing
            idx = np.argmax(data < target)

            # Two points straddling the threshold
            t0, t1 = time_list[idx-1], time_list[idx]
            d0, d1 = data[idx-1], data[idx]

            # Linear interpolation to find exact crossing time
            tau = t0 + (target - d0) * (t1 - t0) / (d1 - d0)
            decay_rate = 1 / tau
            phase_diagram[int(j / 31), int(j % 31)] = decay_rate

            '''
            log_data = np.log(data)
            #print(data)
            slope, intercept = np.polyfit(time_list, log_data, 1)

            decay_rate = -slope

            phase_diagram[int(i / 11), int(i % 11)] = decay_rate
            '''
            if i == 2:
                plt.plot(time_list, data, label='Autocorrelation')
                plt.xscale('log')

        
        else:
            print(f"File not found for index {i}")


#print(phase_diagram)

fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    phase_diagram,
    origin='lower',                  # (0,0) at bottom-left
    extent=[0, 1, 0, 1],        # x and y axes from 0 to 1.1
    aspect='auto',
    cmap=cmap,
    vmin=0, vmax=0.7
)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Decay Rate λ')

ax.set_xlabel('$f/J * cos(φ)$')
ax.set_ylabel('$f/J * sin(φ)$')
ax.set_title(f'Phase Diagram: L = {L_array[0]}, Γ = {gamma_array[0]}, Polar Grid (f/J, φ)')
plt.rcParams["figure.figsize"] = (18, 7)
#ax.clim(0, 1.5)  # Set color limits for better contrast
plt.savefig('phase_diagram.png', dpi=600)
plt.savefig('phase_diagram.pdf', format='pdf')

plt.tight_layout()
plt.show()